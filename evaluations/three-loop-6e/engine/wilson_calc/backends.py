"""Fresh, bounded HyperInt evaluation of centered single-opposite moments.

The public operation accepts four integer endpoint weights and an epsilon
order, never arbitrary Maple source.  It computes coefficients of

 H_ab;cd(e,z) = -(1+z) z**(-e/2) integral_[0,1]^2
   D_z(t,s)**(-1+e) t**(a*e) (1-t)**(b*e)
   s**(c*e) (1-s)**(d*e) dt ds.

No graph/result table or historical moment cache is consulted.  Root-to-moment
and physical normalization bridges belong to exchange/assembly, not here.
"""
from __future__ import annotations

import ast
import ctypes
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from typing import Iterable
from uuid import uuid4

import sympy as sp

from .algebra import from_sympy, rational
from .cas_budget import CASBudgetExhausted


class BackendOutputError(ValueError):
    """The CAS result is outside this adapter's explicitly supported syntax."""


class ProcessLimitUnavailable(RuntimeError):
    """This host cannot enforce the requested owned-process-tree limits."""


def run_bounded_process(command, *, cwd, stdout_path, timeout_seconds, memory_mb) -> dict:
    """Run an owned process tree in a Windows kill-on-close, memory-limited Job.

    A Python launcher waits for a parent message before spawning the command.
    It is attached to the Job before that message is sent, so all descendants
    inherit the Job; a fast child cannot escape assignment. The memory limit
    is total committed memory of the Job, enforced by Windows, not sampled RSS.
    Non-Windows hosts explicitly refuse until an equivalent adapter is added.
    """
    if os.name != "nt":
        raise ProcessLimitUnavailable("aggregate process memory limits require the Windows Job adapter")
    if (not isinstance(command, (list, tuple)) or not command
            or any(not isinstance(v, str) or "\0" in v for v in command)):
        raise ValueError("process command must be a nonempty string argument list")
    if (isinstance(timeout_seconds, bool) or not math.isfinite(timeout_seconds)
            or timeout_seconds <= 0):
        raise ValueError("positive finite process timeout required")
    if type(memory_mb) is not int or memory_mb <= 0:
        raise ValueError("positive integer memory_mb required")
    from ctypes import wintypes as w

    class BasicLimits(ctypes.Structure):
        _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong),
                    ("PerJobUserTimeLimit", ctypes.c_longlong), ("LimitFlags", w.DWORD),
                    ("MinimumWorkingSetSize", ctypes.c_size_t), ("MaximumWorkingSetSize", ctypes.c_size_t),
                    ("ActiveProcessLimit", w.DWORD), ("Affinity", ctypes.c_size_t),
                    ("PriorityClass", w.DWORD), ("SchedulingClass", w.DWORD)]

    class IOCounters(ctypes.Structure):
        _fields_ = [(name, ctypes.c_ulonglong) for name in
                    ("ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
                     "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]

    class ExtendedLimits(ctypes.Structure):
        _fields_ = [("BasicLimitInformation", BasicLimits), ("IoInfo", IOCounters),
                    ("ProcessMemoryLimit", ctypes.c_size_t), ("JobMemoryLimit", ctypes.c_size_t),
                    ("PeakProcessMemoryUsed", ctypes.c_size_t), ("PeakJobMemoryUsed", ctypes.c_size_t)]

    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.CreateJobObjectW.argtypes = [ctypes.c_void_p, w.LPCWSTR]
    kernel.CreateJobObjectW.restype = w.HANDLE
    kernel.SetInformationJobObject.argtypes = [w.HANDLE, ctypes.c_int, ctypes.c_void_p, w.DWORD]
    kernel.SetInformationJobObject.restype = w.BOOL
    kernel.AssignProcessToJobObject.argtypes = [w.HANDLE, w.HANDLE]
    kernel.AssignProcessToJobObject.restype = w.BOOL
    kernel.QueryInformationJobObject.argtypes = [w.HANDLE, ctypes.c_int, ctypes.c_void_p, w.DWORD, ctypes.c_void_p]
    kernel.QueryInformationJobObject.restype = w.BOOL
    kernel.TerminateJobObject.argtypes = [w.HANDLE, w.UINT]
    kernel.TerminateJobObject.restype = w.BOOL
    kernel.CloseHandle.argtypes = [w.HANDLE]
    kernel.CloseHandle.restype = w.BOOL
    job = kernel.CreateJobObjectW(None, None)
    if not job:
        raise ProcessLimitUnavailable("Windows could not create an owned Job")
    limits = ExtendedLimits()
    limits.BasicLimitInformation.LimitFlags = 0x2000 | 0x200  # KILL_ON_JOB_CLOSE | JOB_MEMORY
    limits.JobMemoryLimit = memory_mb * 1024 * 1024
    if not kernel.SetInformationJobObject(job, 9, ctypes.byref(limits), ctypes.sizeof(limits)):
        kernel.CloseHandle(job)
        raise ProcessLimitUnavailable("Windows could not set the Job memory limit")
    launcher = ("import json,subprocess,sys; "
                "p=json.loads(sys.stdin.readline()); "
                "sys.exit(subprocess.call(p['command'],cwd=p['cwd'],stdin=subprocess.DEVNULL))")
    child = None
    started = time.monotonic()
    stop_reason = None
    peak = 0
    try:
        with Path(stdout_path).open("wb") as stream:
            child = subprocess.Popen([sys.executable, "-u", "-c", launcher], cwd=cwd,
                stdin=subprocess.PIPE, stdout=stream, stderr=subprocess.STDOUT,
                creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0))
            if not kernel.AssignProcessToJobObject(job, w.HANDLE(int(child._handle))):
                # The waiting launcher has not received a command and has no descendants.
                child.kill()
                child.wait()
                raise ProcessLimitUnavailable("Windows refused Job assignment; command was not started")
            message = json.dumps({"command": list(command), "cwd": str(Path(cwd).resolve())}) + "\n"
            try:
                child.communicate(message.encode("utf-8"), timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                stop_reason = "time"
                kernel.TerminateJobObject(job, 1)
                child.communicate()
            observed = ExtendedLimits()
            if kernel.QueryInformationJobObject(job, 9, ctypes.byref(observed), ctypes.sizeof(observed), None):
                peak = int(observed.PeakJobMemoryUsed)
    finally:
        # Also kills descendants left behind by a successful/failed leader.
        kernel.CloseHandle(job)
        if child is not None and child.poll() is None:
            child.wait()
    if child.returncode != 0 and stop_reason is None:
        tail = Path(stdout_path).read_bytes()[-16000:].decode("utf-8", errors="replace")
        if re.search(r"MemoryError|out of memory|unable to allocate|memory allocation|memory limit", tail, re.I):
            stop_reason = "memory"
    return {"returncode": child.returncode, "stop_reason": stop_reason,
            "wall_seconds": time.monotonic() - started, "peak_job_memory_bytes": peak,
            "memory_mb": memory_mb, "memory_enforcement": "windows-job-total-committed",
            "owned_process_tree_cleanup": "kill-on-job-close"}


@dataclass(frozen=True)
class MapleConfig:
    executable: Path
    hyperint: Path
    periods: Path
    timeout_seconds: float = 120.0
    max_words: int = 50_000
    max_output_bytes: int = 8_000_000
    memory_mb: int = 1024

    def validate(self) -> None:
        for path in (self.executable, self.hyperint, self.periods):
            if not Path(path).is_file():
                raise ValueError(f"backend file is unavailable: {path}")
            _maple_path(path)
        if (isinstance(self.timeout_seconds, bool)
                or not math.isfinite(self.timeout_seconds)
                or not 0 < self.timeout_seconds <= 120):
            raise ValueError("each fresh Maple process needs a timeout in (0, 120]")
        if type(self.max_words) is not int or self.max_words <= 0:
            raise ValueError("positive integer max_words required")
        if type(self.max_output_bytes) is not int or self.max_output_bytes <= 0:
            raise ValueError("positive integer max_output_bytes required")
        if type(self.memory_mb) is not int or self.memory_mb <= 0:
            raise ValueError("positive integer memory_mb required")

    def fingerprint(self) -> dict:
        self.validate()
        return {name: {"path": str(Path(path).resolve()), "sha256": _sha(path)}
                for name, path in (("executable", self.executable),
                                   ("hyperint", self.hyperint),
                                   ("period_table", self.periods))}


def find_maple_config(**limits) -> MapleConfig | None:
    """Find the configured local backend; discovery never starts a process."""
    executable = os.environ.get("WILSON_MAPLE") or shutil.which("cmaple")
    if not executable:
        candidate = Path("C:/Program Files/Maple 2025/bin.X86_64_WINDOWS/cmaple.exe")
        if candidate.is_file():
            executable = str(candidate)
    library = Path(os.environ.get("WILSON_HYPERINT",
        "C:/Users/Lenovo/Documents/learning/work/hyperint/HyperInt.mpl"))
    periods = Path(os.environ.get("WILSON_HYPERINT_PERIODS",
        str(library.with_name("periodLookups.m"))))
    if not executable or not library.is_file() or not periods.is_file():
        return None
    config = MapleConfig(Path(executable), library, periods, **limits)
    config.validate()
    return config


def _sha(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _maple_path(path) -> str:
    value = Path(path).resolve().as_posix()
    if any(c in value for c in ('"', '\n', '\r', '\0')):
        raise ValueError("unsupported character in a Maple file path")
    return value


def _validate_parameters(parameters, order) -> tuple[int, int, int, int]:
    if (not isinstance(parameters, (list, tuple)) or len(parameters) != 4
            or any(type(p) is not int or p < 0 for p in parameters)):
        raise ValueError("four nonnegative integer endpoint weights required")
    if type(order) is not int or order < 0:
        raise ValueError("nonnegative integer epsilon order required")
    return tuple(parameters)


_HLOG = re.compile(r"Hlog\(z,\s*\[([^\]]*)\]\)")
_ZETA = re.compile(r"zeta\[(\d+)\]")


def parse_hyperint_coefficient(text: str, *, max_words=50_000, max_nodes=200_000) -> dict:
    """Read a linear Hlog(z, {0,-1}*) sum without evaluating source text.

    Single zeta periods of any integer weight >= 2 are retained in the existing
    scalar AST.  A new alphabet or an unrecognized period is a stated adapter
    limitation, never silently discarded.  The positive-z single-opposite
    alphabet is the scope of this backend, not a restriction on HLogSum.
    """
    from .hlog import GPL, HLogSum, HLogTerm

    if not isinstance(text, str) or not text.strip():
        raise BackendOutputError("empty HyperInt coefficient")
    if len(text) > 8_000_000:
        raise BackendOutputError("HyperInt coefficient exceeds parser byte budget")
    if type(max_words) is not int or max_words <= 0:
        raise ValueError("positive max_words required")
    if type(max_nodes) is not int or max_nodes <= 0:
        raise ValueError("positive max_nodes required")
    words: list[tuple[int, ...]] = []
    word_indices: dict[tuple[int, ...], int] = {}
    names: dict[str, sp.Expr] = {"z": sp.Symbol("z")}

    def replace_hlog(match):
        raw = match.group(1).strip()
        entries = raw.split(",") if raw else []
        if any(not re.fullmatch(r"\s*(?:0|-1)\s*", x) for x in entries):
            raise BackendOutputError("single-opposite output has an unsupported alphabet")
        word = tuple(int(x) for x in entries)
        if len(word) > 64:
            raise BackendOutputError("HyperInt word exceeds parser length budget")
        if not word:
            return "1"
        if word not in word_indices:
            word_indices[word] = len(words)
            words.append(word)
            if len(words) > max_words:
                raise BackendOutputError("coefficient exceeds the HLOG word budget")
        name = "hlog_" + str(word_indices[word])
        names[name] = sp.Symbol(name)
        return name

    source = _HLOG.sub(replace_hlog, text.strip())
    if len(words) > max_words:
        raise BackendOutputError("coefficient exceeds the HLOG word budget")

    def replace_zeta(match):
        if len(match.group(1)) > 3:
            raise BackendOutputError("zeta index exceeds parser budget")
        n = int(match.group(1))
        if not 2 <= n <= 64:
            raise BackendOutputError("divergent or unsupported zeta period")
        name = "zeta_" + str(n)
        names[name] = sp.zeta(n)
        return name

    source = _ZETA.sub(replace_zeta, source).replace("^", "**")
    try:
        tree = ast.parse(source, mode="eval")
    except (SyntaxError, MemoryError, RecursionError) as exc:
        raise BackendOutputError("unsupported HyperInt output syntax") from exc
    nodes = list(ast.walk(tree))
    if len(nodes) > max_nodes:
        raise BackendOutputError("HyperInt coefficient exceeds parser node budget")
    if any(isinstance(node, ast.Constant) and type(node.value) is int
           and node.value.bit_length() > 4096 for node in nodes):
        raise BackendOutputError("HyperInt coefficient integer exceeds parser budget")
    stack = [(tree.body, 0)]
    while stack:
        node, depth = stack.pop()
        if depth > 256:
            raise BackendOutputError("HyperInt coefficient exceeds parser depth budget")
        stack.extend((child, depth + 1) for child in ast.iter_child_nodes(node))

    # Reject nonlinear HLOGs and explosive distributive products before SymPy
    # constructs or expands any expression. This intentionally accepts only
    # the linear fibration-basis export contract, not arbitrary identities that
    # might become linear after cancellation.
    expansion_budget = min(max_nodes, 4 * max_words)
    numeric_bit_budget = 65_536

    def integer_exponent(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return node.value
        if (isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd))
                and isinstance(node.operand, ast.Constant) and type(node.operand.value) is int):
            return -node.operand.value if isinstance(node.op, ast.USub) else node.operand.value
        raise BackendOutputError("CAS exponent must be an integer literal")

    def preflight(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return 0, 1, max(1, node.value.bit_length())
        if isinstance(node, ast.Name) and node.id in names:
            # Even zeta values include modest rational Bernoulli coefficients.
            return int(node.id.startswith("hlog_")), 1, 512 if node.id.startswith("zeta_") else 1
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
            return preflight(node.operand)
        if isinstance(node, ast.BinOp):
            left_degree, left_terms, left_bits = preflight(node.left)
            if isinstance(node.op, ast.Pow):
                exponent = integer_exponent(node.right)
                if abs(exponent) > 64:
                    raise BackendOutputError("unsupported CAS power")
                if left_degree and (exponent not in (0, 1) or not isinstance(node.left, ast.Name)):
                    raise BackendOutputError("CAS output was not reduced to a linear HLOG sum")
                degree = left_degree if exponent else 0
                count = left_terms ** abs(exponent)
                bits = left_bits * abs(exponent) + 1
            else:
                right_degree, right_terms, right_bits = preflight(node.right)
                if isinstance(node.op, (ast.Add, ast.Sub)):
                    degree, count = max(left_degree, right_degree), left_terms + right_terms
                    bits = max(left_bits, right_bits) + 1
                elif isinstance(node.op, ast.Mult):
                    degree, count = left_degree + right_degree, left_terms * right_terms
                    bits = left_bits + right_bits
                elif isinstance(node.op, ast.Div):
                    if right_degree:
                        raise BackendOutputError("nonpolynomial HLOG output")
                    degree, count = left_degree, left_terms
                    bits = left_bits + right_bits
                else:
                    raise BackendOutputError("unsupported CAS operation")
            if degree > 1:
                raise BackendOutputError("CAS output was not reduced to a linear HLOG sum")
            if count > expansion_budget:
                raise BackendOutputError("CAS expression exceeds distributive expansion budget")
            if bits > numeric_bit_budget:
                raise BackendOutputError("CAS expression exceeds integer-growth budget")
            return degree, count, bits
        raise BackendOutputError("unsupported executable or mathematical CAS output")

    preflight(tree.body)

    def read(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return sp.Integer(node.value)
        if isinstance(node, ast.Name) and node.id in names:
            return names[node.id]
        if isinstance(node, ast.UnaryOp):
            if isinstance(node.op, ast.USub):
                return -read(node.operand)
            if isinstance(node.op, ast.UAdd):
                return read(node.operand)
        if isinstance(node, ast.BinOp):
            left, right = read(node.left), read(node.right)
            if isinstance(node.op, ast.Add): return left + right
            if isinstance(node.op, ast.Sub): return left - right
            if isinstance(node.op, ast.Mult): return left * right
            if isinstance(node.op, ast.Div):
                if right == 0:
                    raise BackendOutputError("zero denominator in CAS output")
                return left / right
            if isinstance(node.op, ast.Pow) and right.is_Integer and abs(int(right)) <= 64:
                if left == 0 and right < 0:
                    raise BackendOutputError("zero denominator in CAS output")
                return left ** right
        raise BackendOutputError("unsupported executable or mathematical CAS output")

    expression = sp.expand(read(tree.body))
    symbols = [names["hlog_" + str(i)] for i in range(len(words))]
    try:
        if symbols and sp.Poly(expression, *symbols).total_degree() > 1:
            raise BackendOutputError("CAS output was not reduced to a linear HLOG sum")
    except sp.PolynomialError as exc:
        raise BackendOutputError("nonpolynomial HLOG output") from exc
    endpoint = {"kind": "var", "name": "z"}
    terms = []
    constant = expression.subs(dict.fromkeys(symbols, 0))
    for word, coefficient in [((), constant), *[(w, expression.coeff(g))
                                              for w, g in zip(words, symbols)]]:
        if coefficient == 0:
            continue
        try:
            scalar = from_sympy(coefficient)
        except ValueError as exc:
            raise BackendOutputError("unsupported scalar in HLOG output") from exc
        terms.append(HLogTerm(scalar, GPL([rational(x) for x in word], endpoint,
                                        branch="positive-real", endpoint_kind="ordinary")))
    if len(terms) > max_words:
        raise BackendOutputError("coefficient exceeds the HLOG term budget")
    return HLogSum(terms, parameters=("z",)).to_json()


def _script(parameters, order, config: MapleConfig) -> str:
    a, b, c, d = _validate_parameters(parameters, order)
    return f'''restart:
_hyper_autoload_periods := ["{_maple_path(config.periods)}"]:
read "{_maple_path(config.hyperint)}":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
_hyper_ignore_nonlinear_polynomials:=false:
try
  wbBase:=-(1+z)/((1+x)*(1+y)*(1+z*x*y)):
  wbLog:=ln(1+z*x*y)-ln(1+x)-ln(1+y)-ln(z)/2
    +({a})*(ln(x)-ln(1+x))-({b})*ln(1+x)
    +({c})*(ln(y)-ln(1+y))-({d})*ln(1+y):
  wbAnswer:=hyperInt(convert(wbBase*wbLog^({order})/({order})!,HlogRegInf),[x,y]):
  wbRaw:=fibrationBasis(wbAnswer,[z]):
  wbRef:=subs(zz=z,fibrationBasis(convert(subs(z=1/zz,wbRaw),HlogRegInf),[zz])):
  for wbChannel in ["raw","reflected"] do
    if wbChannel="raw" then wbValue:=wbRaw else wbValue:=wbRef fi:
    if type(wbValue,`+`) then wbCount:=nops(wbValue) else wbCount:=1 fi:
    if wbCount>{config.max_words} then error "HLOG word budget exceeded" fi:
    wbFile:=fopen(cat(wbChannel,".txt"),WRITE,TEXT):
    fprintf(wbFile,"%a\\n",wbValue):fclose(wbFile):
  od:
  printf("WILSON_MOMENT_COMPLETE order={order}\\n"):
catch:
  printf("WILSON_MOMENT_FAILED %a\\n",[lastexception]):
end try:
quit:
'''


def compute_single_opposite_moments(parameters, through: int, *, run_dir,
        config: MapleConfig | None = None, orders: Iterable[int] | None = None,
        total_seconds: float | None = None, cas_budget=None) -> dict:
    """Compute requested coefficients in fresh, independent bounded processes.

    ``orders`` supports resume by requesting only missing orders.  It never
    reads previous answers; the caller validates/reuses its own checkpoint.
    Successful lower orders survive a later timeout in this return value and
    in their individual artifact directories.  Each invocation creates a new
    directory and preserves source, raw output, log and machine result.
    """
    parameters = _validate_parameters(parameters, through)
    selected = tuple(range(through + 1)) if orders is None else tuple(orders)
    if (len(set(selected)) != len(selected)
            or any(type(q) is not int or not 0 <= q <= through for q in selected)):
        raise ValueError("orders must be unique nonnegative integers <= through")
    selected = tuple(sorted(selected))
    config = config or find_maple_config()
    if config is None:
        return {"status": "unsupported", "coefficients": {}, "complete_orders": [],
                "requested_orders": list(selected), "diagnostic": "Maple/HyperInt unavailable"}
    fingerprint = config.fingerprint()
    total_seconds = config.timeout_seconds if total_seconds is None else total_seconds
    if (isinstance(total_seconds, bool) or not math.isfinite(total_seconds)
            or total_seconds <= 0):
        raise ValueError("positive finite total_seconds required")
    deadline = time.monotonic() + total_seconds
    root = Path(run_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    batch = root / ("single-opposite-" + uuid4().hex)
    batch.mkdir(exist_ok=False)
    request = {"operation": "centered-single-opposite-moments", "parameters": list(parameters),
               "through": through, "orders": list(selected), "context": "z > 0; expansion at epsilon=0",
               "definition": "-(1+z)*z^(-eps/2)*integral D^(-1+eps)*t^(a eps)*(1-t)^(b eps)*s^(c eps)*(1-s)^(d eps)",
               "integration_order": ["x", "y"], "domain_after_substitution": "positive quadrant",
               "substitution": "t=x/(1+x), s=y/(1+y)",
               "backend": fingerprint, "timeout_seconds_per_order": config.timeout_seconds,
               "total_seconds": total_seconds, "max_words": config.max_words,
               "memory_mb": config.memory_mb,
               "formal_state": "not_run"}
    (batch / "request.json").write_text(json.dumps(request, indent=2) + "\n", encoding="utf-8")
    coefficients, records = {}, []
    status, diagnostic = "computed", None
    for q in selected:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            status, diagnostic = "resource_exhausted", "total backend time budget exhausted"
            break
        target = batch / f"order-{q}"
        target.mkdir(exist_ok=False)
        script = target / "compute.mpl"
        script.write_text(_script(parameters, q, config), encoding="utf-8", newline="\n")
        command = [str(Path(config.executable).resolve()), "-q", str(script)]
        record = {"order": q, "directory": str(target), "command": command,
                  "script_sha256": _sha(script), "started_utc": datetime.now(timezone.utc).isoformat(),
                  "timeout_seconds": min(config.timeout_seconds, remaining)}
        try:
            runner = run_bounded_process if cas_budget is None else (
                lambda *a, **kw: cas_budget.run(run_bounded_process, *a, **kw))
            execution = runner(command, cwd=target, stdout_path=target / "maple.log",
                timeout_seconds=record["timeout_seconds"], memory_mb=config.memory_mb)
            record.update(execution)
            log_text = (target / "maple.log").read_text(encoding="utf-8", errors="replace")
            marker = f"WILSON_MOMENT_COMPLETE order={q}"
            if execution["stop_reason"] in ("time", "memory"):
                status, diagnostic = "resource_exhausted", f"coefficient {q} exhausted {execution['stop_reason']} budget"
            elif (execution["returncode"] != 0 or marker not in log_text
                    or "WILSON_MOMENT_FAILED" in log_text or "Error," in log_text):
                if "HLOG word budget exceeded" in log_text or "alloc" in log_text.lower() or "MemoryError" in log_text:
                    status, diagnostic = "resource_exhausted", "CAS word or memory budget exhausted"
                else:
                    status, diagnostic = "failed", "CAS did not complete the requested coefficient"
                record["log_tail"] = log_text[-4000:]
            else:
                payload = {}
                for channel in ("raw", "reflected"):
                    path = target / (channel + ".txt")
                    if not path.is_file() or path.stat().st_size > config.max_output_bytes:
                        raise BackendOutputError("missing or oversized coefficient artifact")
                    payload[channel] = parse_hyperint_coefficient(path.read_text(encoding="utf-8"),
                                                                  max_words=config.max_words)
                coefficients[str(q)] = payload
        except CASBudgetExhausted as exc:
            status, diagnostic = "resource_exhausted", str(exc)
            record["reason_code"] = exc.code
        except ProcessLimitUnavailable as exc:
            status, diagnostic = "unsupported", str(exc)
        except BackendOutputError as exc:
            status, diagnostic = "unsupported", str(exc)
        except OSError as exc:
            status, diagnostic = "failed", str(exc)
        record["status"] = "computed" if str(q) in coefficients else status
        record["finished_utc"] = datetime.now(timezone.utc).isoformat()
        if str(q) in coefficients:
            record["coefficient"] = coefficients[str(q)]
        elif diagnostic:
            record["diagnostic"] = diagnostic
        (target / "result.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
        records.append(record)
        if record["status"] != "computed":
            break
    result = {"status": status, "coefficients": coefficients,
              "complete_orders": sorted(map(int, coefficients)), "requested_orders": list(selected),
              "provenance": request, "run_directory": str(batch), "records": records,
              "formal_state": "not_run"}
    if diagnostic:
        result["diagnostic"] = diagnostic
    if records and "reason_code" in records[-1]:
        result["reason_code"] = records[-1]["reason_code"]
    (batch / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result
