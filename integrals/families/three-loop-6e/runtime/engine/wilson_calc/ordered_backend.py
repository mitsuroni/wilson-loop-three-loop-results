"""Fresh finite-order integration of convergent ordered D-kernel products.

This adapter reads the actual prepared J and retains its nested finite bounds.
It does not use graph IDs, a Q/R answer dictionary, or an unchecked finite-part
prescription. It admits unweighted products at q=0, and two opposite kernels
with one eliminated adjacent pair through regular-kernel order two. Positive-z
ordinary convergence is a structural precondition; CAS divergence checks
remain enabled in every run.
"""
from __future__ import annotations

from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from uuid import uuid4

from preprocessing.ordered import arithmetic

from .backends import (BackendOutputError, MapleConfig, ProcessLimitUnavailable,
                       _maple_path, find_maple_config, parse_hyperint_coefficient,
                       run_bounded_process)
from .exchange import analyze_exchange
from .cas_budget import CASBudgetExhausted


def _hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _emit_scalar(node, name_map):
    """Emit Maple only from the existing closed scalar AST and fresh binders."""
    if node["kind"] == "rational":
        numerator, denominator = int(node["numerator"]), int(node["denominator"])
        if denominator <= 0:
            raise ValueError("positive rational denominator required")
        return f"({numerator}/{denominator})" if denominator != 1 else str(numerator)
    if node["kind"] == "var":
        if node["name"] not in name_map:
            raise ValueError("unbound variable in ordered backend")
        return name_map[node["name"]]
    args = [_emit_scalar(arg, name_map) for arg in node["args"]]
    operation = node["op"]
    if operation == "neg" and len(args) == 1:
        return f"(-({args[0]}))"
    signs = {"add": "+", "sub": "-", "mul": "*", "div": "/", "rpow": "^"}
    if operation not in signs or len(args) != 2:
        raise ValueError("unsupported arithmetic operation in ordered backend")
    return f"(({args[0]}){signs[operation]}({args[1]}))"


def prepare_ordered_constant(j, *, loop_order=None, endpoint_edges=None,
                             topology=None, variable_bindings=None):
    """Return checked emitted q=0 input, or a precise unsupported remainder."""
    analysis = analyze_exchange(j, loop_order=loop_order, endpoint_edges=endpoint_edges,
                                topology=topology, variable_bindings=variable_bindings)
    if analysis["status"] != "reduced":
        return {"status": "unsupported", "residuals": analysis["residuals"]}
    descriptor = analysis["descriptor"]
    n = descriptor["loop_order"]
    if descriptor["adjacent_count"] or descriptor["opposite_count"] != n or n not in (2, 3):
        return {"status": "unsupported", "residuals": [{"kind": "ordered_constant_rule_scope",
                "detail": "This q=0 rule requires two or three distinct opposite kernels and no adjacent densities.",
                "J": deepcopy(j)}]}
    name_map = {name: f"wox{index + 1}" for index, name in enumerate(j["variables"])}
    name_map.update(z1="z", eps="0")
    allowed = set(j["variables"]) | {"eps", "z1"}

    def emit(text):
        return _emit_scalar(arithmetic(text, allowed), name_map)

    integrand = emit(j["numerator"])
    for polynomial in j["propagator_polynomials"].values():
        integrand += f"/({emit(polynomial)})"
    inner_to_outer = list(reversed(j["integration_order_outer_to_inner"]))
    intervals = [f'{name_map[name]}={emit(j["bounds"][name][0])}..{emit(j["bounds"][name][1])}'
                 for name in inner_to_outer]
    convergence = {
        "kind": "ordinary_integrable_distinct_opposite_product",
        "context": "z1 > 0; epsilon in a fixed disk |eps| < 1/2",
        "argument": "Each endpoint belongs to exactly one opposite propagator. At either zero corner D or D_tilde is comparable to u+v for fixed positive z. The model (u+v)^(-1-eta) is integrable in two variables for eta<1. Enlarging the ordered domain to the cube yields a product of those finite bounds, uniformly on compact epsilon disks. The integral is holomorphic at zero, so its constant coefficient is the original ordinary ordered integral with eps=0.",
        "endpoint_interpretation": "ordinary limits of complete convergent integrals; no independent finite-part replacement",
        "source": "benchmarks/2ne/families/same-channel-three/method.md#2-有限性可换序与解析展开",
        "scope": "The same distinct-endpoint product bound applies to two kernels; this is not a claim about adjacent or vertex kernels.",
        "formal_state": "not_run",
    }
    return {"status": "ready", "operation": "ordered-opposite-constant",
            "quantity": "J", "order": 0, "epsilon_lower_bound": 0,
            "input_J": deepcopy(j), "input_J_sha256": descriptor["input_J_sha256"],
            "loop_order": n, "channel_partition": descriptor["channel_partition"],
            "original_outer_to_inner": list(j["integration_order_outer_to_inner"]),
            "execution_inner_to_outer": inner_to_outer,
            "variable_map": name_map, "maple_integrand": integrand,
            "maple_intervals": intervals, "ordinary_endpoint_admission": convergence,
            "formal_state": "not_run"}


def _script(prepared, config):
    intervals = ",".join(prepared["maple_intervals"])
    return f'''restart:
_hyper_autoload_periods := ["{_maple_path(config.periods)}"]:
read "{_maple_path(config.hyperint)}":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
_hyper_ignore_nonlinear_polynomials:=false:
try
  woInput:={prepared["maple_integrand"]}:
  woAnswer:=hyperInt(woInput,[{intervals}]):
  woRaw:=fibrationBasis(woAnswer,[z]):
  woRef:=subs(wozz=z,fibrationBasis(convert(subs(z=1/wozz,woRaw),HlogRegInf),[wozz])):
  for woChannel in ["raw","reflected"] do
    if woChannel="raw" then woValue:=woRaw else woValue:=woRef fi:
    if type(woValue,`+`) then woCount:=nops(woValue) else woCount:=1 fi:
    if woCount>{config.max_words} then error "HLOG word budget exceeded" fi:
    woFile:=fopen(cat(woChannel,".txt"),WRITE,TEXT):
    fprintf(woFile,"%a\\n",woValue):fclose(woFile):
  od:
  printf("WILSON_ORDERED_CONSTANT_COMPLETE\\n"):
catch:
  printf("WILSON_ORDERED_CONSTANT_FAILED %a\\n",[lastexception]):
end try:
quit:
'''


def compute_ordered_constant(j, *, run_dir, loop_order=None, endpoint_edges=None,
                             topology=None, variable_bindings=None,
                             config: MapleConfig | None = None, cas_budget=None):
    """Compute J_0 and its positive-z reflection in one fresh bounded process."""
    prepared = prepare_ordered_constant(j, loop_order=loop_order, endpoint_edges=endpoint_edges,
                                        topology=topology, variable_bindings=variable_bindings)
    if prepared["status"] != "ready":
        return {**prepared, "coefficients": {}, "complete_orders": [], "formal_state": "not_run"}
    config = config or find_maple_config()
    if config is None:
        return {"status": "unsupported", "coefficients": {}, "complete_orders": [],
                "diagnostic": "Maple/HyperInt unavailable", "provenance": prepared,
                "formal_state": "not_run"}
    fingerprint = config.fingerprint()
    root = Path(run_dir).resolve()
    root.mkdir(parents=True, exist_ok=True)
    target = root / ("ordered-constant-" + uuid4().hex)
    target.mkdir(exist_ok=False)
    request = {**prepared, "backend": fingerprint, "timeout_seconds": config.timeout_seconds,
               "memory_mb": config.memory_mb, "max_words": config.max_words}
    (target / "request.json").write_text(json.dumps(request, indent=2) + "\n", encoding="utf-8")
    script_path = target / "compute.mpl"
    script_path.write_text(_script(prepared, config), encoding="utf-8", newline="\n")
    command = [str(Path(config.executable).resolve()), "-q", str(script_path)]
    record = {"directory": str(target), "command": command, "script_sha256": _hash(script_path),
              "started_utc": datetime.now(timezone.utc).isoformat()}
    status, diagnostic, coefficients = "failed", None, {}
    try:
        runner = run_bounded_process if cas_budget is None else (
            lambda *a, **kw: cas_budget.run(run_bounded_process, *a, **kw))
        execution = runner(command, cwd=target, stdout_path=target / "maple.log",
                           timeout_seconds=config.timeout_seconds, memory_mb=config.memory_mb)
        record.update(execution)
        log_text = (target / "maple.log").read_text(encoding="utf-8", errors="replace")
        if execution["stop_reason"] in ("time", "memory"):
            status, diagnostic = "resource_exhausted", "ordered constant exceeded process " + execution["stop_reason"] + " budget"
        elif (execution["returncode"] != 0 or "WILSON_ORDERED_CONSTANT_COMPLETE" not in log_text
              or "WILSON_ORDERED_CONSTANT_FAILED" in log_text or "Error," in log_text):
            if "HLOG word budget exceeded" in log_text or "alloc" in log_text.lower() or "MemoryError" in log_text:
                status, diagnostic = "resource_exhausted", "ordered constant exceeded word or memory budget"
            elif "divergen" in log_text.lower() or "FAILED: factors" in log_text:
                status, diagnostic = "unsupported", "the current finite-interval route could not establish ordinary endpoints or linear factorization"
            else:
                status, diagnostic = "failed", "CAS did not complete the original ordered integral"
            record["log_tail"] = log_text[-4000:]
        else:
            payload = {}
            for channel in ("raw", "reflected"):
                path = target / (channel + ".txt")
                if not path.is_file() or path.stat().st_size > config.max_output_bytes:
                    raise BackendOutputError("missing or oversized ordered coefficient artifact")
                payload[channel] = parse_hyperint_coefficient(path.read_text(encoding="utf-8"),
                                                               max_words=config.max_words)
            coefficients["0"] = payload
            status = "computed"
    except CASBudgetExhausted as exc:
        status, diagnostic = "resource_exhausted", str(exc)
        record["reason_code"] = exc.code
    except ProcessLimitUnavailable as exc:
        status, diagnostic = "unsupported", str(exc)
    except BackendOutputError as exc:
        status, diagnostic = "unsupported", str(exc)
    except OSError as exc:
        status, diagnostic = "failed", str(exc)
    record["status"] = status
    record["finished_utc"] = datetime.now(timezone.utc).isoformat()
    result = {"status": status, "coefficients": coefficients,
              "complete_orders": [0] if coefficients else [], "epsilon_lower_bound": 0,
              "provenance": request, "run_directory": str(target), "records": [record],
              "formal_state": "not_run"}
    if diagnostic:
        result["diagnostic"] = diagnostic
    if "reason_code" in record:
        result["reason_code"] = record["reason_code"]
    (target / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def prepare_ordered_moments(j, *, max_order=2, loop_order=None, endpoint_edges=None,
                            topology=None, variable_bindings=None):
    """Eliminate one adjacent pair and retain the complete normalized gap weight.

    Output coefficients belong to M(e,z), with J(e,z)=prefactor(e,z)*M(e,z).
    The prefactor is returned as structured powers, never baked into M_q.
    """
    if type(max_order) is not int or not 0 <= max_order <= 2:
        raise ValueError("this ordered moment adapter supports Taylor orders 0 through 2")
    analysis = analyze_exchange(j, loop_order=loop_order, endpoint_edges=endpoint_edges,
                                topology=topology, variable_bindings=variable_bindings)
    if analysis["status"] != "reduced":
        return {"status": "unsupported", "residuals": analysis["residuals"]}
    descriptor = analysis["descriptor"]
    if descriptor["adjacent_count"] != 1 or descriptor["opposite_count"] != 2:
        return {"status": "unsupported", "residuals": [{"kind": "ordered_moment_rule_scope",
                "detail": "This rule requires one adjacent pair and two distinct opposite kernels.",
                "J": deepcopy(j)}]}
    if any(len(word) > 1 for edge in descriptor["edges"] for word in edge["gap_words"]):
        return {"status": "unsupported", "residuals": [{"kind": "ordered_gap_word_not_singleton",
                "detail": "Each eliminated gap must be empty or a singleton A/B.", "J": deepcopy(j)}]}
    names = [name for name in j["variables"] if name in descriptor["remaining_bounds"]]
    name_map = {name: f"wox{index + 1}" for index, name in enumerate(names)}
    name_map.update(z1="z", eps="woEps")
    allowed = set(names) | {"eps", "z1"}

    def emit(text):
        return _emit_scalar(arithmetic(text, allowed), name_map)

    def power(base):
        # These are the analytic continuations of the Re(eps)>0 primitives.
        # In particular, 0**eps is identically zero before continuation;
        # it is never expanded as exp(eps*log(0)).
        if base == "0":
            return "0"
        if base == "1":
            return "1"
        return f"exp(woEps*ln({emit(base)}))"

    def complement(endpoint):
        return {"0": "1", "1": "0"}.get(endpoint, f"1-({endpoint})")

    gap_records, factors = [], []
    for edge in descriptor["edges"]:
        points = ["0"] + edge["marked_endpoints"] + ["1"]
        for index, word in enumerate(edge["gap_words"]):
            if not word:
                continue
            lower, upper = points[index:index + 2]
            if word == "A":
                factor = f"({power(upper)}-{power(lower)})"
            else:
                factor = f"({power(complement(lower))}-{power(complement(upper))})"
            gap_records.append({"edge": edge["edge"], "word": word,
                                "lower": lower, "upper": upper, "maple_factor": factor})
            factors.append(factor)
    polynomials = [emit(p["original_polynomial"]) for p in descriptor["propagators"]
                   if p["kind"] == "opposite"]
    weight = "*".join(factors) or "1"
    integrand = "*".join([f"(exp(woEps*ln({p}))/({p}))" for p in polynomials] + [f"({weight})"])
    execution_order = list(reversed(descriptor["remaining_order_outer_to_inner"]))
    intervals = [f'{name_map[name]}={emit(descriptor["remaining_bounds"][name][0])}..{emit(descriptor["remaining_bounds"][name][1])}'
                 for name in execution_order]
    return {"status": "ready", "operation": "ordered-opposite-gap-moments",
            "quantity": "regular_kernel_M", "max_order": max_order, "epsilon_lower_bound": -2,
            "regular_kernel_epsilon_lower_bound": 0,
            "prefactor_powers": {"opposite_factor_power": 2, "z_epsilon_power": descriptor["adjacent_z_count"],
                                 "epsilon_power": -2},
            "input_J": deepcopy(j), "input_J_sha256": descriptor["input_J_sha256"],
            "loop_order": descriptor["loop_order"], "channel_partition": descriptor["channel_partition"],
            "original_outer_to_inner": list(j["integration_order_outer_to_inner"]),
            "remaining_outer_to_inner": descriptor["remaining_order_outer_to_inner"],
            "remaining_bounds": descriptor["remaining_bounds"],
            "execution_inner_to_outer": execution_order, "variable_map": name_map,
            "normalized_gap_factors": gap_records, "maple_weight": weight,
            "maple_integrand": integrand, "maple_intervals": intervals,
            "ordinary_endpoint_admission": {
                "kind": "one_adjacent_pair_absolute_convergence_bridge",
                "context": "z1 > 0; initially 0 < Re(eps) < 1/8, then continuation to |eps| < 1/8",
                "argument": "On the initial strip, the original adjacent endpoint densities are absolutely integrable. Their two endpoints each occur in one singleton gap; Fubini within this absolute-convergence domain gives exactly the displayed primitives and preserves the relative original order of all opposite endpoints. Each remaining opposite kernel uses a distinct endpoint pair and is comparable to u+v at its two zero corners. The finite primitive differences are bounded by sums of endpoint powers; at most two such powers supplement the kernel exponent. The product cube majorant is integrable uniformly for |eps|<1/8, including logarithmic derivatives. Thus M is holomorphic and its complete Taylor coefficient may be integrated termwise. This also specifies the meromorphic continuation J=prefactor*M from the initial strip.",
                "endpoint_interpretation": "ordinary limits of each complete Taylor coefficient; no isolated finite parts",
                "formal_state": "not_run"}, "formal_state": "not_run"}


def _moments_script(prepared, config):
    intervals = ",".join(prepared["maple_intervals"])
    return f'''restart:
_hyper_autoload_periods := ["{_maple_path(config.periods)}"]:
read "{_maple_path(config.hyperint)}":
_hyper_verbosity:=0:_hyper_check_divergences:=true:_hyper_abort_on_divergence:=true:
_hyper_ignore_nonlinear_polynomials:=false:
try
  woRegular:={prepared["maple_integrand"]}:
  woTaylor:=convert(series(woRegular,woEps=0,{prepared["max_order"] + 1}),polynom):
  for woQ from 0 to {prepared["max_order"]} do
    woInput:=expand(coeff(woTaylor,woEps,woQ)):
    if woInput=0 then woRaw:=0:woRef:=0 else
      woAnswer:=hyperInt(woInput,[{intervals}]):
      woRaw:=fibrationBasis(woAnswer,[z]):
      woRef:=subs(wozz=z,fibrationBasis(convert(subs(z=1/wozz,woRaw),HlogRegInf),[wozz])):
    fi:
    for woChannel in ["raw","reflected"] do
      if woChannel="raw" then woValue:=woRaw else woValue:=woRef fi:
      if type(woValue,`+`) then woCount:=nops(woValue) else woCount:=1 fi:
      if woCount>{config.max_words} then error "HLOG word budget exceeded" fi:
      woFile:=fopen(cat("q",woQ,".",woChannel,".txt"),WRITE,TEXT):
      fprintf(woFile,"%a\\n",woValue):fclose(woFile):
    od:
    printf("WILSON_ORDERED_MOMENT_%d_COMPLETE\\n",woQ):
    gc():
  od:
  printf("WILSON_ORDERED_MOMENTS_COMPLETE\\n"):
catch:
  printf("WILSON_ORDERED_MOMENTS_FAILED %a\\n",[lastexception]):
end try:
quit:
'''


def compute_ordered_moments(j, *, run_dir, max_order=2, loop_order=None,
                            endpoint_edges=None, topology=None, variable_bindings=None,
                            config: MapleConfig | None = None, cas_budget=None):
    """Compute complete regular-kernel M_0..M_q in a fresh bounded process.

    Multiply raw by prefactor(e,z) and reflected by prefactor(e,1/z)
    before applying any physical I normalization. The output is not raw J_q.
    """
    prepared = prepare_ordered_moments(j, max_order=max_order, loop_order=loop_order,
        endpoint_edges=endpoint_edges, topology=topology, variable_bindings=variable_bindings)
    if prepared["status"] != "ready":
        return {**prepared, "coefficients": {}, "complete_orders": [], "formal_state": "not_run"}
    config = config or find_maple_config()
    if config is None:
        return {"status": "unsupported", "coefficients": {}, "complete_orders": [],
                "diagnostic": "Maple/HyperInt unavailable", "provenance": prepared, "formal_state": "not_run"}
    target = Path(run_dir).resolve() / ("ordered-moments-" + uuid4().hex)
    target.mkdir(parents=True, exist_ok=False)
    request = {**prepared, "backend": config.fingerprint(), "timeout_seconds": config.timeout_seconds,
               "memory_mb": config.memory_mb, "max_words": config.max_words}
    (target / "request.json").write_text(json.dumps(request, indent=2) + "\n", encoding="utf-8")
    script_path = target / "compute.mpl"
    script_path.write_text(_moments_script(prepared, config), encoding="utf-8", newline="\n")
    command = [str(Path(config.executable).resolve()), "-q", str(script_path)]
    record = {"directory": str(target), "command": command, "script_sha256": _hash(script_path),
              "started_utc": datetime.now(timezone.utc).isoformat()}
    status, diagnostic, coefficients = "failed", None, {}
    try:
        runner = run_bounded_process if cas_budget is None else (
            lambda *a, **kw: cas_budget.run(run_bounded_process, *a, **kw))
        execution = runner(command, cwd=target, stdout_path=target / "maple.log",
                           timeout_seconds=config.timeout_seconds, memory_mb=config.memory_mb)
        record.update(execution)
        log_text = (target / "maple.log").read_text(encoding="utf-8", errors="replace")
        if execution["stop_reason"] in ("time", "memory"):
            status, diagnostic = "resource_exhausted", "ordered moments exceeded process " + execution["stop_reason"] + " budget"
        elif (execution["returncode"] != 0 or "WILSON_ORDERED_MOMENTS_COMPLETE" not in log_text
              or "WILSON_ORDERED_MOMENTS_FAILED" in log_text or "Error," in log_text):
            if "HLOG word budget exceeded" in log_text or "alloc" in log_text.lower() or "MemoryError" in log_text:
                status, diagnostic = "resource_exhausted", "ordered moments exceeded word or memory budget"
            elif "divergen" in log_text.lower() or "FAILED: factors" in log_text:
                status, diagnostic = "unsupported", "the current finite-interval route could not establish ordinary endpoints or linear factorization"
            else:
                status, diagnostic = "failed", "CAS did not complete the ordered moment request"
            record["log_tail"] = log_text[-4000:]
        else:
            candidate = {}
            for order in range(max_order + 1):
                payload = {}
                for channel in ("raw", "reflected"):
                    path = target / f"q{order}.{channel}.txt"
                    if not path.is_file() or path.stat().st_size > config.max_output_bytes:
                        raise BackendOutputError("missing or oversized ordered moment artifact")
                    payload[channel] = parse_hyperint_coefficient(path.read_text(encoding="utf-8"),
                                                                   max_words=config.max_words)
                candidate[str(order)] = payload
            coefficients, status = candidate, "computed"
    except CASBudgetExhausted as exc:
        status, diagnostic = "resource_exhausted", str(exc)
        record["reason_code"] = exc.code
    except ProcessLimitUnavailable as exc:
        status, diagnostic = "unsupported", str(exc)
    except BackendOutputError as exc:
        status, diagnostic = "unsupported", str(exc)
    except OSError as exc:
        status, diagnostic = "failed", str(exc)
    record.update(status=status, finished_utc=datetime.now(timezone.utc).isoformat())
    result = {"status": status, "coefficients": coefficients,
              "complete_orders": list(range(max_order + 1)) if coefficients else [],
              "epsilon_lower_bound": -2, "regular_kernel_epsilon_lower_bound": 0,
              "quantity": "regular_kernel_M",
              "prefactor_powers": prepared["prefactor_powers"],
              "provenance": request, "run_directory": str(target), "records": [record],
              "formal_state": "not_run"}
    if diagnostic:
        result["diagnostic"] = diagnostic
    if "reason_code" in record:
        result["reason_code"] = record["reason_code"]
    (target / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result
