"""Parse QGRAF mathematica.sty .out file into Diagram dataclasses.

Supports arbitrary n-loop by auto-detecting loop momentum symbols
and using Dict[str, int] for momentum arithmetic.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List


@dataclass
class Prop:
    field1: str
    field2: str
    mom: str


@dataclass
class Leg:
    field: str
    mom: str


@dataclass
class Vert:
    legs: List[Leg]


@dataclass
class Diagram:
    idx: int
    sign: int
    sym_factor: str
    props: List[Prop] = field(default_factory=list)
    verts: List[Vert] = field(default_factory=list)
    # Optional tag set by the combinatorial generator recording the
    # emission-type name this diagram was built for.  When present,
    # classify_orbits() uses it directly instead of re-deriving the
    # emission type from the (synthetic) momentum structure, which is
    # essential for loop-type vertices (Pi/Gamma) that carry no loop
    # momentum symbols in the synthetic representation.
    emission_type: str = ""


def mom_dict(mom: str) -> Dict[str, int]:
    """Parse signed momentum string into {symbol: coefficient}.

    "k1+k2-p1" -> {"k1": 1, "k2": 1, "p1": -1}
    """
    s = mom.replace(" ", "")
    if not s or s == "0":
        return {}
    if not (s.startswith("+") or s.startswith("-")):
        s = "+" + s
    coeffs: Dict[str, int] = {}
    i = 0
    while i < len(s):
        sign = 1
        if s[i] == "+":
            sign = 1; i += 1
        elif s[i] == "-":
            sign = -1; i += 1
        j = i
        while j < len(s) and (s[j].isalnum() or s[j] == "_"):
            j += 1
        if j > i:
            sym = s[i:j]
            coeffs[sym] = coeffs.get(sym, 0) + sign
            i = j
        else:
            i += 1
    return {k: v for k, v in coeffs.items() if v != 0}


def mom_eq(a: str, b: str) -> bool:
    return mom_dict(a) == mom_dict(b)


def mom_neg(a: str, b: str) -> bool:
    ad = mom_dict(a)
    bd = mom_dict(b)
    if set(ad.keys()) != set(bd.keys()):
        return False
    return all(ad[k] == -bd[k] for k in ad)


_RE_PROP = re.compile(
    r"prop\(\s*([A-Za-z0-9]+)\s*,\s*([A-Za-z0-9]+)\s*,\s*([^)]+)\)")
_RE_AMP_HEAD = re.compile(
    r"amp\((\d+)\)\s*=\s*\((-?\+?\-?1)\)\s*\*\s*\(\s*([^)]+)\s*\)\s*\*")
_RE_DIAG_START = re.compile(r"^%\s*diagram\s+(\d+)\s*$", re.MULTILINE)


def _split_top_commas(text: str) -> List[str]:
    out: List[str] = []
    depth = 0
    cur = []
    for ch in text:
        if ch == "(":
            depth += 1; cur.append(ch)
        elif ch == ")":
            depth -= 1; cur.append(ch)
        elif ch == "," and depth == 0:
            out.append("".join(cur).strip()); cur = []
        else:
            cur.append(ch)
    tail = "".join(cur).strip()
    if tail:
        out.append(tail)
    return out


def _extract_vert_bodies(text: str) -> List[str]:
    bodies: List[str] = []
    i = 0; n = len(text)
    while i < n:
        j = text.find("vert(", i)
        if j == -1:
            break
        start = j + 5; depth = 1; k = start
        while k < n and depth > 0:
            if text[k] == "(":
                depth += 1
            elif text[k] == ")":
                depth -= 1
            k += 1
        if depth == 0:
            bodies.append(text[start:k - 1]); i = k
        else:
            break
    return bodies


def _parse_leg(text: str) -> Leg:
    m = re.match(r"\s*([A-Za-z0-9]+)\s*\(\s*(.+?)\s*\)\s*$", text)
    if not m:
        return Leg(field=text.strip(), mom="0")
    return Leg(field=m.group(1), mom=m.group(2).strip())


def parse_file(path: str | Path) -> List[Diagram]:
    text = Path(path).read_text(encoding="utf-8")
    headers = [(m.start(), int(m.group(1)))
               for m in _RE_DIAG_START.finditer(text)]
    headers.append((len(text), -1))
    diagrams: List[Diagram] = []
    for (start, idx), (stop, _) in zip(headers[:-1], headers[1:]):
        block = text[start:stop]
        m_amp = re.search(
            r"amp\(\d+\)\s*=\s*\(\s*(-?\+?\-?1)\s*\)\s*\*\s*\(\s*([^)]+)\s*\)\s*\*",
            block)
        if not m_amp:
            continue
        sign = -1 if m_amp.group(1).strip().startswith("-") else 1
        sym = m_amp.group(2).strip()
        body = block[m_amp.end():].split("1;", 1)[0]
        flat = re.sub(r"\s+", "", body)

        props = [Prop(field1=m.group(1), field2=m.group(2), mom=m.group(3))
                 for m in _RE_PROP.finditer(flat)]

        verts = []
        for vbody in _extract_vert_bodies(flat):
            legs = [_parse_leg(t) for t in _split_top_commas(vbody)]
            verts.append(Vert(legs=legs))

        diagrams.append(Diagram(
            idx=idx, sign=sign, sym_factor=sym, props=props, verts=verts))
    return diagrams
