"""Run QGRAF binary and return the output file path."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def run_qgraf(control_path: Path, qgraf_bin: str = "qgraf") -> Path:
    control_path = Path(control_path).resolve()
    work_dir = control_path.parent

    result = subprocess.run(
        [qgraf_bin],
        capture_output=True, text=True,
        cwd=str(work_dir),
        input=control_path.read_text(),
        timeout=300,
    )

    if result.returncode != 0:
        print(f"QGRAF stderr:\n{result.stderr}", file=sys.stderr)
        raise RuntimeError(f"QGRAF failed with code {result.returncode}")

    control_text = control_path.read_text()
    out_name = "qgraf.out"
    for line in control_text.splitlines():
        line = line.strip()
        if line.startswith("output"):
            out_name = line.split("=")[1].strip().strip("'").strip('"')
            break

    out_path = work_dir / out_name
    if not out_path.exists():
        raise FileNotFoundError(f"QGRAF output not found: {out_path}")
    return out_path
