"""Run the shared same-loop calculation; read this integral's member."""
from pathlib import Path
import runpy
root=Path(__file__).resolve().parents[3]
runpy.run_path(str(root/'families/low-loop-runtime/integrals/one-loop/reproduce.py'),run_name='__main__')
