"""Run the shared same-loop calculation; read this integral's member."""
from pathlib import Path
import runpy,sys
root=Path(__file__).resolve().parents[3]
runtime=root/'families/low-loop-runtime/integrals/two-loop'
sys.path.insert(0,str(runtime))
runpy.run_path(str(runtime/'reproduce.py'),run_name='__main__')
