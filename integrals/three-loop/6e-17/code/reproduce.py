"""Run the whole39-entry same-batch generation, then select this member."""
from pathlib import Path
import runpy,sys
runtime=Path(__file__).resolve().parents[3]/'families/three-loop-6e/runtime'
sys.path.insert(0,str(runtime/'scripts'))
runpy.run_path(str(runtime/'scripts/evaluate_batch.py'),run_name='__main__')
