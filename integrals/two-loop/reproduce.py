from pathlib import Path
import sys,runpy
runtime=Path(__file__).resolve().parents[1]/"families/low-loop-runtime/integrals/two-loop"
sys.path.insert(0,str(runtime))
runpy.run_path(str(runtime/"reproduce.py"),run_name="__main__")
