import argparse,json,time
from .prepare import prepare

def main():
    parser=argparse.ArgumentParser(description='Prepare color- and Euler-normalized C/J inputs without Lean')
    parser.add_argument('--output')
    args=parser.parse_args();start=time.perf_counter()
    result=prepare(output=args.output)
    print(json.dumps(dict(status='prepared',recipe=result['recipe'],euler_compensation=result['euler_compensation'],total_orbits=result['total_orbits'],unique_colors=result['unique_color_expressions'],processing_seconds=time.perf_counter()-start,lean_invoked=False,output=args.output or 'benchmarks/2ne/prepared'),indent=2))
if __name__=='__main__':main()
