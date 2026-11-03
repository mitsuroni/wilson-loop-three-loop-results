"""Thin JSON CLI. Logs never share stdout with the structured result."""
import argparse
import json
import sys
from interfaces.validate import load
from . import solve
from .contracts import result_shell


def main(argv=None):
    parser=argparse.ArgumentParser(description='Solve an ordered integral to a declared epsilon order')
    parser.add_argument('request',help='versioned JSON request file')
    args=parser.parse_args(argv)
    try:result=solve(load(args.request))
    except (OSError,ValueError,TypeError,KeyError) as exc:
        result=result_shell({});result.update(status='failed',reason={'code':'invalid_input','message':str(exc)})
    print(json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False))
    return {'computed':0,'partial':2,'unsupported':2,'resource_exhausted':3,'failed':1}[result['status']]


if __name__=='__main__':sys.exit(main())
