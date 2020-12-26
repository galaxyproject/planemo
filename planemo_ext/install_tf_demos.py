import argparse
import os
import subprocess
import sys
import urllib.request

from bioblend import galaxy

WF = 'https://drive.google.com/uc?export=download&id=13xE8o7tucHGNA0qYkEP98FfUGl2wdOU5'
HIST = 'https://drive.google.com/uc?export=download&id=1V0ZN9ZBuqcGJvt2AP7s3g0q11uYEhdDB'
WF_FILE = 'tf_workflow.ga'
HIST_FILE = 'tf_history.tgz'

def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--galaxy", help='URL of target galaxy', default="http://localhost:9090")
    parser.add_argument("-a", "--key", help='Galaxy admin key', default=None)
    return parser

def main():
    """
    load the planemo tool_factory demonstration history and tool generating workflow
    """
    args = _parser().parse_args()
    urllib.request.urlretrieve(WF, WF_FILE)
    urllib.request.urlretrieve(HIST,HIST_FILE)
    assert args.key, 'Need an administrative key for the target Galaxy supplied please'
    wfp = os.path.abspath(WF_FILE)
    hp = os.path.abspath(HIST_FILE)
    gi = galaxy.GalaxyInstance(url=args.galaxy, key=args.key, email="planemo@galaxyproject.org")
    x = gi.workflows.import_workflow_from_local_path(WF_FILE, publish=True)
    print(f'installed {WF_FILE} Returned = {x}\n')
    x = gi.histories.import_history(file_path=HIST_FILE)
    print(f'installed {HIST_FILE} Returned = {x}\n')


if __name__ == "__main__":
    main()
