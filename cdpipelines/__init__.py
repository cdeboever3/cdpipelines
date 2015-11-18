import os

import atacseq, general, prepare, rnaseq

_root = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]

_scripts = os.path.join(
    os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-2]),
    'scripts')
