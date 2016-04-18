import os

import atacseq, general, prepare, rnaseq
import convert_bed_to_saf

_root = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]

_scripts = os.path.join(os.path.split(os.path.abspath(__file__))[0], 'scripts')
# _scripts = os.path.join(os.path.abspath(__file__), 'scripts')
# _scripts = os.path.join(
#     os.path.sep.join(os.path.abspath(__file__).split(os.path.sep)[0:-2]),
#     'scripts')
