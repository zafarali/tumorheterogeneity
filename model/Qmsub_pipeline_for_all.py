from __future__ import print_function
import sys
from glob import glob
import os 
## in super run:
# print(sys.argv)
valid_file_paths = map(os.path.dirname, glob(sys.argv[1]+'*/all_*'))
# print(valid_file_paths)
TEMPLATE = "Qmsub -h 16 -n 4 -q sw pipeline_for_one.sh {filepath}/"
for filepath in valid_file_paths:
    print(TEMPLATE.format(filepath=filepath))
