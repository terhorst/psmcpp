#!/bin/bash
# qsub -V -cwd -pe openmp (number of cores)
python ./test/unit/test_stitch.py