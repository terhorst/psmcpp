#!/bin/bash

source /export/home/users/jchan/opt/terhorst/bin/activate
export SCRM_PATH="/export/home/users/jchan/scrm"
python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf_variable.py 10 50000 5 -theta 68.75 -rho 60.5 -o "tenx" -sr 24000 605 -sr 26000 60.5

