#!/bin/bash

source /export/home/users/jchan/opt/terhorst_new/bin/activate
export SCRM_PATH="/export/home/users/jchan/scrm"

#python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 2200.0 100000 -o "flat_5x" --demography "human" -t 500.0
#python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 8800.0 100000 -o "flat_20x" --demography "human" -t 500.0
#python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 22000.0 100000 -o "flat_50x" --demography "human" -t 500.0
#python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 33000.0 100000 -o "flat_75x" --demography "human" -t 500.0
#python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 44000.0 100000 -o "flat_100x" --demography "human" -t 500.0
#python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 88000.0 100000 -o "flat_200x" --demography "human" -t 500.0
#python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 220000.0 100000 -o "flat_500x" --demography "human" -t 500.0
for i in {0..19}
do
   python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 440.0 100000 -o "tenx_5k_$i" --demography "human" -t 500.0 -sr 47500 4400.0 -sr 52500 440.0
done
#python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 440.0 100000 -o "fiftyx_15k" --demography "human" -t 500.0 -sr 42500 22000.0 -sr 57500 440.0
#python /export/home/users/jchan/terhorst_psmcpp/psmcpp/util/scrm2vcf.py 25 440.0 100000 -o "hundredx_15k" --demography "human" -t 500.0 -sr 42500 44000.0 -sr 57500 440.0

