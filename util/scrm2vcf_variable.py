#!/usr/bin/env python
from __future__ import print_function, division
import shutil
import sh
import os
import sys
import argparse
# try:
#     from shutil import which
# except ImportError:
#     from backports.shutil_which import which

#import smcpp.util
#from scrm import demography_from_params

SCRM = os.environ.get('SCRM_PATH', False) #or which('scrm')

if __name__ == "__main__":
    if not SCRM:
        sys.exit("Can't find scrm. Please set SCRM_PATH.")
    scrm = sh.Command(SCRM)
    parser = argparse.ArgumentParser()
    parser.add_argument("--contig", default="contig1", help="name of contig in VCF")
    #parser.add_argument("--demography", choices=["human", "sawtooth"])
    parser.add_argument("-o", help="output location (default: stdout)")
    parser.add_argument("n", type=int, help="diploid sample size")
    parser.add_argument("-rho", type=float, help="initial recombination rate.")
    #parser.add_argument("--var_rho", type=str, help = "Space separated variable rho. Using scrm format.")
    parser.add_argument("length", type=int, help="length of chromosome to simulate")
    parser.add_argument("-theta", type=float, help = "mutation rate")
    parser.add_argument("replicates", type = int, help = "number of independent replicates")

    # --contig tenx --demography 1 -o "./data" -n 10 --length 50000 --theta 5e-8 --replicates 5 --rho 4.4e-8 --var_rho 24000 4.4e-7 26000 4.4e-8

    args, scrm_extra_args = parser.parse_known_args()

    # if args.demography is not None:
    #     demo = getattr(smcpp.util, args.demography)
    #     a = demo['a']
    #     b = demo['b']
    #     s = demo['s'] * 0.5
    #     scrm_extra_args += demography_from_params((a, b, s))

    #rho_changes = args.var_rho.split()
    #if len(rho_changes) % 2 == 1:
    #    sys.exit("Incorrect variable rho params.")

    scrm_args = [2 * args.n, args.replicates]
    scrm_args.append("--transpose-segsites")
    scrm_args += ["-SC", "abs", "-p", 14]
    scrm_args += ["-r", args.rho, args.length]
    scrm_args += ["-t", args.theta]
    #for i in range(int(len(rho_changes) / 2)):
    #    scrm_args += ["-sr", float(rho_changes[i]), float(rho_changes[i+1])] #TOFIX
    scrm_args += scrm_extra_args

    # Iterate over scrm output
    print("Calling scrm with args: %s" % str(scrm_args), file=sys.stderr)
    it = scrm(*scrm_args, _iter=True)
    line = next(it)
    print(line, file=sys.stdout)
    while not line.startswith("//"): # leave iterator at //
        line = next(it)
        print(line, file=sys.stdout)

    # Create new files for each replicate
    for rep in range(args.replicates):
        print("Processing rep %d" % rep, file=sys.stderr)
        # define files
        out_name = args.o + "_rep" + str(rep)
        out = open(out_name, "wt")
        sites = out_name + "_sites"
        loc = out_name + "_loc"
        out_sites = open(sites, "wt")
        out_loc = open(loc, "wt")

        # Preprocess whitespace
        while not line.startswith("transposed segsites"):
            line = next(it)
            print(line, file=sys.stdout)
        segsites = float(line.strip().split()[-1])
        line = next(it) # skip the position line
        print(line, file=sys.stdout)

        # Create sites and locs header
        sites_header = [str(args.n), str(segsites), "2"]
        loc_header = [str(segsites), str(args.length), "L"] #crossing-over model
        print(" ".join(sites_header), file=out_sites)
        print(" ".join(loc_header), file=out_loc)

        # Create a minimal VCF header
        header = ["##fileformat=VCFv4.0", """##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""]
        header.append("##contig=<ID={},length={}>".format(args.contig, args.length))
        h = "#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split()
        h += ["sample%d" % i for i in range(1, args.n + 1)]
        header.append("\t".join(h))
        print("\n".join(header), file=out)

        print("Starting vcf and loc...", file=sys.stderr)
        seqs = [""]*args.n
        for line in it:
            print(line, file=sys.stdout)
            ary = line.strip().split()
            if len(ary) < 2:
                break
            pos, time = ary[:2]

            ## Write to VCF file
            gts = ary[2:]
            cols = [args.contig, str(int(float(pos))), ".", "A", "C", ".", "PASS", ".", "GT"]
            cols += ["|".join(gt) for gt in zip(gts[::2], gts[1::2])]
            print("\t".join(cols), file=out)

            ## Write to locs file
            print(pos + "\t", file=out_loc)

            ## Collect sites information
            gts = [int(i) for i in gts]
            dip = [sum(gt) for gt in zip(gts[::2], gts[1::2])]
            dip = [1 if i == 2 else 2 if i == 1 else 0 for i in dip] # ldhat represents heterozygotes as 2
            seqs = [i[0] + str(i[1]) for i in zip(seqs, dip)]

        print("Starting sites...", file=sys.stderr)
        ## Write to sites file
        for i in range(1, args.n + 1):
            print(">sample%d" % i, file=out_sites)
            print(seqs[i-1], file=out_sites)
        print("Finished with rep %d" % rep, file=sys.stderr)
        out.close()
        out_sites.close()
        out_loc.close()

