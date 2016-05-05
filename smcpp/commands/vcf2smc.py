import argparse
import warnings
import itertools as it
from logging import getLogger
import numpy as np
import sys
from pysam import VariantFile, TabixFile
logger = getLogger(__name__)

from ..logging import init_logging
from ..util import optional_gzip, RepeatingWriter


def init_parser(parser):
    parser.add_argument("--ignore-missing", default=False, action="store_true",
            help="ignore samples which are missing in the data")
    parser.add_argument("-i", "--distinguished_index", type=int, default=0, 
            help="index of distinguished lineage in sample ids (default: 0)")
    parser.add_argument("--missing-cutoff", "-c", metavar="c", type=int, default=None,
            help="treat runs of homozygosity longer than <c> base pairs as missing")
    parser.add_argument("--mask", "-m", help="BED-formatted mask of missing regions", widget="FileChooser")
    parser.add_argument("vcf", metavar="vcf[.gz]", help="input VCF file", widget="FileChooser")
    parser.add_argument("out", metavar="out[.gz]", help="output SMC++ file", widget="FileChooser")
    parser.add_argument("region", help="SAM-style region to parse")
    parser.add_argument("sample_ids", nargs="+", help="Column(s) to pull from the VCF, or file containing the same.")

def validate(args):
    if args.missing_cutoff and args.mask:
        raise RuntimeError("--missing-cutoff and --mask are mutually exclusive")

def main(args):
    init_logging(".", False)
    validate(args)
    if len(args.sample_ids) == 1:
        try:
            with open(args.sample_ids[0], "rt") as f:
                args.sample_ids = [line.strip() for line in f]
        except IOError:
            pass
    dist = args.sample_ids[args.distinguished_index]
    undist = [sid for j, sid in enumerate(args.sample_ids) if j != args.distinguished_index]
    logger.info("Distinguished sample: " + dist)
    logger.info("Undistinguished samples: " + ",".join(undist))
    vcf = VariantFile(args.vcf)
    with optional_gzip(args.out, "wt") as out:
        samples = list(vcf.header.samples)
        if dist not in samples:
            raise RuntimeError("Distinguished lineage not found in data?")
        missing = [u for u in undist if u not in samples]
        if missing:
            msg = "The following samples were not found in the data: %s. " % ", ".join(missing)
            if args.ignore_missing:
                logger.warn(msg)
            else:
                msg += "If you want to continue without these samples, use --ignore-missing."
                raise RuntimeError(msg)
        undist = [u for u in undist if u not in missing]
        nb = 2 * len(undist)

        # function to convert a VCF record to our format <span, dist gt, undist gt, # undist>
        def rec2gt(rec):
            ref = rec.alleles[0]
            if None in rec.samples[dist].alleles:
                a = -1
            else:
                a = sum(allele != ref for allele in rec.samples[dist].alleles)
            bs = [allele != ref for u in undist for allele in rec.samples[u].alleles if allele is not None]
            b = sum(bs)
            nb = len(bs)
            return [a, b, nb]

        region_iterator = vcf.fetch(region=args.region)
        if args.mask:
            mask_iterator = TabixFile(args.mask).fetch(region=args.region)
            args.missing_cutoff = np.inf
        else:
            mask_iterator = iter([])
        mask_iterator = (x.split("\t") for x in mask_iterator)
        mask_iterator = ((x[0], int(x[1]), int(x[2])) for x in mask_iterator)
        snps_only = (rec for rec in region_iterator if len(rec.alleles) == 2 and set(rec.alleles) <= set("ACTG"))

        def interleaved():
            cmask = next(mask_iterator, None)
            csnp = next(snps_only, None)
            while cmask or csnp:
                if cmask is None:
                    yield "snp", csnp
                    csnp = next(snps_only, None)
                elif csnp is None:
                    yield "mask", cmask
                    cmask = next(mask_iterator, None)
                else:
                    if csnp.pos < cmask[1]:
                        yield "snp", csnp
                        csnp = next(snps_only, None)
                    elif csnp.pos <= cmask[2]:
                        while csnp is not None and csnp.pos <= cmask[2]:
                            csnp = next(snps_only, None)
                        yield "mask", cmask
                        cmask = next(mask_iterator, None)
                    else:
                        yield "mask", cmask
                        cmask = next(mask_iterator, None)

        with RepeatingWriter(out) as rw:
            records = interleaved()
            ty, rec = next(records)
            if ty == "snp":
                last_pos = rec.pos
                rw.write([1] + rec2gt(rec))
            else:
                last_pos = rec[2]
                rw.write([rec[2] - rec[1] + 1, -1, 0, 0])
            for ty, rec in records:
                if ty == "mask":
                    span = rec[1] - last_pos
                    rw.write([span, 0, 0, nb])
                    rw.write([rec[2] - rec[1] + 1, -1, 0, 0])
                    last_pos = rec[2]
                    continue
                abnb = rec2gt(rec)
                span = rec.pos - last_pos - 1
                if 1 <= span <= args.missing_cutoff:
                    rw.write([span, 0, 0, nb])
                elif span > args.missing_cutoff:
                    rw.write([span, -1, 0, 0])
                rw.write([1] + abnb)
                last_pos = rec.pos
