import argparse
import warnings
import itertools as it
from logging import getLogger
import numpy as np
import sys
import pysam
from pysam import VariantFile, TabixFile
import json
from collections import Counter, namedtuple
import progressbar
logger = getLogger(__name__)

from ..logging import setup_logging
from ..util import optional_gzip, RepeatingWriter
from ..version import __version__

print 'HELLO!!!'
print pysam.__file__
#print VariantFile.__file__
#print TabixFile.__file__
SampleList = namedtuple("SampleList", "pid samples")


def sample_list(x):
    try:
        x1, x2 = x.split(":")
        return SampleList(x1, x2.split(","))
    except:
        raise argparse.ArgumentTypeError(
                "%r should be a comma-separated list of sample ids preceded by a "
                "population identifier. See 'smc++ vcf2smc -h'." % x)


def init_parser(parser):
    parser.add_argument("-d", nargs=2, help=argparse.SUPPRESS)
    parser.add_argument("--ignore-missing", default=False, action="store_true",
            help="ignore samples which are missing in the data")
    parser.add_argument("--missing-cutoff", "-c", metavar="c", type=int, default=None,
            help="treat runs of homozygosity longer than <c> base pairs as missing")
    parser.add_argument("--mask", "-m",
            help="BED-formatted mask of missing regions",
            widget="FileChooser")
    parser.add_argument("vcf", metavar="vcf[.gz]",
            help="input VCF file", widget="FileChooser")
    parser.add_argument("out", metavar="out[.gz]",
            help="output SMC++ file", widget="FileChooser")
    parser.add_argument("contig", help="contig to parse")
    parser.add_argument("pop1", type=sample_list,
            help="List of sample ids from population 1. "
                 "Format: <pop_id>:<sample_id_1>,<sample_id_2>,...,<sample_id_N>") 
    parser.add_argument("pop2", type=sample_list, nargs="?", default=SampleList(None, []),
            help="List of sample ids from population 2, in same format.")


def validate(args):
    if args.missing_cutoff and args.mask:
        raise RuntimeError("--missing-cutoff and --mask are mutually exclusive")


def main(args):
    setup_logging(0)
    validate(args)
    for i in [1, 2]:
        attr = "pop%d" % i
        pid, ary = getattr(args, attr)
        if len(ary) == 1 and ary[0][0] == "@":
            setattr(args, attr, SampleList(pid, open(ary[0][1:], "rt").read().strip().split("\n")))
    pop_d = dict([args.pop1, args.pop2])
    for pid in pop_d:
        if pop_d[pid]:
            c = Counter(pop_d[pid])
            if max(c.values()) > 1:
                raise RuntimeError(
                        "Population %s has duplicated samples: %s" %
                        (pid, [item for item in c.items() if item[1] > 1]))
    dist = [[], []]
    if not args.d:
        first_sid = args.pop1.samples[0]
        args.d = [first_sid + ":0", first_sid + ":1"]
    all_samples = set(args.pop1.samples) | set(args.pop2.samples)
    for sid_i in args.d:
        sid, i = sid_i.split(":")
        i = int(i)
        if sid not in all_samples:
            raise RuntimeError("%s is not in the sample list" % sid)
        if sid in args.pop1.samples:
            d = dist[0]
        else:
            assert sid in args.pop2.samples
            d = dist[1]
        d.append((sid, i))
    undist = [[(k, i) for k in p.samples for i in (0, 1) if (k, i) not in d]
              for p, d in zip((args.pop1, args.pop2), dist)]
    npop = 1
    def print_pop(i):
        logger.info("Population %d:" % i)
        logger.info("Distinguished lineages: " + ", ".join("%s:%d" % t for t in dist[i - 1]))
        logger.info("Undistinguished lineages: " + ", ".join("%s:%d" % t for t in undist[i - 1]))
    print_pop(1)
    if args.pop2.pid is not None:
        npop = 2
        common = set(args.pop1.samples) & set(args.pop2.samples)
        if common:
            logger.error("Populations 1 and 2 should be disjoint, "
                         "but both contain " + ", ".join(common))
            sys.exit(1)
        print_pop(2)

    ## Start parsing
    vcf = VariantFile(args.vcf)
    with optional_gzip(args.out, "wt") as out:
        samples = list(vcf.header.samples)
        dist = dist[:npop]
        undist = undist[:npop]
        if not set([dd[0] for d in dist for dd in d]) <= set(samples):
            raise RuntimeError("Distinguished lineages not found in data?")
        missing = [s for u in undist for s, _ in u if s not in samples]
        if missing:
            msg = "The following samples were not found in the data: %s. " % ", ".join(missing)
            if args.ignore_missing:
                logger.warn(msg)
            else:
                msg += "If you want to continue without these samples, use --ignore-missing."
                raise RuntimeError(msg)
        undist = [[t for t in u if t[0] not in missing] for u in undist]

        # Write header
        pids = [a.pid for a in (args.pop1, args.pop2)[:npop]]
        out.write("# SMC++ ")
        json.dump({"__version__": __version__, "pids": pids, "undist": undist, "dist": dist}, out)
        out.write("\n")
        na = list(map(len, dist))
        nb = list(map(len, undist))

        # function to convert a VCF record to our format <span, dist gt, undist gt, # undist>
        def rec2gt(rec):
            ref = rec.alleles[0]
            da = [[rec.samples[d].alleles[i] for d, i in di] for di in dist]
            a = [sum([x != ref for x in d]) if None not in d else -1 for d in da]
            bs = [[rec.samples[d].alleles[i] != ref for d, i in un if rec.samples[d].alleles[i] is not None] 
                  for un in undist]
            b = [sum(_) for _ in bs]
            nb = [len(_) for _ in bs]
            # Fold non-polymorphic (in subsample) sites
            if np.array_equal(b, nb) and np.array_equal(a, na):
                a = [0] * len(a)
                b = [0] * len(b)
            return list(sum(zip(a, b, nb), tuple()))
        region_iterator = vcf.fetch()
        #region_iterator = vcf.fetch(contig=args.contig, until_eof=True)
        contig_length = vcf.header.contigs[args.contig].length
        if args.mask:
            mask_iterator = TabixFile(args.mask).fetch(reference=args.contig)
            args.missing_cutoff = np.inf
        else:
            mask_iterator = iter([])
            if args.missing_cutoff is None:
                args.missing_cutoff = np.inf
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

        abnb_miss = [-1, 0, 0] * len(nb)
        abnb_nonseg = sum([[0, 0, x] for x in nb], [])
        multiples = set()
        with RepeatingWriter(out) as rw, progressbar.ProgressBar(max_value=contig_length) as bar:
            last_pos = 0
            for ty, rec in interleaved():
                if ty == "mask":
                    span = rec[1] - last_pos
                    rw.write([span] + abnb_nonseg)
                    rw.write([rec[2] - rec[1] + 1] + abnb_miss)
                    last_pos = rec[2]
                    continue
                bar.update(rec.pos)
                abnb = rec2gt(rec)
                if rec.pos == last_pos:
                    multiples.add(rec.pos)
                    continue
                span = rec.pos - last_pos - 1
                if 1 <= span <= args.missing_cutoff:
                    rw.write([span] + abnb_nonseg)
                elif span > args.missing_cutoff:
                    rw.write([span] + abnb_miss)
                rw.write([1] + abnb)
                last_pos = rec.pos
            rw.write([contig_length - last_pos] + abnb_nonseg)
        if multiples:
            # FIXME: what to do with multiple records at same site
            logger.warn("Multiple entries found at %d positions; skipped all but the first", len(multiples))
