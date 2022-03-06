import shlex
import os
import sys
from argparse import ArgumentParser, SUPPRESS
from collections import defaultdict
from shared.intervaltree.intervaltree import IntervalTree

import random
import trio.param_t as param
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import subprocess_popen, IUPAC_base_to_ACGT_base_dict as BASE2BASE, IUPAC_base_to_num_dict as BASE2NUM
from shared.utils import str2bool





def Check(args):
    call_vcf = args.call_vcf
    true_vcf = args.true_vcf
    contig_name = args.ctgName
    bed_fn = args.bed_fn
    tree = bed_tree_from(bed_fn, contig_name=contig_name)
    is_tree_empty = len(tree.keys()) == 0

    print(' *** checking the de novo variants sites, required merged VCF, and trio in order of child, parent1, parent2 at last 3 columns')
    c_col, p1_col, p2_col = 9, 10, 11

    call_de_novo_s = set()
    if call_vcf and os.path.exists(call_vcf):
        f = subprocess_popen(shlex.split("gzip -fdc %s" % (call_vcf)))
        for row in f.stdout:
            if row[0] == '#':
                continue
            columns =row.rstrip().split('\t')
            ctg_name = columns[0]
            if contig_name and contig_name != ctg_name:
                continue
            pos = int(columns[1])
            ref_base = columns[3]
            alt_base = columns[4]
            qual = float(columns[5])

            if not (is_tree_empty or is_region_in(tree=tree, contig_name=ctg_name, \
                        region_start=pos - 1, region_end=pos)):
                continue

            def get_gt(x):
                tmp_s = x.split(':')[0]
                gt_1, gt_2 = int(tmp_s[0]), int(tmp_s[2])
                return gt_1 + gt_2
            try:
                gt_c, gt_p1, gt_p2 = get_gt(columns[c_col]), get_gt(columns[p1_col]), get_gt(columns[p2_col])
            except:
                print('error, please input merged trio VCF')
                return 0
            if gt_c == 1 and gt_p1 == 0 and gt_p2 == 0:
                call_de_novo_s.add(pos)
            # else:
            #     continue
            # print(gt_c, gt_p1, gt_p2)
            
        f.stdout.close()
        f.wait()


    true_de_novo_s = set()
    true_s_sites = dict()
    if true_vcf and os.path.exists(true_vcf):
        f = subprocess_popen(shlex.split("gzip -fdc %s" % (true_vcf)))
        for row in f.stdout:
            if row[0] == '#':
                continue
            columns =row.rstrip().split('\t')
            ctg_name = columns[0]
            if contig_name and contig_name != ctg_name:
                continue
            pos = int(columns[1])
            ref_base = columns[3]
            alt_base = columns[4]
            qual = float(columns[5])

            if not (is_tree_empty or is_region_in(tree=tree, contig_name=ctg_name, \
                        region_start=pos - 1, region_end=pos)):
                continue

            def get_gt(x):
                tmp_s = x.split(':')[0]
                gt_1, gt_2 = int(tmp_s[0]), int(tmp_s[2])
                return gt_1 + gt_2
            try:
                gt_c, gt_p1, gt_p2 = get_gt(columns[c_col]), get_gt(columns[p1_col]), get_gt(columns[p2_col])
            except:
                print('error, please input merged trio VCF')
                return 0
            if gt_c == 1 and gt_p1 == 0 and gt_p2 == 0:
                true_de_novo_s.add(pos)
            else:
                continue
            true_s_sites[pos] = columns
            # print(row, gt_c, gt_p1, gt_p2)
            
        f.stdout.close()
        f.wait()
    called_n = len(call_de_novo_s)
    true_n = len(true_de_novo_s)
    TP_s = true_de_novo_s & call_de_novo_s
    TP_n = len(TP_s)
    FP_n = len(call_de_novo_s - TP_s)
    FN_n = len(true_de_novo_s - TP_s)
    # import pdb; pdb.set_trace()

    # print(true_de_novo_s & call_de_novo_s)
    print('in %s, called de novo sites %s' % (contig_name, called_n))
    print('in %s, true de novo sites %s' % (contig_name, true_n))
    print('FN, TP, FP: %s,%s,%s' % (FN_n, TP_n, FP_n))
    print('FN sites: ', true_de_novo_s - TP_s)
    # for i in true_de_novo_s - TP_s:
    #     print(true_s_sites[i])



def main():

    parser = ArgumentParser(description="Select heterozygous snp candidates, require input one chrosome a time")

    parser.add_argument('--call_vcf', type=str, default=None,
                        help="Path to the merged call VCf. (default: %(default)s)")

    parser.add_argument('--true_vcf', type=str, default=None,
                        help="Path to the true call VCf. (default: %(default)s)")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, default: %(default)s") 

    parser.add_argument('--bed_fn', type=str, default=None,
        help="constrain select region within bed file, (default: %(default)s)")


    args = parser.parse_args()
    Check(args)

if __name__ == "__main__":
    main()
