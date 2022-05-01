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

def variant_from(vcf_fn, tree, is_tree_empty, contig_name = None):
    Y = set()

    if vcf_fn and os.path.exists(vcf_fn):
        f = subprocess_popen(shlex.split("gzip -fdc %s" % (vcf_fn)))
        for row in f.stdout:
            if row[0] == '#':
                continue
            columns = row.rstrip().split('\t')
            ctg_name = columns[0]
            if contig_name and contig_name != ctg_name:
                continue
            pos = int(columns[1])
            ref_base = columns[3]
            alt_base = columns[4]

            if not (is_tree_empty or is_region_in(tree, ctg_name, pos)):
                continue

            Y.add(pos)
        f.stdout.close()
        f.wait()
    return Y



def candidate_from(alt_fn, tree, is_tree_empty, contig_name = None):
    ref_call_pos_list = []
    var_call_pos_list = []

    if alt_fn and os.path.exists(alt_fn):
        f = subprocess_popen(shlex.split("gzip -fdc %s" % (alt_fn)))
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

            #ref_call
            if ref_base == alt_base  or alt_base == '.':
                ref_call_pos_list.append((pos,qual))
            else:
                var_call_pos_list.append((pos,qual))
        f.stdout.close()
        f.wait()
    return ref_call_pos_list, var_call_pos_list


def get_low_qual_rec(candidate, ref_pct_full, var_pct_full):
    ref_call_pos_list, var_call_pos_list = candidate[0], candidate[1]
    low_qual_ref_list = sorted(ref_call_pos_list, key=lambda x: x[1])[:int(ref_pct_full * len(ref_call_pos_list))]
    low_qual_variant_list = sorted(var_call_pos_list, key=lambda x: x[1])[:int(var_pct_full * len(var_call_pos_list))]
    return low_qual_ref_list, low_qual_variant_list

def FiterHeteSnp(args):

    """
    Filter heterozygous snp variant for training, if there are too many candidates for full alignment training, we
    would select more in low quality variants, which is more challenging for pileup model to predict and using more
    information will benefit calling those variants.
    """
    
    var_pct_full = args.var_pct_full
    ref_var_max_ratio = args.ref_var_max_ratio
    ref_pct_full = args.ref_pct_full if args.ref_pct_full is not None else var_pct_full
    chr_prefix = args.chr_prefix
    contig_name = args.ctgName
    phasing_window_size = param.phasing_window_size
    chunk_id = args.chunk_id - 1 if args.chunk_id else None # 1-base to 0-base
    DEPTH = args.depth
    chunk_num = args.chunk_num
    sample_name = args.sampleName
    split_bed_size = args.split_bed_size
    split_folder = args.split_folder
    extend_bp = param.extend_bp
    phasing_info_in_bam = args.phasing_info_in_bam
    need_phasing_list = []
    need_phasing_set = set()
    ref_call_pos_list = []
    chr_prefix_length = len(chr_prefix)
    variant_dict = defaultdict(str)
    realign_window_size = args.realign_window_size if args.realign_window_size is not None else param.flankingBaseNum
    candidate_positions = set()

    bed_fn = args.bed_fn
    tree = bed_tree_from(bed_fn, contig_name=contig_name)
    is_tree_empty = len(tree.keys()) == 0

    # legacy
    # vcf_fn_c=args.vcf_fn_c
    # vcf_fn_p1=args.vcf_fn_p1
    # vcf_fn_p2=args.vcf_fn_p2

    alt_fn_c = args.alt_fn_c
    alt_fn_p1 = args.alt_fn_p1
    alt_fn_p2 = args.alt_fn_p2


    # # legacy
    # vcf_fn = args.vcf_fn_c

    # legacy
    alt_fn = args.alt_fn_c

    # import pdb; pdb.set_trace()

    # Y_c = variant_from(vcf_fn_c, tree, is_tree_empty, contig_name)
    # Y_p1 = variant_from(vcf_fn_p1, tree, is_tree_empty, contig_name)
    # Y_p2 = variant_from(vcf_fn_p2, tree, is_tree_empty, contig_name) 
    # Y_all = Y_c | Y_p1 | Y_p2
    # print('in %s, readed true %d, %d, %d, union %d' % (contig_name, len(Y_c), len(Y_p1), len(Y_p2), len(Y_all)))

    candidate_c = candidate_from(alt_fn_c, tree, is_tree_empty, contig_name)
    candidate_p1 = candidate_from(alt_fn_p1, tree, is_tree_empty, contig_name)
    candidate_p2 = candidate_from(alt_fn_p2, tree, is_tree_empty, contig_name)

    # print(len(candidate_c[0]), len(candidate_c[1]))
    # print(len(candidate_p1[0]), len(candidate_p1[1]))
    # print(len(candidate_p2[0]), len(candidate_p2[1]))

    lq_can_c = get_low_qual_rec(candidate_c, ref_pct_full, var_pct_full)
    lq_can_p1 = get_low_qual_rec(candidate_p1, ref_pct_full, var_pct_full)
    lq_can_p2 = get_low_qual_rec(candidate_p2, ref_pct_full, var_pct_full)

    low_qual_ref_list, low_qual_variant_list = lq_can_c[0], lq_can_c[1]
    print('DEPTH %s' % (args.depth))
    print('[ORI C] {} {} select ref calling (cutoff {}) to process: {}'.format(sample_name, contig_name, low_qual_ref_list[-1][1], len(low_qual_ref_list)))
    print('[ORI C] {} {} select variant calling (cutoff {}) to process: {}'.format(sample_name, contig_name, low_qual_variant_list[-1][1], len(low_qual_variant_list)))
    
    low_qual_ref_list, low_qual_variant_list = lq_can_p1[0], lq_can_p1[1]
    print('[ORI P1] {} {} select ref calling (cutoff {}) to process: {}'.format(sample_name, contig_name, low_qual_ref_list[-1][1], len(low_qual_ref_list)))
    print('[ORI P1] {} {} select variant calling (cutoff {}) to process: {}'.format(sample_name, contig_name, low_qual_variant_list[-1][1], len(low_qual_variant_list)))
    
    low_qual_ref_list, low_qual_variant_list = lq_can_p2[0], lq_can_p2[1]
    print('[ORI P2] {} {} select ref calling (cutoff {}) to process: {}'.format(sample_name, contig_name, low_qual_ref_list[-1][1], len(low_qual_ref_list)))
    print('[ORI P2] {} {} select variant calling (cutoff {}) to process: {}'.format(sample_name, contig_name, low_qual_variant_list[-1][1], len(low_qual_variant_list)))
    
    can_set_c = set([item[0] for item in lq_can_c[1]])
    can_set_p1 = set([item[0] for item in lq_can_p1[1]])
    can_set_p2 = set([item[0] for item in lq_can_p2[1]])

    can_set_all = can_set_c | can_set_p1 | can_set_p2

    _ori_ref_n = int(sum([len(lq_can_c[0]), len(lq_can_p1[0]), len(lq_can_p2[0])]) / 3)
    _ori_v_n = int(sum([len(lq_can_c[1]), len(lq_can_p1[1]), len(lq_can_p2[1])]) / 3)


    ref_set_c = set([item[0] for item in lq_can_c[0]])
    ref_set_p1 = set([item[0] for item in lq_can_p1[0]])
    ref_set_p2 = set([item[0] for item in lq_can_p2[0]])

    ref_set_all = ref_set_c | ref_set_p1 | ref_set_p2
    n_ref_set_size = len(ref_set_all)

    _r_n_ref_set_size_ind = int(_ori_v_n * ref_var_max_ratio)

    _tmp_ref_set_c = set([item[0] for item in lq_can_c[0][:_r_n_ref_set_size_ind]])
    _tmp_ref_set_p1 = set([item[0] for item in lq_can_p1[0][:_r_n_ref_set_size_ind]])
    _tmp_ref_set_p2 = set([item[0] for item in lq_can_p2[0][:_r_n_ref_set_size_ind]])
    _tmp_ref_set_all = _tmp_ref_set_c | _tmp_ref_set_p1 | _tmp_ref_set_p2
    _r_n_ref_set_size = len(_tmp_ref_set_all)


    print('check ref site, %s or maximum %s (ind %s)' % (len(ref_set_all), _r_n_ref_set_size, _r_n_ref_set_size_ind))
    if _r_n_ref_set_size < len(ref_set_all):
        print('using maximum ref/var ratio %s' % ref_var_max_ratio) 
        n_ref_set_size = _r_n_ref_set_size
    # n_ref_set_size = min(int(_ori_ref_n / _ori_v_n * len(can_set_all)), len(ref_set_all))
    # n_ref_set_size = _ori_ref_n

    print("ori avg v: %d, ref: %d, new union: v %d, ref %d (%.2f%% /%d)" % (_ori_v_n, _ori_ref_n, len(can_set_all), \
                n_ref_set_size, n_ref_set_size * 100 / len(ref_set_all), len(ref_set_all)))

    tmp_ref_l = list(ref_set_all)
    random.seed(0)
    random.shuffle(tmp_ref_l)
    n_ref_set_all = set(tmp_ref_l[:n_ref_set_size])

    if _r_n_ref_set_size < len(ref_set_all):
        n_ref_set_all = set(_tmp_ref_set_all)

    # import pdb; pdb.set_trace()

    need_phasing_row_list = sorted(list(can_set_all | n_ref_set_all))

    if args.get_subtract:
        tmp_s_v = set(can_set_p1 | can_set_p2) - can_set_c
        tmp_s_r = set(ref_set_p1 | ref_set_p2) - ref_set_c
        print('get trio subtract candidate v: {}, non v: {}'.format(len(tmp_s_v), len(tmp_s_r)))
        need_phasing_row_list = sorted(list(tmp_s_v | tmp_s_r))
        # import pdb; pdb.set_trace()

    print('total site %d' % (len(need_phasing_row_list)))
    all_candidate_size = len(need_phasing_row_list)
    chunk_size = all_candidate_size // chunk_num + 1 if all_candidate_size % chunk_num else all_candidate_size // chunk_num

    for chunk_idx in range(chunk_num):
        start_pos = chunk_idx * chunk_size
        end_pos = min(start_pos + chunk_size, all_candidate_size)
        split_output = need_phasing_row_list[start_pos:end_pos]
        split_output = [(item - realign_window_size, item + realign_window_size + 2) for item in
                        split_output]  # a windows region for create tensor # samtools mpileup not include last position

        split_output = sorted(split_output, key=lambda x: x[0])
        with open(os.path.join(split_folder,
            '{}_{}_{}_{}'.format(sample_name, DEPTH, contig_name[chr_prefix_length:], chunk_idx+1)), # zero-base to one-base
            'w') as output_file:
            output_file.write('\n'.join(
                ['\t'.join([contig_name, str(x[0] - 1), str(x[1] - 1), ]) for x in
                 split_output]) + '\n')  # bed format
       

def main():

    parser = ArgumentParser(description="Select heterozygous snp candidates, require input one chrosome a time")

    parser.add_argument('--split_folder', type=str, default=None,
                        help="Path to directory that stores small bed region for raw alignment. (default: %(default)s)")

    # legacy
    parser.add_argument('--alt_fn', type=str, default=None, help=SUPPRESS)

    ## Path of provided alternative file
    parser.add_argument('--alt_fn_c', type=str, default=None, help=SUPPRESS)
    parser.add_argument('--alt_fn_p1', type=str, default=None, help=SUPPRESS)
    parser.add_argument('--alt_fn_p2', type=str, default=None, help=SUPPRESS)

    parser.add_argument('--var_pct_full', type=float, default=0.3,
                        help="Default variant call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--ref_pct_full', type=float, default=None,
                        help="Default reference call proportion for raw alignment or remove low quality proportion for whatshap phasing. (default: %(default)f)")

    parser.add_argument('--ref_var_max_ratio', type=float, default=5,
                        help="Default variant call proportion for ref sites. (default: %(default)f)")


    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, default: %(default)s")

    parser.add_argument('--sampleName', type=str, default="",
                        help="Define the sample name to be shown in the VCF file, optional")

    # options for debug purpose
    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="DEBUG: Input bam or sam have phasing info in HP tag, default: False")

    parser.add_argument('--split_bed_size', type=int, default=1000,
                        help="DEBUG: Default split bed size for parallel excution, default: %(default)s")

    parser.add_argument('--realign_window_size', type=int, default=None,
                        help="DEBUG: The window size of read realignment, work with need_realignment")

    parser.add_argument('--split_region_size', type=int, default=40000000,
                        help="DEBUG: Vcf phasing split_region_size default: %(default)s")

    # options for internal process control
    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None, help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None, help=SUPPRESS)

    ## Output all alternative candidates path
    parser.add_argument('--all_alt_fn', type=str, default=None, help=SUPPRESS)

    ## Default chr prefix for contig name
    parser.add_argument('--chr_prefix', type=str, default='chr', help=SUPPRESS)

    ## Default subsample depth for subsample bam file, 1000 means no subsampling
    parser.add_argument('--depth', type=int, default=1000, help=SUPPRESS)

    parser.add_argument('--bed_fn', type=str, default=None,
        help="constrain select region within bed file, (default: %(default)s)")

    parser.add_argument('--get_subtract', type=str2bool, default=False,
                        help=SUPPRESS)

    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    FiterHeteSnp(args)

if __name__ == "__main__":
    main()
