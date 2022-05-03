import sys
import gc
import shlex
import os
import tables
import numpy as np
from random import random
import logging
from argparse import ArgumentParser, SUPPRESS

from shared.utils import subprocess_popen
import trio.param_t as param
from subprocess import PIPE

logging.basicConfig(format='%(message)s', level=logging.INFO)

def get_all_candidate_sites(candidate_fn):
    sites = {}
    for row in open(candidate_fn, 'r'):
        _info=row.strip().split('\t')

        # note that some alt may have not alt info [XDL]
        if len(_info) == 2:
            _info.append('')
        if len(_info) != 3:
            continue

        chr_pos = _info[0]
        key = chr_pos.replace(' ', ':')
        _alt = _info[1] + ':' + _info[2]
        sites[key] = _alt
    return sites


def reader_generator_from_tensors(f_n):
    for row_idx, row in enumerate(f_n.stdout):
        chrom, coord, seq, string, alt_info = row.split("\t")
        alt_info = alt_info.rstrip()
        key = chrom+":"+str(coord)
        yield [key, chrom, int(coord), seq, alt_info, string]
    f_n.stdout.close()
    f_n.wait()
    yield None

# get all tensors in same position (add_no_phasing_data_training)
def _generator_get_lst(g):
    _key, _lst = None, []
    while True:
        _info = next(g)
        if _info is None or len(_info) <= 1:
            break
        if _key is not None and _info[0] != _key:
            yield _lst
            _lst = []
        _lst.append(_info)
        _key = _info[0]
    yield _lst
    yield []


def Merge_tensors(args):
    _trio_info_n = args.candidate_fn.strip().split('/')[-1]
    logging.info("%s: begin merging trio tensors" % (_trio_info_n))


    # import pdb; pdb.set_trace()


    #read candidates info
    can_c = get_all_candidate_sites(args.candidate_fn_c)
    can_p1 = get_all_candidate_sites(args.candidate_fn_p1)
    can_p2 = get_all_candidate_sites(args.candidate_fn_p2)

    s_c = set(can_c.keys())
    s_p1 = set(can_p1.keys())
    s_p2 = set(can_p2.keys())
    # import pdb; pdb.set_trace()

    trio_sites = s_c & s_p2 & s_p1
    non_trio_sites = s_c - trio_sites
    _trio_s, _c_s = len(trio_sites), len(s_c)
    logging.info("%s: trio sites [%d, %d, %d], valid %d/%d (%.2f%%)" % (_trio_info_n, len(can_c), \
        len(can_p1), len(can_p2), _trio_s, _c_s, 100.*_trio_s/_c_s))

    # import pdb; pdb.set_trace()


    subprocess_c = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, args.tensor_fn_c)))
    subprocess_p1 = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, args.tensor_fn_p1)))
    subprocess_p2 = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, args.tensor_fn_p2)))


    _generator_c = reader_generator_from_tensors(subprocess_c)
    _generator_p1 = reader_generator_from_tensors(subprocess_p1)
    _generator_p2 = reader_generator_from_tensors(subprocess_p2)

    lst_g_c = _generator_get_lst(_generator_c)
    lst_g_p1 = _generator_get_lst(_generator_p1)
    lst_g_p2 = _generator_get_lst(_generator_p2)

    tensor_can_fpo = open(args.tensor_fn, "wb")
    tensor_can_fp = subprocess_popen(shlex.split("{} -c".format(param.zstd)), stdin=PIPE, stdout=tensor_can_fpo)

    alt_fp = open(args.candidate_fn, 'w')


    # read and merge tensors
    _hit_cnt = 0
    _hit_len_cnt = {}
    while True:
        _info_c = next(lst_g_c)
        if len(_info_c) > 0:
            _tar_key = _info_c[0][0]
            # get tensor constrain to called at all sites
            # if _tar_key not in trio_sites:
            #     continue
            if _tar_key in non_trio_sites:
                continue

            _info_p1 = next(lst_g_p1)
            while _info_p1[0][0] != _tar_key:
                _info_p1 = next(lst_g_p1)

            _info_p2 = next(lst_g_p2)
            while _info_p2[0][0] != _tar_key:
                _info_p2 = next(lst_g_p2)

            _hit_l = (len(_info_c), len(_info_p1), len(_info_p2))
            if _hit_l[1] > _hit_l[2]:
                _hit_l = (_hit_l[0], _hit_l[2], _hit_l[1])
            if _hit_l not in _hit_len_cnt:
                _hit_len_cnt[_hit_l] = 1
            else:
                _hit_len_cnt[_hit_l] += 1
            # if len(_info_c) != len(_info_p1) or len(_info_c) != len(_info_p2):
            #     print('no match', _trio_info_n, _info_c[0][0], len(_info_c), len(_info_p1), len(_info_p2))
            #     # return 1
            #     # import pdb; pdb.set_trace()
            _hit_cnt += 1


            # output alt
            for _idx, _c_tensor in enumerate(_info_c):
                tensor_str = "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    _c_tensor[1], _c_tensor[2], _c_tensor[3], _c_tensor[5], _c_tensor[4], \
                    _info_p1[-1*_idx][5], _info_p1[-1*_idx][4], \
                    _info_p2[-1*_idx][5], _info_p2[-1*_idx][4]
                    )
                tensor_can_fp.stdin.write(tensor_str)

                # print(_idx, -1*_idx, _c_tensor[5][:20], ' | ', _info_p1[-1*_idx][5][:20], ' | ', _info_p2[-1*_idx][5][:20])
                # import pdb; pdb.set_trace()

            alt_info_c = can_c[_tar_key].replace(':', '\t')
            alt_info_p1 = can_p1[_tar_key].replace(':', '\t')
            alt_info_p2 = can_p2[_tar_key].replace(':', '\t')
            ctg_name, pos = _tar_key.split(':')
            alt_fp.write('\t'.join([ctg_name + ' ' + str(pos), alt_info_c, alt_info_p1, alt_info_p2]) + '\n')

            # import pdb; pdb.set_trace()
        else:
            break
    print('%s: hit trio %d' % (_trio_info_n, _hit_cnt))
    for _hit_l in _hit_len_cnt:
        print(_hit_l, _hit_len_cnt[_hit_l])
    # import pdb; pdb.set_trace()

    tensor_can_fp.stdin.close()
    tensor_can_fp.wait()
    tensor_can_fpo.close()

    alt_fp.close()

    logging.info("Finish!")



# combine trio tensors and convert to single tensor
def main():
    parser = ArgumentParser(description="Combine Trio tensors")

    parser.add_argument('--tensor_fn_c', type=str, default=None, required=True, help="Tensor input, required")
    parser.add_argument('--tensor_fn_p1', type=str, default=None, required=True, help="Tensor input, required")
    parser.add_argument('--tensor_fn_p2', type=str, default=None, required=True, help="Tensor input, required")

    parser.add_argument('--candidate_fn_c',type=str, default=None, required=True, help="Tensor alt info, required")
    parser.add_argument('--candidate_fn_p1', type=str, default=None, required=True, help="Tensor alt info, required")
    parser.add_argument('--candidate_fn_p2', type=str, default=None, required=True, help="Tensor alt info, required")

    parser.add_argument('--tensor_fn', type=str, default=None, required=True, help="output tensor name, required")
    parser.add_argument('--candidate_fn', type=str, default=None, required=True, help="output alt name, required")


    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Merge_tensors(args)


if __name__ == "__main__":
    main()
