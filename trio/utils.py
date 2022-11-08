import sys
import gc
import copy
import shlex
import os
import tables
import numpy as np
from random import random
from collections import Counter

from clair3.task.main import *
from shared.interval_tree import bed_tree_from, is_region_in
from shared.utils import subprocess_popen, IUPAC_base_to_ACGT_base_dict as BASE2BASE, IUPAC_base_to_num_dict as BASE2NUM


from clair3.task.gt21 import GT21_Type
from clair3.task.variant_length import variant_length_index_offset


FILTERS = tables.Filters(complib='blosc:lz4hc', complevel=5)
shuffle_bin_size = 50000
output_bin_size = 1000
# shuffle_bin_size = 100
PREFIX_CHAR_STR = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"


def setup_environment():
    gc.enable()


def batches_from(iterable, item_from, batch_size=1):
    iterable = iter(iterable)
    while True:
        chunk = []
        for _ in range(batch_size):
            try:
                chunk.append(item_from(next(iterable)))
            except StopIteration:
                yield chunk
                return
        yield chunk


def tensor_generator_from(tensor_file_path, batch_size, pileup, platform):
    global param
    float_type = 'int32'
    import trio.param_t as param
    float_type = 'int8'

    if tensor_file_path != "PIPE":
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tensor_file_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    processed_tensors = 0

    #legacy
    tensor_shape = param.ont_input_shape_trio

    tensor_shape_one = param.ont_input_shape
    tensor_shape_trio = param.ont_input_shape_trio

    prod_tensor_shape = np.prod(tensor_shape_trio)

    def item_from(row):
        chrom, coord, seq, string_c, alt_info_c, string_p1, alt_info_p1, string_p2, alt_info_p2 = row.rstrip().split("\t")

        new_position_matrix_c = padded_tensor(string_c, tensor_shape_one)
        new_position_matrix_p1 = padded_tensor(string_p1, tensor_shape_one)
        new_position_matrix_p2 = padded_tensor(string_p2, tensor_shape_one)
        tensor = new_position_matrix_c + new_position_matrix_p1 + new_position_matrix_p2
        tensor = np.array(tensor, dtype=np.dtype(float_type))

        pos = chrom + ":" + coord + ":" + seq
        alt_info = '%s\t%s\t%s' % (alt_info_c, alt_info_p1, alt_info_p2)
        return tensor, pos, seq, alt_info

    for batch in batches_from(fo, item_from=item_from, batch_size=batch_size):
        #import pdb; pdb.set_trace()
        tensors = np.empty(([batch_size, prod_tensor_shape]), dtype=np.dtype(float_type))
        positions = []
        alt_info_list = []
        for tensor, pos, seq, alt_info in batch:
            if seq[param.flankingBaseNum] not in BASE2NUM:
                continue
            tensors[len(positions)] = tensor
            positions.append(pos)
            alt_info_list.append(alt_info)

        current_batch_size = len(positions)
        X = np.reshape(tensors, ([batch_size] + tensor_shape_trio))

        if processed_tensors > 0 and processed_tensors % 20000 == 0:
            print("Processed %d tensors" % processed_tensors, file=sys.stderr)

        processed_tensors += current_batch_size

        if current_batch_size <= -1:
            continue
        yield X[:current_batch_size], positions[:current_batch_size], alt_info_list[:current_batch_size]

    if tensor_file_path != "PIPE":
        fo.close()
        f.wait()

def ori_tensor_generator_from(tensor_file_path, batch_size, pileup, platform):
    global param
    float_type = 'int32'
    if pileup:
        import shared.param_p as param
    else:
        import shared.param_f as param
        float_type = 'int8'

    if tensor_file_path != "PIPE":
        f = subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tensor_file_path)))
        fo = f.stdout
    else:
        fo = sys.stdin

    processed_tensors = 0
    tensor_shape = param.ont_input_shape if platform == 'ont' else param.input_shape
    prod_tensor_shape = np.prod(tensor_shape)

    def item_from(row):
        chrom, coord, seq, tensor, alt_info = row.split("\t")
        if pileup:
            tensor = np.array(tensor.split(), dtype=np.dtype(float_type))
            depth = int(alt_info.split('-', maxsplit=1)[0])
            max_depth = param.max_depth_dict[platform]
            # for extreme high coverage data, make sure we could have a truncated coverage
            if depth > 0 and depth > max_depth * 1.5:
                scale_factor = depth / max_depth
                tensor = tensor / scale_factor
        else:
            # need add padding if depth is lower than maximum depth.
            tensor = [int(item) for item in tensor.split()]
            tensor_depth = len(tensor) // tensor_shape[1] // tensor_shape[2]
            padding_depth = tensor_shape[0] - tensor_depth
            prefix_padding_depth = int(padding_depth / 2)
            suffix_padding_depth = padding_depth - int(padding_depth / 2)
            prefix_zero_padding = [0] * prefix_padding_depth * tensor_shape[1] * tensor_shape[2]
            suffix_zero_padding = [0] * suffix_padding_depth * tensor_shape[1] * tensor_shape[2]
            tensor = prefix_zero_padding + tensor + suffix_zero_padding
            tensor = np.array(tensor, dtype=np.dtype(float_type))

        pos = chrom + ":" + coord + ":" + seq
        return tensor, pos, seq, alt_info

    for batch in batches_from(fo, item_from=item_from, batch_size=batch_size):
        tensors = np.empty(([batch_size, prod_tensor_shape]), dtype=np.dtype(float_type))
        positions = []
        alt_info_list = []
        for tensor, pos, seq, alt_info in batch:
            if seq[param.flankingBaseNum] not in BASE2NUM:
                continue
            tensors[len(positions)] = tensor
            positions.append(pos)
            alt_info_list.append(alt_info)

        current_batch_size = len(positions)
        X = np.reshape(tensors, ([batch_size] + tensor_shape))

        if processed_tensors > 0 and processed_tensors % 20000 == 0:
            print("Processed %d tensors" % processed_tensors, file=sys.stderr)

        processed_tensors += current_batch_size

        if current_batch_size <= 0:
            continue
        yield X[:current_batch_size], positions[:current_batch_size], alt_info_list[:current_batch_size]

    if tensor_file_path != "PIPE":
        fo.close()
        f.wait()



def remove_common_suffix(ref_base, alt_base):
    min_length = min(len(ref_base) - 1, min([len(item) - 1 for item in alt_base]))  # keep at least one base
    prefix = ref_base[::-1]
    for string in alt_base:
        string = string[::-1]
        while string[:len(prefix)] != prefix and prefix:
            prefix = prefix[:len(prefix) - 1]
        if not prefix:
            break
    res_length = len(prefix)
    if res_length > min_length:
        return ref_base, alt_base
    return ref_base[:len(ref_base) - res_length], [item[:len(item) - res_length] for item in alt_base]

    return ref_base[-min_length], [item[-min_length] for item in alt_base]


def decode_alt(ref_base, alt_base):
    if ',' not in alt_base:
        return [ref_base], [alt_base]
    alt_base = alt_base.split(',')
    ref_base_list, alt_base_list = [], []
    for ab in alt_base:
        rb,ab = remove_common_suffix(ref_base, [ab])
        ref_base_list.append(rb)
        alt_base_list.append(ab[0])
    return ref_base_list, alt_base_list


ALL_Y_INFO = {}
def variant_map_from(var_fn, tree, is_tree_empty, check_mcv_id = 0):
    Y = {}
    truth_alt_dict = {}
    miss_variant_set = set()
    if var_fn is None:
        return Y, miss_variant_set, truth_alt_dict

    f = subprocess_popen(shlex.split("gzip -fdc %s" % (var_fn)))
    for row in f.stdout:
        if row[0] == "#":
            continue
        columns = row.strip().split()
        # ctg_name, position_str, ref_base = columns[0], columns[1], columns[2]
        # genotype1, genotype2 = columns[-2], columns[-1]
        ctg_name, position_str, ref_base, alt_base, genotype1, genotype2 = columns
        key = ctg_name + ":" + position_str
        if genotype1 == '-1' or genotype2 == '-1':
            miss_variant_set.add(key)
            continue
        if not (is_tree_empty or is_region_in(tree, ctg_name, int(position_str))):
            continue

        Y[key] = output_labels_from_vcf_columns(columns)
        ref_base_list, alt_base_list = decode_alt(ref_base, alt_base)
        truth_alt_dict[int(position_str)] = (ref_base_list, alt_base_list)

        if check_mcv_id != 0:
            if key not in ALL_Y_INFO:
                ALL_Y_INFO[key] = dict()
            ALL_Y_INFO[key][check_mcv_id] = columns

    f.stdout.close()
    f.wait()
    # return Y, miss_variant_set
    return Y, miss_variant_set, truth_alt_dict

def padded_tensor(string, tensor_shape, add_padding=False, padding_value='100'):
    if len(string) == 1:
        string = string[0]
    position_matrix = string
    position_matrix = position_matrix.split()

    tensor_depth = len(position_matrix) // tensor_shape[1] // tensor_shape[2]
    padding_depth = tensor_shape[0] - tensor_depth
    prefix_padding_depth = int(padding_depth / 2)
    suffix_padding_depth = padding_depth - int(padding_depth / 2)

    if add_padding:
        # adding a new channel, 1 for no padding, 0 for padding
        new_position_matrix = []
        i = 0
        while i < len(position_matrix):
            new_position_matrix += position_matrix[i: (i + tensor_shape[2])]
            new_position_matrix += [padding_value]
            i += tensor_shape[2]
        prefix_zero_padding = ['0'] * prefix_padding_depth * tensor_shape[1] * (tensor_shape[2] + 1)
        suffix_zero_padding = ['0'] * suffix_padding_depth * tensor_shape[1] * (tensor_shape[2] + 1)
        new_position_matrix = prefix_zero_padding + new_position_matrix + suffix_zero_padding

    else:
        prefix_zero_padding = ['0'] * prefix_padding_depth * tensor_shape[1] * tensor_shape[2]
        suffix_zero_padding = ['0'] * suffix_padding_depth * tensor_shape[1] * tensor_shape[2]
        new_position_matrix = prefix_zero_padding + position_matrix + suffix_zero_padding

    return new_position_matrix


def write_table_dict(table_dict, pos, \
                string_c, alt_info_c, label_c, \
                string_p1, alt_info_p1, label_p1, \
                string_p2, alt_info_p2, label_p2, \
                total, tensor_shape, add_padding=False):
    """
    Write pileup or full alignment tensor into a dictionary.compressed bin file.
    table_dict: dictionary include all training information (tensor position, label, altnative bases).
    string: input tensor string, need add padding to meet the depth requirement.
    label: include gt21 genotype, indel length 1, indel length 2.
    alt_info: altnative information for querying variant.
    """
    global param

    # legacy
    string, alt_info, label = string_c, alt_info_c, label_c

    new_position_matrix_c = padded_tensor(string_c, tensor_shape, add_padding, param.padding_value_c)
    new_position_matrix_p1 = padded_tensor(string_p1, tensor_shape, add_padding, param.padding_value_p1)
    new_position_matrix_p2 = padded_tensor(string_p2, tensor_shape, add_padding, param.padding_value_p2)

    table_dict['position_matrix'].append(new_position_matrix_c + new_position_matrix_p1 + new_position_matrix_p2)
    table_dict['label'].append(label_c + label_p1 + label_p2)
    table_dict['alt_info'].append('%s\t%s\t%s' % (alt_info_c, alt_info_p1, alt_info_p2))

    if len(table_dict['alt_info'][-1]) > param.max_alt_info_length:
        print('[W], in pos %s, alt len overflow  %d < %d, %s' % (pos, param.max_alt_info_length, len(table_dict['alt_info'][-1]), table_dict['alt_info'][-1][:20]))

    table_dict['position'].append(pos)

    # print(alt_info_c, len(table_dict['position_matrix'][-1]), table_dict['alt_info'])
    # import pdb; pdb.set_trace()
    
    return total + 1


def update_table_dict():
    table_dict = {}
    table_dict['position_matrix'] = []
    table_dict['alt_info'] = []
    table_dict['position'] = []
    table_dict['label'] = []
    return table_dict


def write_table_file(table_file, table_dict, tensor_shape, label_size, float_type):
    """
    Write pileup or full alignment tensor into compressed bin file.
    table_dict: dictionary include all training information (tensor position, label, altnative bases).
    string: input tensor string, need add padding to meet the depth requirement.
    tree: dictionary(contig name : intervaltree) for quick region querying.
    miss_variant_set:  sometimes there will have true variant missing after downsampling reads.
    is_allow_duplicate_chr_pos: whether allow duplicate positions when training, if there exists downsampled data, lower depth will add a random prefix character.
    non_variant_subsample_ratio: define a maximum non variant ratio for training, we always expect use more non variant data, while it would greatly increase training
    time, especially in ont data, here we usually use 1:1 or 1:2 for variant candidate: non variant candidate.
    """

    position_matrix = np.array(table_dict['position_matrix'], np.dtype(float_type)).reshape([-1] + tensor_shape)
    table_file.root.position_matrix.append(position_matrix)

    table_file.root.alt_info.append(np.array(table_dict['alt_info']).reshape(-1, 1))
    table_file.root.position.append(np.array(table_dict['position']).reshape(-1, 1))
    table_file.root.label.append(np.array(table_dict['label'], np.dtype(float_type)).reshape(-1, label_size))
    
    # print(len(table_dict['alt_info']), tensor_shape, label_size)
    # import pdb; pdb.set_trace() 
    table_dict = update_table_dict()
    return table_dict


def print_bin_size(path, prefix=None):
    import tables
    import os
    total = 0
    for file_name in os.listdir(path):
        if prefix and not file_name.startswith(prefix):
            continue
        table = tables.open_file(os.path.join(path, file_name), 'r')
        print("[INFO] {} size is: {}".format(file_name, len(table.root.label)))
        total += len(table.root.label)
    print('[INFO] total: {}'.format(total))

  
def find_read_support(pos, truth_alt_dict, alt_info):
    alt_info = alt_info.rstrip().split('-')
    seqs = alt_info[1].split(' ') if len(alt_info) > 1 else ''
    seq_alt_bases_dict = dict(zip(seqs[::2], [int(item) for item in seqs[1::2]])) if len(seqs) else {}

    pos = int(pos)
    if pos not in truth_alt_dict:
        # candidate position not in the truth vcf or unified truth vcf
        return None
    ref_base_list, alt_base_list = truth_alt_dict[pos]
    found = 0
    for alt_type in seq_alt_bases_dict:
        if '*' in alt_type or '#' in alt_type or 'R' in alt_type:
            continue
        if alt_type[0] == 'X':
            if alt_type[1] in alt_base_list:
                found += 1
        elif alt_type[0] == 'I':
            if alt_type[1:] in alt_base_list:
                found += 1
        elif alt_type[0] == 'D':
            del_cigar = alt_type[1:]
            for rb, ab in zip(ref_base_list, alt_base_list):
                if rb[1:] == del_cigar and len(ab) == 1:
                    found += 1
    if found >= len(alt_base_list):
        return True
    # return False if we find any alternative bases missed in subsampled bam, then remove the position from training
    return False                                                   

def bin_reader_generator_from(subprocess_list, Y_true_var_c, Y_true_var_p1, Y_true_var_p2, Y_c, Y_p1, Y_p2, \
    is_tree_empty, tree, miss_variant_set_c, miss_variant_set_p1, miss_variant_set_p2, \
    truth_alt_dict_c, truth_alt_dict_p1, truth_alt_dict_p2, \
    is_allow_duplicate_chr_pos=False, maximum_non_variant_ratio=1.0):

    """
    Bin reader generator for bin file generation.
    subprocess_list: a list includes all tensor generator of each tensor file.
    Y (c, p1, p2): dictionary (contig name: label information) to store all variant and non variant information.
    tree: dictionary(contig name : intervaltree) for quick region querying.
    miss_variant_set (c, p1, p2):  sometimes there will have true variant missing after downsampling reads.
    is_allow_duplicate_chr_pos: whether allow duplicate positions when training, if there exists downsampled data, lower depth will add a random prefix character.
    maximum_non_variant_ratio: define a maximum non variant ratio for training, we always expect use more non variant data, while it would greatly increase training
    time, especially in ont data, here we usually use 1:1 or 1:2 for variant candidate: non variant candidate.
    """

    # legacy
    Y = Y_c
    miss_variant_set = miss_variant_set_c

    X = {}
    total = 0
    find_rst = []
    ref_list = []
    variant_set_with_read_support = set()
    variants_without_read_support = 0
    for f in subprocess_list:
        for row_idx, row in enumerate(f.stdout):
            chrom, coord, seq, string_c, alt_info_c, string_p1, alt_info_p1, string_p2, alt_info_p2 = row.rstrip().split("\t")
            # print(coord, alt_info_c, alt_info_p1, alt_info_p2)
            
            # legacy
            string, alt_info = string_c, alt_info_c

            if not (is_tree_empty or is_region_in(tree, chrom, int(coord))):
                continue

            seq = seq.upper()
            if seq[param.flankingBaseNum] not in 'ACGT':
                continue

            key = chrom + ":" + coord

            if key in miss_variant_set_c:
                continue

            is_reference = (key not in Y_true_var_c) and (key not in Y_true_var_p1) and (key not in Y_true_var_p2)

            have_read_support_c = find_read_support(pos=coord, truth_alt_dict=truth_alt_dict_c, alt_info=alt_info_c)
            have_read_support_p1 = find_read_support(pos=coord, truth_alt_dict=truth_alt_dict_p1, alt_info=alt_info_p1)
            have_read_support_p2 = find_read_support(pos=coord, truth_alt_dict=truth_alt_dict_p2, alt_info=alt_info_p2)

            find_rst.append((have_read_support_c, have_read_support_p1, have_read_support_p2))
           
            if have_read_support_c is not None and not have_read_support_c:
                miss_variant_set_c.add(key)
                variants_without_read_support += 1
                continue 

            if have_read_support_p1 is not None and not have_read_support_p1:
                miss_variant_set_p1.add(key)
                variants_without_read_support += 1
                continue 

            if have_read_support_p2 is not None and not have_read_support_p2:
                miss_variant_set_p2.add(key)
                variants_without_read_support += 1
                continue 
            variant_set_with_read_support.add(key)


            if key not in X:
                X[key] = (string_c, alt_info_c, string_p1, alt_info_p1, string_p2, alt_info_p2, seq)
                if is_reference:
                    # import pdb; pdb.set_trace()
                    ref_list.append(key)
            elif is_allow_duplicate_chr_pos:
                new_key = ""
                for character in PREFIX_CHAR_STR:
                    tmp_key = character + key
                    if tmp_key not in X:
                        new_key = tmp_key
                        break
                if len(new_key) > 0:
                    X[new_key] = (string_c, alt_info_c, string_p1, alt_info_p1, string_p2, alt_info_p2, seq)
                if is_reference:
                    # import pdb; pdb.set_trace()
                    ref_list.append(new_key)

            # update all Y info if is reference, Y my be del for not allow duplicate_chr_pos reason !!!
            if key not in Y_c:
                Y_c[key] = output_labels_from_reference(BASE2BASE[seq[param.flankingBaseNum]])
            if key not in Y_p1:
                Y_p1[key] = output_labels_from_reference(BASE2BASE[seq[param.flankingBaseNum]])
            if key not in Y_p2:
                Y_p2[key] = output_labels_from_reference(BASE2BASE[seq[param.flankingBaseNum]])

            if len(X) == shuffle_bin_size:
                if maximum_non_variant_ratio is not None:
                    _filter_non_variants(X, ref_list, maximum_non_variant_ratio)
                yield X, total
                X = {}
                ref_list = []
            total += 1
            if total % 100000 == 0:
                print("[INFO] Processed %d tensors" % total, file=sys.stderr)
        f.stdout.close()
        f.wait()


    _cnt = Counter(find_rst)
    _all_none_cnt = sum([_cnt[_i] for _i in _cnt if (_i == (None, None, None))])
    _false_cnt = sum([_cnt[_i] for _i in _cnt if False in _i])
    _true_cnt = len(find_rst) - _all_none_cnt - _false_cnt
    print('[INFO] read all sites, None, True, False hit', len(find_rst), _all_none_cnt, _true_cnt, _false_cnt)
    for key, value in _cnt.most_common():
        print('[INFO] ', key, value) 
    print("[INFO] Variants with read support/variants without read support: {}/{}".format(len(variant_set_with_read_support), variants_without_read_support))
    # import pdb; pdb.set_trace()
    if maximum_non_variant_ratio is not None:
        _filter_non_variants(X, ref_list, maximum_non_variant_ratio)
    yield X, total
    yield None, total


def _filter_non_variants(X, ref_list, maximum_non_variant_ratio):
    non_variant_num = len(ref_list)
    variant_num = len(X) - non_variant_num
    if non_variant_num > variant_num * maximum_non_variant_ratio:
        non_variant_keep_fraction = maximum_non_variant_ratio * variant_num / (1. * non_variant_num)
        probabilities = np.random.random_sample((non_variant_num,))
        for key, p in zip(ref_list, probabilities):
            if p > non_variant_keep_fraction:
                X.pop(key)


def get_training_array(tensor_fn, var_fn_c, var_fn_p1, var_fn_p2, bed_fn, bin_fn, shuffle=True, is_allow_duplicate_chr_pos=True, chunk_id=None,
                       chunk_num=None, platform='ont', pileup=False, add_padding=False, check_mcv=False, maximum_non_variant_ratio=None, candidate_details_fn_prefix=None):

    """
    Generate training array for training. here pytables with blosc:lz4hc are used for extreme fast compression and decompression,
    which can meet the requirement of gpu utilization. lz4hc decompression allows speed up training array decompression 4~5x compared
    with tensorflow tfrecord file format, current gpu utilization could reach over 85% with only 10G memory.
    tensor_fn: string format tensor acquired from CreateTensorPileup or CreateTensorFullAlign, include contig name position, tensor matrix, alternative information.
    var_fn (_c, _p1, _p2): simplified variant(vcf) format from GetTruths, which include contig name, position, reference base, alternative base, genotype.
    bin_fn: pytables format output bin file name.
    shuffle: whether apply index shuffling when generating training data, default True, which would promote robustness.
    is_allow_duplicate_chr_pos: whether allow duplicate positions when training, if there exists downsampled data, lower depth will add a random prefix character.
    chunk_id: specific chunk id works with total chunk_num for parallel execution. Here will merge all tensor file with sampe prefix.
    chunk_num: total chunk number for parallel execution. Each chunk refer to a smaller reference regions.
    platform: platform for tensor shape, ont give a larger maximum depth compared with pb and illumina.
    pileup: whether in pileup mode. Define two calling mode, pileup or full alignment.
    maximum_non_variant_ratio: define a maximum non variant ratio for training, we always expect use more non variant data, while it would greatly increase training
    time, especially in ont data, here we usually use 1:1 or 1:2 for variant candidate: non variant candidate.
    candidate_details_fn_prefix: a counter to calculate total variant and non variant from the information in alternative file.
    """
    # legacy
    var_fn = var_fn_c

    # read bed file
    tree = bed_tree_from(bed_file_path=bed_fn)
    is_tree_empty = len(tree.keys()) == 0

    # read true Y info
    Y_true_var_c, miss_variant_set_c, truth_alt_dict_c = variant_map_from(var_fn_c, tree, is_tree_empty, 'c' if check_mcv else 0)
    Y_true_var_p1, miss_variant_set_p1, truth_alt_dict_p1 = variant_map_from(var_fn_p1, tree, is_tree_empty, 'p1' if check_mcv else 0)
    Y_true_var_p2, miss_variant_set_p2, truth_alt_dict_p2 = variant_map_from(var_fn_p2, tree, is_tree_empty, 'p2' if check_mcv else 0)
    Y_c = copy.deepcopy(Y_true_var_c)
    Y_p1 = copy.deepcopy(Y_true_var_p1)
    Y_p2 = copy.deepcopy(Y_true_var_p2)
    # import pdb; pdb.set_trace()


    print('read true ', len(Y_c), len(Y_p1), len(Y_p2))
    # print('read true ', len(miss_variant_set_c), len(miss_variant_set_p1), len(miss_variant_set_p2))

    # legacy
    # Y, miss_variant_set = variant_map_from(var_fn, tree, is_tree_empty)

    global param
    import trio.param_t as param

    float_type = 'int8' #iinfo(min=-128, max=127, dtype=int8)
    tensor_shape = param.ont_input_shape_trio 

    # legacy
    tensor_shape_one = param.ont_input_shape 
    tensor_shape_trio = param.ont_input_shape_trio 

    if add_padding:
        # change trio shape when using add padding channel
        print("[INFO] add padding tensors %s, %s, %s" % \
            (param.padding_value_c, param.padding_value_p1, param.padding_value_p2), file=sys.stderr)
        tensor_shape_trio = param.p_ont_input_shape_trio 

    # import pdb; pdb.set_trace()

    label_size_one = param.label_size
    label_size_trio = param.label_size_trio

    if maximum_non_variant_ratio != None:
        print("[INFO] non variants/ variants subsample ratio set to: {}".format(maximum_non_variant_ratio))

    
    # get all tensors files name and split into different chunck set
    # select all match prefix if file path not exists
    subprocess_list = []
    if os.path.exists(tensor_fn):
        subprocess_list.append(subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, tensor_fn))))
    else:
        tensor_fn = tensor_fn.split('/')
        directry, file_prefix = '/'.join(tensor_fn[:-1]), tensor_fn[-1]
        all_file_name = []
        for file_name in os.listdir(directry):
            if file_name.startswith(file_prefix + '_') or file_name.startswith(
                    file_prefix + '.'):  # add '_.' to avoid add other prefix chr
                all_file_name.append(file_name)
        all_file_name = sorted(all_file_name)
        if chunk_id is not None:
            chunk_size = len(all_file_name) // chunk_num if len(all_file_name) % chunk_num == 0 else len(
                all_file_name) // chunk_num + 1
            chunk_start = chunk_size * chunk_id
            chunk_end = chunk_start + chunk_size
            all_file_name = all_file_name[chunk_start:chunk_end]
        if not len(all_file_name):
            print("[INFO] chunk_id exceed total file number, skip chunk", file=sys.stderr)
            return 0
        for file_name in all_file_name:
            subprocess_list.append(subprocess_popen(shlex.split("{} -fdc {}".format(param.zstd, os.path.join(directry, file_name)))))


    tables.set_blosc_max_threads(64)
    int_atom = tables.Atom.from_dtype(np.dtype(float_type))
    string_atom = tables.StringAtom(itemsize=param.no_of_positions + 50)
    long_string_atom = tables.StringAtom(itemsize=param.max_alt_info_length)  # max alt_info length
    table_file = tables.open_file(bin_fn, mode='w', filters=FILTERS)
    table_file.create_earray(where='/', name='position_matrix', atom=int_atom, shape=[0] + tensor_shape_trio,
                             filters=FILTERS)
    table_file.create_earray(where='/', name='position', atom=string_atom, shape=(0, 1), filters=FILTERS)
    table_file.create_earray(where='/', name='label', atom=int_atom, shape=(0, label_size_trio), filters=FILTERS)
    table_file.create_earray(where='/', name='alt_info', atom=long_string_atom, shape=(0, 1), filters=FILTERS)

    table_dict = update_table_dict()
    # import pdb; pdb.set_trace()

    # generator to avoid high memory occupy
    bin_reader_generator = bin_reader_generator_from(subprocess_list=subprocess_list,
                                                     Y_true_var_c=Y_true_var_c,
                                                     Y_true_var_p1=Y_true_var_p1,
                                                     Y_true_var_p2=Y_true_var_p2,
                                                     Y_c=Y_c,
                                                     Y_p1=Y_p1,
                                                     Y_p2=Y_p2,
                                                     is_tree_empty=is_tree_empty,
                                                     tree=tree,
                                                     miss_variant_set_c=miss_variant_set_c,
                                                     miss_variant_set_p1=miss_variant_set_p1,
                                                     miss_variant_set_p2=miss_variant_set_p2,
                                                     truth_alt_dict_c=truth_alt_dict_c,
                                                     truth_alt_dict_p1=truth_alt_dict_p1,
                                                     truth_alt_dict_p2=truth_alt_dict_p2,
                                                     is_allow_duplicate_chr_pos=is_allow_duplicate_chr_pos,
                                                     maximum_non_variant_ratio=maximum_non_variant_ratio)

    def get_label(key, Y, is_allow_duplicate_chr_pos):
        label, new_key = None, ''
        if key in Y:
            label = Y[key]
            new_key = key
            if not is_allow_duplicate_chr_pos:
                del Y[key]
        elif is_allow_duplicate_chr_pos:
            tmp_key = key[1:]
            label = Y[tmp_key]
            new_key = tmp_key
            pos = tmp_key + ':' + seq
        return label, new_key

        
    cnt_mvc, cnt_all, cnt_pass = 0, 0, 0
    total_compressed = 0
    while True:
        X, total = next(bin_reader_generator)

        if X is None or not len(X):
            break
        all_chr_pos = sorted(X.keys())
        if shuffle == True:
            np.random.shuffle(all_chr_pos)
        for key in all_chr_pos:
            string_c, alt_info_c, string_p1, alt_info_p1, string_p2, alt_info_p2, seq = X[key]

            #legacy
            string, alt_info = string_c, alt_info_c

            del X[key]

            # get label
            label_c, _key = get_label(key, Y_c, is_allow_duplicate_chr_pos)
            label_p1, _ = get_label(key, Y_p1, is_allow_duplicate_chr_pos)
            label_p2, _ = get_label(key, Y_p2, is_allow_duplicate_chr_pos)

            pos = _key + ':' + seq

            # legacy
            label = label_c

            if (label_c is None) or (label_p1 is None) or (label_p2 is None):
                print('warning, key empty', key)
                continue
            cnt_all += 1

            # check mcv

            def get_info_from_label(label):
                idx = np.argmax(label[:21])
                gt_21 = GT21_LABELS[idx]
                indel_1 = np.argmax(label[24:57]) - variant_length_index_offset # 21 + 3 : 21 + 3 + 33
                indel_2 = np.argmax(label[57:]) - variant_length_index_offset # 21 + 3 + 33 :
                return gt_21, indel_1, indel_2

            def get_ind_gt(_ref, _alt):
                if len(_ref) == len(_alt):
                    return _alt[0]
                if len(_ref) > len(_alt):
                    return 'D-' + _ref[1: 1 + len(_ref) - len(_alt)]
                if len(_ref) < len(_alt):
                    return 'I-' + _alt[1: 1 + len(_alt) - len(_ref)]

            def get_sep_gt(gt_21, indel_1, indel_2, _key, _t, ALL_Y_INFO):
                ar1, ar2 = '', ''
                if len(gt_21) == 2:
                    ar1, ar2 = gt_21[0], gt_21[1]
                    return ar1, ar2
                if _key in ALL_Y_INFO and _t in ALL_Y_INFO[_key]:
                    _ref_s, _alt_s, _gt1, _gt2 = \
                    ALL_Y_INFO[_key][_t][2], ALL_Y_INFO[_key][_t][3].split(','), int(ALL_Y_INFO[_key][_t][4]), int(ALL_Y_INFO[_key][_t][5])
                    _alt_s = [_ref_s] + _alt_s
                    ar1 = get_ind_gt(_ref_s, _alt_s[_gt1])
                    ar2 = get_ind_gt(_ref_s, _alt_s[_gt2])
                    return ar1, ar2
                # if is indel but no alter infor
                print('err', gt_21, indel_1, indel_2, _key)


            if check_mcv:
                
                c_gt_21, c_indel_1, c_indel_2 = get_info_from_label(label_c)
                p1_gt_21, p1_indel_1, p1_indel_2 = get_info_from_label(label_p1)
                p2_gt_21, p2_indel_1, p2_indel_2 = get_info_from_label(label_p2)

                c_ar1, c_ar2 = get_sep_gt(c_gt_21, c_indel_1, c_indel_2, _key, 'c', ALL_Y_INFO)
                p1_ar1, p1_ar2 = get_sep_gt(p1_gt_21, p1_indel_1, p1_indel_2, _key, 'p1', ALL_Y_INFO)
                p2_ar1, p2_ar2 = get_sep_gt(p2_gt_21, p2_indel_1, p2_indel_2, _key, 'p2', ALL_Y_INFO)

                _MC = 0
                for _ar1 in [p1_ar1, p1_ar2]:
                    for _ar2 in [p2_ar1, p2_ar2]:
                        if ((c_ar1 == _ar1) and (c_ar2 == _ar2)) or ((c_ar2 == _ar1) and (c_ar1 == _ar2)):
                            _MC = 1
                            break
                    if _MC == 1:
                        break

                if _key in ALL_Y_INFO and _MC==0:
                    print(_key)
                    print(ALL_Y_INFO[_key])
                    print(c_gt_21, c_indel_1, c_indel_2, c_ar1, c_ar2)
                    print(p1_gt_21, p1_indel_1, p1_indel_2, p1_ar1, p1_ar2)
                    print(p2_gt_21, p2_indel_1, p2_indel_2, p2_ar1, p2_ar2)
                    print(_MC)
                    cnt_mvc += 1
                    # import pdb; pdb.set_trace()
                    continue

            # continue
            cnt_pass += 1

            total_compressed = write_table_dict(table_dict, pos, \
                string_c, alt_info_c, label_c, \
                string_p1, alt_info_p1, label_p1, \
                string_p2, alt_info_p2, label_p2, \
                total_compressed, tensor_shape_one, add_padding)

            # import pdb; pdb.set_trace()
            # write to table on every output_bin_size records parsed
            if (total_compressed % output_bin_size == 0) and total_compressed > 0:
                # print(total_compressed)
                # import pdb; pdb.set_trace()
                table_dict = write_table_file(table_file, table_dict, tensor_shape_trio, label_size_trio, float_type)

            if total_compressed % 50000 == 0:
                print("[INFO] Compressed %d tensor" % (total_compressed), file=sys.stderr)

    if (total_compressed % output_bin_size != 0) and total_compressed > 0:
        table_dict = write_table_file(table_file, table_dict, tensor_shape_trio, label_size_trio, float_type)

    table_file.close()
    print("[INFO] Compressed %d/%d tensor" % (total_compressed, total), file=sys.stderr)
    print("[INFO] total sites %d, mcv %d, pass %d" % (cnt_all, cnt_mvc, cnt_pass), file=sys.stderr)
