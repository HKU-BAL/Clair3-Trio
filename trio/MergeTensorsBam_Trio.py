import sys
import shlex
import subprocess
import multiprocessing
import signal
import random
import os
from os.path import dirname
from time import sleep
from argparse import ArgumentParser, SUPPRESS
import logging

logging.getLogger().setLevel(logging.INFO)


from shared.command_options import (
    CommandOption,
    CommandOptionWithNoValue,
    ExecuteCommand,
    command_string_from,
    command_option_from
)
from shared.utils import file_path_from, executable_command_string_from, subprocess_popen, str2bool, log_warning
# import shared.param_p as param
import trio.param_t as param


class InstancesClass(object):
    def __init__(self):
        # legacy
        self.create_tensor = None

        self.create_tensor_c = None
        self.create_tensor_p1 = None
        self.create_tensor_p2 = None

    def poll(self):
        # legacy
        self.create_tensor.poll()

        self.create_tensor_c.poll()
        self.create_tensor_p1.poll()
        self.create_tensor_p2.poll()


ins_c = InstancesClass()


def check_return_code(signum, frame):
    ins_c.poll()
    if ins_c.create_tensor.returncode != None and ins_c.create_tensor.returncode != 0:
        ins_c.create_tensor.kill()
        sys.exit("CreateTensor.py exited with exceptions. Exiting...")

    if ins_c.create_tensor_c.returncode != None and ins_c.create_tensor_c.returncode != 0:
        ins_c.create_tensor_c.kill()
        sys.exit("CreateTensor.py exited with exceptions. Exiting...")

    if ins_c.create_tensor_p1.returncode != None and ins_c.create_tensor_p1.returncode != 0:
        ins_c.create_tensor_p1.kill()
        sys.exit("CreateTensor.py exited with exceptions. Exiting...")

    if ins_c.create_tensor_p2.returncode != None and ins_c.create_tensor_p2.returncode != 0:
        ins_c.create_tensor_p2.kill()
        sys.exit("CreateTensor.py exited with exceptions. Exiting...")

    # if (
    #         ins_c.create_tensor_c.returncode == None or
    #         ins_c.create_tensor_p1.returncode == None or
    #         ins_c.create_tensor_p2.returncode
    # ):
        # signal.alarm(5)

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



def get_info_k(_tar):
    _info = next(_tar)
    if _info != None and len(_info) > 0:
        return  _info, int(_info[0][0].split(':')[1])
    else:
        return None, None

def _get_min_k(k1, k2, k3):
    try:
        _cur_k = min(k1, k2, k3)
    except:
        _cur_k = None
    return _cur_k


class TensorStdout(object):
    def __init__(self, handle):
        self.stdin = handle

    def __del__(self):
        self.stdin.close()


def Run(args):
    logging.info('run MergeTensorsBAM_Trio')

    if args.tensor_can_fn != "PIPE":
        tensor_can_fpo = open(args.tensor_can_fn, "wb")
        tensor_can_fp = subprocess_popen(shlex.split("{} -c".format(args.zstd)), stdin=PIPE, stdout=tensor_can_fpo)
    else:
        tensor_can_fp = TensorStdout(sys.stdout)

    basedir = dirname(__file__)

    CT_Bin = basedir + "/../clair3.py CreateTensorFullAlignment"

    pypyBin = executable_command_string_from(args.pypy, exit_on_not_found=True)
    pythonBin = executable_command_string_from(args.python, exit_on_not_found=True)
    samtoolsBin = executable_command_string_from(args.samtools, exit_on_not_found=True)

    bam_fn_c = file_path_from(args.bam_fn_c)
    bam_fn_p1 = file_path_from(args.bam_fn_p1)
    bam_fn_p2 = file_path_from(args.bam_fn_p2)
    # import pdb; pdb.set_trace()
    if (bam_fn_c is None or bam_fn_c == "") or \
       (bam_fn_p1 is None or bam_fn_p2 == "") or \
       (bam_fn_p2 is None or bam_fn_p2 == ""):
        logging.info(log_warning(
            "[WARNING] Skip clair3-trio variant calling for empty bam file"))
        return

    ref_fn = file_path_from(args.ref_fn, exit_on_not_found=True)
    bed_fn = file_path_from(args.bed_fn)
    vcf_fn = file_path_from(args.vcf_fn)
    extend_bed = file_path_from(args.extend_bed)
    full_aln_regions = file_path_from(args.full_aln_regions)

    platform = args.platform
    if not platform or platform not in param.support_platform:
        sys.exit("[ERROR] Provided platform are not in support platform list [ont]")

    # call_fn = args.call_fn
    sampleName = args.sampleName
    ctgName = args.ctgName
    min_af = args.min_af if args.min_af else param.min_af_dict[platform]
    # snp_min_af = args.snp_min_af
    # indel_min_af = args.indel_min_af

    if ctgName is None:
        sys.exit("--ctgName must be specified. You can call variants on multiple chromosomes simultaneously.")

    # haploid_precise_mode = command_option_from(args.haploid_precise, 'haploid_precise')
    # haploid_sensitive_mode = command_option_from(args.haploid_sensitive, 'haploid_sensitive')
    # output_for_ensemble = command_option_from(args.output_for_ensemble, 'output_for_ensemble')

    showRef_mode = command_option_from(args.showRef, 'showRef')
    qual = command_option_from(args.qual, 'qual', option_value=args.qual)

    # add_indel_length_mode = CommandOption('add_indel_length', args.add_indel_length)
    add_indel_length_mode = True
    phasing_info_in_bam_mode = command_option_from(args.phasing_info_in_bam, 'phasing_info_in_bam')
    need_phasing_mode = command_option_from(args.need_phasing, 'need_phasing')
    # is_from_tables_mode = command_option_from(args.is_from_tables, 'is_from_tables')
    # pileup_mode = command_option_from(args.pileup, 'pileup')
    # fast_mode = CommandOption('fast_mode', args.fast_mode)
    # call_snp_only_mode = CommandOption('call_snp_only', args.call_snp_only)
    gvcf_mode = CommandOption('gvcf', args.gvcf)

    ctgStart = None
    ctgEnd = None
    chunk_id = None
    chunk_num = None
    if args.ctgStart is not None and args.ctgEnd is not None and int(args.ctgStart) <= int(args.ctgEnd):
        ctgStart = CommandOption('ctgStart', args.ctgStart)
        ctgEnd = CommandOption('ctgEnd', args.ctgEnd)

    if args.chunk_id is not None and args.chunk_num is not None and int(args.chunk_id) <= int(args.chunk_num):
        chunk_id = CommandOption('chunk_id', args.chunk_id)
        chunk_num = CommandOption('chunk_num', args.chunk_num)


    # sched_getaffinity_list = list(os.sched_getaffinity(0))
    # maxCpus = len(sched_getaffinity_list)
    # if args.tensorflow_threads is None:
    #     numCpus = maxCpus
    # else:
    #     numCpus = args.tensorflow_threads if args.tensorflow_threads < maxCpus else maxCpus

    # _cpuSet = ",".join(str(x) for x in random.sample(sched_getaffinity_list, numCpus))

    # taskSet = "taskset -c %s" % (_cpuSet)
    # try:
    #     subprocess.check_output("which %s" % ("taskset"), shell=True)
    # except:
    #     taskSet = ""

    create_tensor_command_options = [
        pypyBin,
        CT_Bin,
        CommandOption('ref_fn', ref_fn),
        CommandOption('vcf_fn', vcf_fn),
        CommandOption('ctgName', ctgName),
        CommandOption('min_af', min_af),
        CommandOption('platform', platform),
        CommandOption('samtools', samtoolsBin),
        CommandOption('bed_fn', bed_fn),
        CommandOption('extend_bed', extend_bed),
        CommandOption('sampleName', args.sampleName),
        ctgStart,
        ctgEnd,
        chunk_id,
        chunk_num,
        gvcf_mode,
        phasing_info_in_bam_mode,
        need_phasing_mode,
        CommandOption('full_aln_regions', full_aln_regions),
    ]

    create_tensor_command_options_c = \
        create_tensor_command_options[:] + [CommandOption('bam_fn', bam_fn_c)]
    create_tensor_command_options_p1 = \
        create_tensor_command_options[:] + [CommandOption('bam_fn', bam_fn_p1)]
    create_tensor_command_options_p2 = \
        create_tensor_command_options[:] + [CommandOption('bam_fn', bam_fn_p2)]

    #print(' '.join(shlex.split(command_string_from(create_tensor_command_options_c))))
    #print(' '.join(shlex.split(command_string_from(create_tensor_command_options_p1))))
    #print(' '.join(shlex.split(command_string_from(create_tensor_command_options_p2))))
    #import pdb; pdb.set_trace()

    try:
        ins_c.create_tensor_c = subprocess_popen(
                shlex.split(command_string_from(create_tensor_command_options_c)), )

        ins_c.create_tensor_p1 = subprocess_popen(
                shlex.split(command_string_from(create_tensor_command_options_p1)), )

        ins_c.create_tensor_p2 = subprocess_popen(
                shlex.split(command_string_from(create_tensor_command_options_p2)), )
        # import pdb; pdb.set_trace()

        # ins_c.call_variant = subprocess_popen(
        #     shlex.split(command_string_from(call_variant_command_options)),
        #     stdin=c.create_tensor.stdout, stdout=sys.stderr
        # )
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit("Failed to start required processes. Exiting...")

    # merge tensors
    _generator_c = reader_generator_from_tensors(ins_c.create_tensor_c)
    _generator_p1 = reader_generator_from_tensors(ins_c.create_tensor_p1)
    _generator_p2 = reader_generator_from_tensors(ins_c.create_tensor_p2)

    lst_g_c = _generator_get_lst(_generator_c)
    lst_g_p1 = _generator_get_lst(_generator_p1)
    lst_g_p2 = _generator_get_lst(_generator_p2)

    # init
    _info_c, _key_c = get_info_k(lst_g_c)
    _info_p1, _key_p1 = get_info_k(lst_g_p1)
    _info_p2, _key_p2 = get_info_k(lst_g_p2)
    _cur_k = _get_min_k(_key_c, _key_p1, _key_p2)

    _a_key_c, _a_key_p1, _a_key_p2 = [], [], []
    _a_key_c.append(_key_c)
    _a_key_p1.append(_key_p2)
    _a_key_p2.append(_key_p2)

    while _cur_k != None and _key_c != None and _key_p1 != None and _key_p2 != None:
        if _key_c == _key_p1 and _key_p1 == _key_p2:
            # yeild hit
            # output alt
            for _idx, _c_tensor in enumerate(_info_c):
                tensor_str = "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    _c_tensor[1], _c_tensor[2], _c_tensor[3], _c_tensor[5], _c_tensor[4], \
                    _info_p1[-1*_idx][5], _info_p1[-1*_idx][4], \
                    _info_p2[-1*_idx][5], _info_p2[-1*_idx][4]
                    )
                tensor_can_fp.stdin.write(tensor_str)
                #import pdb; pdb.set_trace()

            # update
            # logging.info('hit %s' % (_key_c))
            _info_c, _key_c = get_info_k(lst_g_c)
            _info_p1, _key_p1 = get_info_k(lst_g_p1)
            _info_p2, _key_p2 = get_info_k(lst_g_p2)

            _a_key_c.append(_key_c)
            _a_key_p1.append(_key_p2)
            _a_key_p2.append(_key_p2)

        while _cur_k != None and _key_c != None and _key_c <= _cur_k:
            _info_c, _key_c = get_info_k(lst_g_c)
            _a_key_c.append(_key_c)
            # logging.info('update c %s' % (_key_c))

        while _cur_k != None and _key_p1 != None and _key_p1 <= _cur_k:
            _info_p1, _key_p1 = get_info_k(lst_g_p1)
            _a_key_p1.append(_key_p1)
            # logging.info('update p1 %s' % (_key_p1))

        while _cur_k != None and _key_p2 != None and _key_p2 <= _cur_k:
            _info_p2, _key_p2 = get_info_k(lst_g_p2)
            _a_key_p2.append(_key_p2)
            # logging.info('update p2 %s' % (_key_p2))

        # update all key to have a minumum of _cur_k
        _cur_k = _get_min_k(_key_c, _key_p1, _key_p2)

 

    _a_key_c = [i for i in _a_key_c if i != None]
    _a_key_p1 = [i for i in _a_key_p1 if i != None]
    _a_key_p2 = [i for i in _a_key_p2 if i != None]
    s_c = set(_a_key_c)
    s_p1 = set(_a_key_p1)
    s_p2 = set(_a_key_p2)
    trio_sites = s_c & s_p2 & s_p1
    non_trio_sites = s_c - trio_sites
    _trio_s, _c_s = len(trio_sites), len(s_c)
    #import pdb; pdb.set_trace()
    if _c_s != 0:
        logging.info("%s trio sites [%d, %d, %d], valid %d/%d (%.2f%%)" % (full_aln_regions, len(s_c), \
            len(s_p1), len(s_p2), _trio_s, _c_s, 100.*_trio_s/_c_s))
    else:
        logging.info("%s trio sites [%d, %d, %d], valid %d/%d" % (full_aln_regions, len(s_c), \
            len(s_p1), len(s_p2), _trio_s, _c_s))


    # signal.signal(signal.SIGALRM, check_return_code)
    # signal.alarm(2)

    # tensor_can_fp.stdin.close()
    # tensor_can_fp.wait()
    # tensor_can_fpo.close()

    # closs fullalignment generation
    while _key_c != None:
        _info_c, _key_c = get_info_k(lst_g_c)

    while _key_p1 != None:
        _info_p1, _key_p1 = get_info_k(lst_g_p1)

    while _key_p2 != None:
        _info_p2, _key_p2 = get_info_k(lst_g_p2)

    logging.info("Finish!")

def main():
    # create tensors options
    parser = ArgumentParser(description="Generate variant candidate tensors using phased full-alignment")

    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    # legacy
    parser.add_argument('--bam_fn', type=str, default="input.bam", required=False,
                        help="Sorted BAM file input, required")

    parser.add_argument('--bam_fn_c', type=str, default="input.bam", required=True,
                        help="Sorted BAM file input, required")

    parser.add_argument('--bam_fn_p1', type=str, default="input.bam", required=True,
                        help="Sorted BAM file input, required")

    parser.add_argument('--bam_fn_p2', type=str, default="input.bam", required=True,
                        help="Sorted BAM file input, required")

    parser.add_argument('--ref_fn', type=str, default="ref.fa", required=True,
                        help="Reference fasta file input, required")

    parser.add_argument('--tensor_can_fn', type=str, default="PIPE",
                        help="Tensor output, stdout by default, default: %(default)s")

    parser.add_argument('--vcf_fn', type=str, default=None,
                        help="Candidate sites VCF file input, if provided, variants will only be called at the sites in the VCF file,  default: %(default)s")

    parser.add_argument('--min_af', type=float, default=0.08,
                        help="Minimum allele frequency for both SNP and Indel for a site to be considered as a condidate site, default: %(default)f")

    parser.add_argument('--ctgName', type=str, default=None,
                        help="The name of sequence to be processed, required if --bed_fn is not defined")

    parser.add_argument('--ctgStart', type=int, default=None,
                        help="The 1-based starting position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--ctgEnd', type=int, default=None,
                        help="The 1-based inclusive ending position of the sequence to be processed, optional, will process the whole --ctgName if not set")

    parser.add_argument('--bed_fn', type=str, default=None,
                        help="Call variant only in the provided regions. Will take an intersection if --ctgName and/or (--ctgStart, --ctgEnd) are set")

    parser.add_argument('--gvcf', type=str2bool, default=False,
                        help="Enable GVCF output, default: disabled")

    # legacy
    parser.add_argument('--sampleName', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the GVCF file")

    parser.add_argument('--sampleName_c', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the GVCF file")

    parser.add_argument('--sampleName_p1', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the GVCF file")

    parser.add_argument('--sampleName_p2', type=str, default="SAMPLE",
                        help="Define the sample name to be shown in the GVCF file")

    # options for advanced users
    parser.add_argument('--minCoverage', type=float, default=2,
                        help="EXPERIMENTAL: Minimum coverage required to call a variant, default: %(default)f")

    parser.add_argument('--minMQ', type=int, default=5,
                        help="EXPERIMENTAL: If set, reads with mapping quality with <$minMQ are filtered, default: %(default)d")

    parser.add_argument('--minBQ', type=int, default=0,
                        help="EXPERIMENTAL: If set, bases with base quality with <$minBQ are filtered, default: %(default)d")

    parser.add_argument('--max_depth', type=int, default=144,
                        help="EXPERIMENTAL: Maximum full alignment depth to be processed. default: %(default)s")

    # options for debug purpose
    parser.add_argument('--phasing_info_in_bam', action='store_true',
                        help="DEBUG: Skip phasing and use the phasing info provided in the input BAM (HP tag), default: False")

    parser.add_argument('--phasing_window_size', type=int, default=param.phasing_window_size,
                        help="DEBUG: The window size for read phasing")

    parser.add_argument('--extend_bed', nargs='?', action="store", type=str, default=None,
                        help="DEBUG: Extend the regions in the --bed_fn by a few bp for tensor creation, default extend 16bp")

    parser.add_argument('--indel_fn', type=str, default=None,
                        help="DEBUG: Output all alternative indel cigar for debug purpose")

    parser.add_argument('--base_err', default=0.001, type=float,
                        help='DEBUG: Estimated base error rate in gvcf option, default: %(default)f')

    parser.add_argument('--gq_bin_size', default=5, type=int,
                        help='DEBUG: Default gq bin size for merge non-variant block in gvcf option, default: %(default)d')

    parser.add_argument('--bp_resolution', action='store_true',
                        help="DEBUG: Enable bp resolution for GVCF, default: disabled")
    
    # options for internal process control
    ## Path to the 'zstd' compression
    parser.add_argument('--zstd', type=str, default=param.zstd,
                        help=SUPPRESS)

    ## Test in specific candidate position. Only for testing
    parser.add_argument('--test_pos', type=int, default=0,
                        help=SUPPRESS)

    ## The number of chucks to be divided into for parallel processing
    parser.add_argument('--chunk_num', type=int, default=None,
                        help=SUPPRESS)

    ## The chuck ID to work on
    parser.add_argument('--chunk_id', type=int, default=None,
                        help=SUPPRESS)

    ## Only call variant in phased vcf file
    parser.add_argument('--phased_vcf_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Apply no phased data in training. Only works in data training, default: False
    parser.add_argument('--add_no_phasing_data_training', action='store_true',
                        help=SUPPRESS)

    ## Output representation unification infos, which refines training labels
    parser.add_argument('--unify_repre', action='store_true',
                        help=SUPPRESS)

    ## Path of representation unification output
    parser.add_argument('--unify_repre_fn', type=str, default=None,
                        help=SUPPRESS)

    ## Provide the regions to be included in full-alignment based calling
    parser.add_argument('--full_aln_regions', type=str, default=None,
                        help=SUPPRESS)

    ## Use Clair3's own phasing module for read level phasing when creating tensor, compared to using Whatshap, speed is faster but has higher memory footprint, default: False
    parser.add_argument('--need_phasing', action='store_true',
                        help=SUPPRESS)

    ## Apply read realignment for illumina platform. Greatly boost indel performance in trade of running time
    parser.add_argument('--need_realignment', action='store_true',
                        help=SUPPRESS)




    parser.add_argument('--samtools', type=str, default="samtools",
                        help="Path to the 'samtools', samtools version >= 1.10 is required, default: %(default)s")

    parser.add_argument('--pypy', type=str, default="pypy3",
                        help="Path to the 'pypy', pypy3 version >= 3.6 is required, default: %(default)s")

    parser.add_argument('--python', type=str, default="python3",
                        help="Path to the 'python3', default: %(default)s")


    parser.add_argument('--tensorflow_threads', type=int, default=4,
                        help="DEBUG: Number of threads per tensorflow job. Tune if you are building your own pipeline")

    ## Wait a short while for no more than a few seconds to start the job. This is to avoid starting multiple jobs simultaneously
    ## that might use up the maximum number of threads allowed, because Tensorflow will create more threads than needed at the beginning of running the program
    ## Obseleted after adding --tensorflow_threads defaulted at 4
    parser.add_argument('--delay', type=int, default=5,
                        help=SUPPRESS)


    # output options
    ## Output reference calls
    parser.add_argument('--showRef', action='store_true',
                        help=SUPPRESS)

    parser.add_argument('--qual', type=int, default=2,
                        help="If set, variants with >=$qual will be marked 'PASS', or 'LowQual' otherwise, optional")
    
    args = parser.parse_args()
    



    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
