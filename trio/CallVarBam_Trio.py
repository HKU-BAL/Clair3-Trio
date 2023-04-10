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
from distutils.util import strtobool
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
        self.create_tensor = None
        self.call_variant = None

    def poll(self):
        self.create_tensor.poll()
        self.call_variant.poll()


ins_c = InstancesClass()


def check_return_code(signum, frame):
    ins_c.poll()
    if ins_c.create_tensor.returncode != None and ins_c.create_tensor.returncode != 0:
        ins_c.create_tensor.kill()
        sys.exit("CreateTensor.py exited with exceptions. Exiting...")

    if ins_c.call_variant.returncode != None and ins_c.call_variant.returncode != 0:
        ins_c.create_tensor.kill()
        sys.exit("call_variant.py exited with exceptions. Exiting...")

    if (ins_c.create_tensor.returncode == None or ins_c.call_variant.returncode == None):
        signal.alarm(5)


def Run(args):

    print('run callvarbam')
    basedir = dirname(__file__)

    CT_Bin = basedir + "/../clair3.py MergeTensorsBam_Trio"
    CVBin = basedir + "/../clair3.py CallVariants_Trio"


    if args.delay > 0:
        delay = random.randrange(0, args.delay)
        print("[INFO] Delay %d seconds before starting variant calling ..." % (delay))
        sleep(delay)

    pypyBin = executable_command_string_from(args.pypy, exit_on_not_found=True)
    pythonBin = executable_command_string_from(args.python, exit_on_not_found=True)
    samtoolsBin = executable_command_string_from(args.samtools, exit_on_not_found=True)


    bam_fn_c = file_path_from(args.bam_fn_c)
    bam_fn_p1 = file_path_from(args.bam_fn_p1)
    bam_fn_p2 = file_path_from(args.bam_fn_p2)
    if (bam_fn_c is None or bam_fn_c == "") or \
       (bam_fn_p1 is None or bam_fn_p2 == "") or \
       (bam_fn_p2 is None or bam_fn_p2 == ""):
        print(log_warning(
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
    sampleName_c = args.sampleName_c
    sampleName_p1 = args.sampleName_p1
    sampleName_p2 = args.sampleName_p2

    ctgName = args.ctgName
    min_af = args.min_af if args.min_af else param.min_af_dict[platform]
    # snp_min_af = args.snp_min_af
    # indel_min_af = args.indel_min_af

    if ctgName is None:
        sys.exit("--ctgName must be specified. You can call variants on multiple chromosomes simultaneously.")

    chkpnt_fn = args.chkpnt_fn
    qual = command_option_from(args.qual, 'qual', option_value=args.qual)
    temp_file_dir = command_option_from(args.temp_file_dir, 'temp_file_dir')
    haploid_precise_mode = CommandOption('haploid_precise', args.haploid_precise)
    haploid_sensitive_mode = CommandOption('haploid_sensitive', args.haploid_sensitive)
    # showRef_mode = command_option_from(None, None)
    # if args.showRef:
    #     showRef_mode = command_option_from(True, 'showRef')
    showRef_mode = CommandOption('showRef', args.showRef)
    add_indel_length_mode = CommandOption('add_indel_length', True)
    phasing_info_in_bam_mode = command_option_from(args.phasing_info_in_bam, 'phasing_info_in_bam')
    need_phasing_mode = command_option_from(args.need_phasing, 'need_phasing')
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

    create_tensor_command_options = [
        pypyBin,
        CT_Bin,
        CommandOption('bam_fn_c', bam_fn_c),
        CommandOption('bam_fn_p1', bam_fn_p1),
        CommandOption('bam_fn_p2', bam_fn_p2),
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



    sched_getaffinity_list = list(os.sched_getaffinity(0))
    maxCpus = len(sched_getaffinity_list)
    if args.tensorflow_threads is None:
        numCpus = maxCpus
    else:
        numCpus = args.tensorflow_threads if args.tensorflow_threads < maxCpus else maxCpus

    _cpuSet = ",".join(str(x) for x in random.sample(sched_getaffinity_list, numCpus))

    taskSet = "taskset -c %s" % (_cpuSet)
    try:
        subprocess.check_output("which %s" % ("taskset"), shell=True)
    except:
        taskSet = ""

    call_variant_command_options = [
        taskSet,
        pythonBin,
        CVBin,
        CommandOption('chkpnt_fn', chkpnt_fn),
        CommandOption('call_fn_c', args.call_fn_c),
        CommandOption('call_fn_p1', args.call_fn_p1),
        CommandOption('call_fn_p2', args.call_fn_p2),
        CommandOption('sampleName_c', args.sampleName_c),
        CommandOption('sampleName_p1', args.sampleName_p1),
        CommandOption('sampleName_p2', args.sampleName_p2),
        CommandOption('ref_fn', ref_fn),
        CommandOption('platform', platform),
        CommandOption('ctgName', ctgName),
        CommandOption('temp_file_dir', args.temp_file_dir),
        qual,
        add_indel_length_mode,
        showRef_mode,
        chunk_id,
        chunk_num,
        gvcf_mode,
    ]
    #print(' '.join(shlex.split(command_string_from(create_tensor_command_options))))
    #print(' '.join(shlex.split(command_string_from(call_variant_command_options))))
    #import pdb; pdb.set_trace()

    try:
        ins_c.create_tensor= subprocess_popen(
                shlex.split(command_string_from(create_tensor_command_options)), )

        #import pdb; pdb.set_trace()
        ins_c.call_variant = subprocess_popen(
            shlex.split(command_string_from(call_variant_command_options)),
            stdin=ins_c.create_tensor.stdout, stdout=sys.stderr)
    except Exception as e:
        print(e, file=sys.stderr)
        sys.exit("Failed to start required processes. Exiting...")

    #import pdb; pdb.set_trace()

    signal.signal(signal.SIGALRM, check_return_code)
    signal.alarm(2)

    try:
        ins_c.call_variant.wait()
        ins_c.create_tensor.stdout.close()
        ins_c.create_tensor.wait()
    except KeyboardInterrupt as e:
        print("KeyboardInterrupt received when waiting at CallVarBam, terminating all scripts.")
        try:
            ins_c.call_variant.terminate()
            ins_c.create_tensor.terminate()
        except Exception as e:
            print(e)
        raise KeyboardInterrupt
    except Exception as e:
        print("Exception received when waiting at CallVarBam, terminating all scripts.")
        print(e)
        try:
            ins_c.call_variant.terminate()
            ins_cc.create_tensor.terminate()
        except Exception as e:
            print(e)

        raise e


def main():
    # create tensors options
    parser = ArgumentParser(description="Generate variant candidate tensors using phased full-alignment")

    # legacy
    parser.add_argument('--platform', type=str, default='ont',
                        help="Sequencing platform of the input. Options: 'ont', default: %(default)s")

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

    # call variants options
    parser.add_argument('--chkpnt_fn', type=str, default=None, required=True,
                        help="Input a trained model for variant calling, required")

    parser.add_argument('--sampleName_c', type=str, default="Child",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--sampleName_p1', type=str, default="Parent1",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--sampleName_p2', type=str, default="Parent2",
                        help="Define the sample name to be shown in the VCF file, optional")

    parser.add_argument('--call_fn_c', type=str, default="clair3_c.vcf",
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--call_fn_p1', type=str, default="clair3_p1.vcf",
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--call_fn_p2', type=str, default="clair3_p2.vcf",
                        help="VCF output filename, or stdout if not set")

    parser.add_argument('--temp_file_dir', type=str, default='./',
                        help="EXPERIMENTAL: The cache directory for storing temporary non-variant information if --gvcf is enabled, default: %(default)s")

    parser.add_argument('--haploid_precise', action='store_true',
                        help="EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant")

    parser.add_argument('--haploid_sensitive', action='store_true',
                        help="EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant")

    # options for debug purpose
    parser.add_argument('--use_gpu', type=str2bool, default=False,
                        help="DEBUG: Use GPU for calling. Speed up is mostly insignficiant. Only use this for building your own pipeline")
   
    # output options
    ## Output reference calls
    # parser.add_argument('--showRef', type=str2bool, default=False)
    parser.add_argument('--showRef', dest='showRef', 
                    type=lambda x: bool(strtobool(x)))

    parser.add_argument('--qual', type=int, default=2,
                        help="If set, variants with >=$qual will be marked 'PASS', or 'LowQual' otherwise, optional")
    

    args = parser.parse_args()
    



    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    Run(args)


if __name__ == "__main__":
    main()
