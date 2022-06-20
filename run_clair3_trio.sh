#!/bin/bash
SCRIPT_NAME=$(basename "$0")
SCRIPT_PATH=`dirname "$0"`
Usage="Usage: ./${SCRIPT_NAME} --bam_fn_c=BAM --bam_fn_p1=BAM --bam_fn_p2=BAM --ref_fn=REF --output=OUTPUT_DIR --threads=THREADS --model_path_clair3=MODEL_PREFIX --model_path_clair3_trio=MODEL_PREFIX [--bed_fn=BED] [options]"


# ENTRANCE SCRIPT FOR CLAIR3-TRIO, SETTING VARIABLE AND CALL TRIO
VERSION='v0.3'

set -e
#./run_clair3_trio.sh --bam_fn_c=child_bam --bam_fn_p1=parent1 --bam_fn_p2=parent2 -f ref.fasta -t 32 -o tmp -p --model_path_clair3=model_path --model_path_clair3_trio=model_path
print_help_messages()
{
    echo $''
    echo ${Usage}
    echo $''
    echo $'Required parameters:'
    echo $'--bam_fn_c=FILE                BAM file input, for child. The input file must be samtools indexed.'
    echo $'--bam_fn_p1=FILE               BAM file input, for parent1. The input file must be samtools indexed.'
    echo $'--bam_fn_p2=FILE               BAM file input, for parent1. The input file must be samtools indexed.'
    echo $'--ref_fn=FILE                  FASTA reference file input. The input file must be samtools indexed.'
    echo $'--model_path_clair3=STR        The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00002, pileup.data-00001-of-00002 pileup.index, full_alignment.data-00000-of-00002, full_alignment.data-00001-of-00002  and full_alignment.index).'
    echo $'--model_path_clair3_trio=STR   The folder path containing a Clair3-Trio model (files structure same as Clair3).'
    echo $'-t, --threads=INT              Max #threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The #chunks being processed simultaneously is ceil(#threads/4)*3. 3 is the overloading factor.'
    echo $'-o, --output=PATH              VCF/GVCF output directory.'
    echo $''
    echo $''
    echo $"Optional parameters (Use \"=value\" instead of \" value\". E.g., \"--bed_fn=fn.bed\" instead of \"--bed_fn fn.bed\".):"
    echo $'-v, --version                  Check Clair3-Trio version'
    echo $'-h, --help                     Check Clair3-Trio help page'
    echo $'--bed_fn=FILE                  Call variants only in the provided bed regions.'
    echo $'--vcf_fn=FILE                  Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.'
    echo $'--ctg_name=STR                 The name of the sequence to be processed.'
    echo $'--pileup_only                  Use the pileup model only when calling, default: disable.'
    echo $'--pileup_phasing               Use the pileup model calling and phasing, default: disable.'
    echo $'--sample_name_c=STR            Define the sample name for Child to be shown in the VCF file.[Child]'
    echo $'--sample_name_p1=STR           Define the sample name for Parent1 to be shown in the VCF file.[Parent1]'
    echo $'--sample_name_p2=STR           Define the sample name for Parent2 to be shown in the VCF file.[Parent2]'
    echo $'--qual=INT                     If set, variants with >=$qual will be marked PASS, or LowQual otherwise.'
    echo $'--samtools=STR                 Path of samtools, samtools version >= 1.10 is required.'
    echo $'--python=STR                   Path of python, python3 >= 3.6 is required.'
    echo $'--pypy=STR                     Path of pypy3, pypy3 >= 3.6 is required.'
    echo $'--parallel=STR                 Path of parallel, parallel >= 20191122 is required.'
    echo $'--whatshap=STR                 Path of whatshap, whatshap >= 1.0 is required.'
    echo $'--chunk_size=INT               The size of each chuck for parallel processing, default: 5000000.'
    echo $'--print_ref_calls              Show reference calls (0/0) in VCF file, default: disable.'
    echo $'--include_all_ctgs             Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}, default: disable.'
    echo $'--gvcf                         Enable GVCF output, default: disable.'
    echo $'--snp_min_af=FLOAT             Minimum SNP AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.08,hifi:0.08,ilmn:0.08.'
    echo $'--indel_min_af=FLOAT           Minimum Indel AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.15,hifi:0.08,ilmn:0.08.'
    echo $'--remove_intermediate_dir      Remove intermediate directory, including intermediate phased BAM, pileup and full-alignment results. default: disable.'
    echo $'--enable_phasing               Output phased variants using whatshap, default: disable.'
    echo $'--var_pct_full=FLOAT           EXPERIMENTAL: Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.3.'
    echo $'--ref_pct_full=FLOAT           EXPERIMENTAL: Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for ilmn and hifi, 0.1 for ont.'
    echo $'--var_pct_phasing=FLOAT        EXPERIMENTAL: Specify an expected percentage of high quality 0/1 variants used in WhatsHap phasing, default: 0.8 for ont guppy5 and 0.7 for other platforms.'
    echo $'--pileup_model_prefix=STR      EXPERIMENTAL: Model prefix in pileup calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index. default: pileup.'
    echo $'--trio_model_prefix=STR        EXPERIMENTAL: Model prefix in trio calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index, default: trio.'
    echo $'--fast_mode                    EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: disable.'
    echo $'--haploid_precise              EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant, default: disable.'
    echo $'--haploid_sensitive            EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant, default: disable.'
    echo $'--call_snp_only                EXPERIMENTAL: Call candidates pass SNP minimum AF only, ignore Indel candidates, default: disable.'
    echo $'--enable_long_indel            EXPERIMENTAL: Call long Indel variants(>50 bp), default: disable.'
    echo $'--no_phasing_for_fa            EXPERIMENTAL: Call variants without whatshap phasing in full alignment calling, default: disable.'
    echo $'-p, --platform=STR             EXPERIMENTAL: [DO NOT CHANGE] Select the sequencing platform of the input. Possible options: {ont}.'
    echo $''
}

print_version()
{
    echo "Clair3-Trio $1"
    exit 0
}

ERROR="\\033[31m[ERROR]"
WARNING="\\033[33m[WARNING]"
NC="\\033[0m"

ARGS=`getopt -o b:f:t:p:o:hv \
-l bam_fn_c:,bam_fn_p1:,bam_fn_p2:,ref_fn:,threads:,model_path_clair3:,model_path_clair3_trio:,platform:,output:,\
bed_fn::,vcf_fn::,ctg_name::,sample_name_c::,sample_name_p1::,sample_name_p2::,qual::,samtools::,python::,pypy::,parallel::,whatshap::,chunk_num::,chunk_size::,var_pct_full::,ref_pct_full::,var_pct_phasing::,\
resumn::,snp_min_af::,indel_min_af::,pileup_model_prefix::,trio_model_prefix::,fast_mode,gvcf,pileup_only,pileup_phasing,print_ref_calls,haploid_precise,haploid_sensitive,include_all_ctgs,no_phasing_for_fa,call_snp_only,remove_intermediate_dir,enable_phasing,enable_long_indel,help,version -n 'run_clair3_trio.sh' -- "$@"`

if [ $? != 0 ] ; then echo"No input. Terminating...">&2 ; exit 1 ; fi
eval set -- "${ARGS}"

# default options
SAMPLE_C="SAMPLE_C"
SAMPLE_P1="SAMPLE_P1"
SAMPLE_P2="SAMPLE_P2"
BED_FILE_PATH="EMPTY"
VCF_FILE_PATH='EMPTY'
CONTIGS="EMPTY"
SAMTOOLS="samtools"
PYPY="pypy3"
PYTHON='python3'
PARALLEL='parallel'
WHATSHAP='whatshap'
PLATFORM=ont
CHUNK_NUM=0
CHUNK_SIZE=5000000
QUAL=2
PHASING_PCT="0"
PRO=0.3
REF_PRO=0
GVCF=False
RESUMN=0
PILEUP_ONLY=False
PILEUP_PHASING=False
FAST_MODE=False
SHOW_REF=False
SNP_AF=0
INDEL_AF=0
HAP_PRE=False
HAP_SEN=False
SNP_ONLY=False
INCLUDE_ALL_CTGS=False
NO_PHASING=False
RM_TMP_DIR=False
ENABLE_PHASING=False
ENABLE_LONG_INDEL=False
PILEUP_PREFIX="pileup"
TRIO_PREFIX="trio"

while true; do
   case "$1" in
    --bam_fn_c ) BAM_FILE_PATH_C="$2"; shift 2 ;;
    --bam_fn_p1 ) BAM_FILE_PATH_P1="$2"; shift 2 ;;
    --bam_fn_p2 ) BAM_FILE_PATH_P2="$2"; shift 2 ;;
    -f|--ref_fn ) REFERENCE_FILE_PATH="$2"; shift 2 ;;
    -t|--threads ) THREADS="$2"; shift 2 ;;
    --model_path_clair3 ) MODEL_PATH_C3="$2"; shift 2 ;;
    --model_path_clair3_trio ) MODEL_PATH_C3T="$2"; shift 2 ;;
    -p|--platform ) PLATFORM="$2"; shift 2 ;;
    -o|--output ) OUTPUT_FOLDER="$2"; shift 2 ;;
    --bed_fn ) BED_FILE_PATH="$2"; shift 2 ;;
    --vcf_fn ) VCF_FILE_PATH="$2"; shift 2 ;;
    --ctg_name ) CONTIGS="$2"; shift 2 ;;
    --sample_name_c ) SAMPLE_C="$2"; shift 2 ;;
    --sample_name_p1 ) SAMPLE_P1="$2"; shift 2 ;;
    --sample_name_p2 ) SAMPLE_P2="$2"; shift 2 ;;
    --chunk_num ) CHUNK_NUM="$2"; shift 2 ;;
    --chunk_size ) CHUNK_SIZE="$2"; shift 2 ;;
    --qual ) QUAL="$2"; shift 2 ;;
    --samtools ) SAMTOOLS="$2"; shift 2 ;;
    --python ) PYTHON="$2"; shift 2 ;;
    --pypy ) PYPY="$2"; shift 2 ;;
    --parallel ) PARALLEL="$2"; shift 2 ;;
    --whatshap ) WHATSHAP="$2"; shift 2 ;;
    --var_pct_full ) PRO="$2"; shift 2 ;;
    --ref_pct_full ) REF_PRO="$2"; shift 2 ;;
    --var_pct_phasing ) PHASING_PCT="$2"; shift 2 ;;
    --snp_min_af ) SNP_AF="$2"; shift 2 ;;
    --indel_min_af ) INDEL_AF="$2"; shift 2 ;;
    --pileup_model_prefix ) PILEUP_PREFIX="$2"; shift 2 ;;
    --trio_model_prefix ) TRIO_PREFIX="$2"; shift 2 ;;
    --gvcf ) GVCF=True; shift 1 ;;
    --resumn ) RESUMN="$2"; shift 2 ;;
    --pileup_only ) PILEUP_ONLY=True; shift 1 ;;
    --pileup_phasing ) PILEUP_PHASING=True; shift 1 ;;
    --fast_mode ) FAST_MODE=True; shift 1 ;;
    --call_snp_only ) SNP_ONLY=True; shift 1 ;;
    --print_ref_calls ) SHOW_REF=True; shift 1 ;;
    --haploid_precise ) HAP_PRE=True; shift 1 ;;
    --haploid_sensitive ) HAP_SEN=True; shift 1 ;;
    --include_all_ctgs ) INCLUDE_ALL_CTGS=True; shift 1 ;;
    --no_phasing_for_fa ) NO_PHASING=True; shift 1 ;;
    --remove_intermediate_dir ) RM_TMP_DIR=True; shift 1 ;;
    --enable_long_indel ) ENABLE_LONG_INDEL=True; shift 1 ;;
    --enable_phasing ) ENABLE_PHASING=True; shift 1 ;;

    -- ) shift; break; ;;
    -h|--help ) print_help_messages; exit 0 ;;
    -v|--version ) print_version ${VERSION}; exit 0 ;;
    * ) print_help_messages; break ;;
   esac
done

if [ -z ${BAM_FILE_PATH_C} ] || [ -z ${BAM_FILE_PATH_P1} ] || [ -z ${BAM_FILE_PATH_P2} ] || [ -z ${REFERENCE_FILE_PATH} ] || [ -z ${THREADS} ] || [ -z ${OUTPUT_FOLDER} ] ||  [ -z ${MODEL_PATH_C3} ] ||[ -z ${MODEL_PATH_C3T} ]; then
      if [ -z ${BAM_FILE_PATH_C} ] && [ -z ${BAM_FILE_PATH_P1} ] && [ -z ${BAM_FILE_PATH_P2} ] && [ -z ${REFERENCE_FILE_PATH} ] && [ -z ${THREADS} ] && [ -z ${OUTPUT_FOLDER} ] && [ -z ${MODEL_PATH_C3} ] && [ -z ${MODEL_PATH_C3T} ]; then print_help_messages; exit 0; fi
      if [ -z ${BAM_FILE_PATH_C} ]; then echo -e "${ERROR} Require to define index BAM input from child by --bam_fn_c=BAM${NC}"; fi
      if [ -z ${BAM_FILE_PATH_P1} ]; then echo -e "${ERROR} Require to define index BAM input from parent1 by --bam_fn_p1=BAM${NC}"; fi
      if [ -z ${BAM_FILE_PATH_P2} ]; then echo -e "${ERROR} Require to define index BAM input from parent2 by --bam_fn_p2=BAM${NC}"; fi
      if [ -z ${REFERENCE_FILE_PATH} ]; then echo -e "${ERROR} Require to define FASTA reference file input by --ref_fn=REF${NC}"; fi
      if [ -z ${THREADS} ]; then echo -e "${ERROR} Require to define max threads to be used by --threads=THREADS${NC}"; fi
      if [ -z ${OUTPUT_FOLDER} ]; then echo -e "${ERROR} Require to define output folder by --output=OUTPUT_DIR${NC}"; fi
      if [ -z ${MODEL_PATH_C3} ]; then echo -e "${ERROR} Require to define model path for cliar3 by --model_path_clair3=MODEL_PREFIX${NC}"; fi
      if [ -z ${MODEL_PATH_C3T} ]; then echo -e "${ERROR} Require to define model path for clair3-trio by --model_path_clair3_trio=MODEL_PREFIX${NC}"; fi
      exit 1;
fi

# force to use absolute path when in docker or singularity environment
if [ `pwd` = "/opt/bin" ]; then
    if [[ ! "${BAM_FILE_PATH_C}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --bam_fn_c=FILE${NC}"; exit 1; fi
    if [[ ! "${BAM_FILE_PATH_P1}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --bam_fn_p1=FILE${NC}"; exit 1; fi
    if [[ ! "${BAM_FILE_PATH_P2}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --bam_fn_p2=FILE${NC}"; exit 1; fi
    if [[ ! "${REFERENCE_FILE_PATH}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --ref_fn=FILE${NC}"; exit 1; fi
    if [[ ! "${MODEL_PATH_C3}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --model_path_clair3=PATH${NC}"; exit 1; fi
    if [[ ! "${MODEL_PATH_C3T}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --model_path_clair3_trio=PATH${NC}"; exit 1; fi
    if [[ ! "${OUTPUT_FOLDER}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --output=PATH${NC}"; exit 1; fi
    if [ "${BED_FILE_PATH}" != "EMPTY" ] &&  [ ! -z ${BED_FILE_PATH} ] && [[ ! "${BED_FILE_PATH}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --bef_fn=FILE${NC}"; exit 1; fi
    if [ "${VCF_FILE_PATH}" != "EMPTY" ] &&  [ ! -z ${VCF_FILE_PATH} ] && [[ ! "${VCF_FILE_PATH}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --vcf_fn=FILE${NC}"; exit 1; fi
fi

# relative path support
if [[ ! "${BAM_FILE_PATH_C}" = /* ]] && [ -f ${BAM_FILE_PATH_C} ]; then BAM_FILE_PATH_C=`pwd`/${BAM_FILE_PATH_C}; fi
if [[ ! "${BAM_FILE_PATH_P1}" = /* ]] && [ -f ${BAM_FILE_PATH_P1} ]; then BAM_FILE_PATH_P1=`pwd`/${BAM_FILE_PATH_P1}; fi
if [[ ! "${BAM_FILE_PATH_P2}" = /* ]] && [ -f ${BAM_FILE_PATH_P2} ]; then BAM_FILE_PATH_P2=`pwd`/${BAM_FILE_PATH_P2}; fi
if [[ ! "${REFERENCE_FILE_PATH}" = /* ]] && [ -f ${REFERENCE_FILE_PATH} ]; then REFERENCE_FILE_PATH=`pwd`/${REFERENCE_FILE_PATH}; fi
if [[ ! "${MODEL_PATH_C3}" = /* ]] && [ -d ${MODEL_PATH_C3} ]; then MODEL_PATH_C3=`pwd`/${MODEL_PATH_C3}; fi
if [[ ! "${MODEL_PATH_C3T}" = /* ]] && [ -d ${MODEL_PATH_C3T} ]; then MODEL_PATH_C3T=`pwd`/${MODEL_PATH_C3T}; fi
if [ "${BED_FILE_PATH}" != "EMPTY" ] && [ ! -z ${BED_FILE_PATH} ] && [[ ! "${BED_FILE_PATH}" = /* ]] && [ -f ${BED_FILE_PATH} ]; then BED_FILE_PATH=`pwd`/${BED_FILE_PATH}; fi
if [ "${VCF_FILE_PATH}" != "EMPTY" ] && [ ! -z ${VCF_FILE_PATH} ] && [[ ! "${VCF_FILE_PATH}" = /* ]] && [ -f ${VCF_FILE_PATH} ]; then VCF_FILE_PATH=`pwd`/${VCF_FILE_PATH}; fi
if [[ ! "${OUTPUT_FOLDER}" = /* ]]; then echo -e "${WARNING} No absolute output path provided, using current directory as prefix${NC}"; OUTPUT_FOLDER=`pwd`/${OUTPUT_FOLDER}; fi

mkdir -p ${OUTPUT_FOLDER}
if [ ! -d ${OUTPUT_FOLDER} ]; then echo -e "${ERROR} Cannot create output folder ${OUTPUT_FOLDER}${NC}"; exit 1; fi

# show default reference proportion 0.3 for ilmn and hifi, 0.1 for ont
if [ "${PLATFORM}" = "ont" ] && [ ! "${REF_PRO}" -gt 0 ]; then REF_PRO=0.1; fi
if [ "${PLATFORM}" != "ont" ] && [ ! "${REF_PRO}" -gt 0 ]; then REF_PRO=0.3; fi

# show default high quality hete variant proportion for whatshap phasing, 0.8 for ont guppy5 and 0.7 for others
if [ "${PHASING_PCT}" = "0" ]; then PHASING_PCT=0.7; fi
BASE_MODEL=$(basename ${MODEL_PATH_C3})
if [ "${BASE_MODEL}" = "r941_prom_sup_g5014" ] || [ "${BASE_MODEL}" = "r941_prom_hac_g5014" ] || [ "${BASE_MODEL}" = "ont_guppy5" ]; then PHASING_PCT=0.8; fi

# remove the last '/' character in directory input
OUTPUT_FOLDER=$(echo ${OUTPUT_FOLDER%*/})
MODEL_PATH_C3=$(echo ${MODEL_PATH_C3%*/})
MODEL_PATH_C3T=$(echo ${MODEL_PATH_C3T%*/})

# optional parameters should use "="
(time (
echo -e "\033[0;94mClair3-Trio\\033[0m (${VERSION})"
echo "-----"
echo "[INFO] OUTPUT FOLDER: ${OUTPUT_FOLDER}"
echo "[INFO] BAM_C FILE PATH:  ${BAM_FILE_PATH_C}"
echo "[INFO] BAM_P1 FILE PATH: ${BAM_FILE_PATH_P1}"
echo "[INFO] BAM_P2 FILE PATH: ${BAM_FILE_PATH_P2}"
echo "[INFO] SAMPLE NAME CHILD:   ${SAMPLE_C}"
echo "[INFO] SAMPLE NAME PARENT1: ${SAMPLE_P1}"
echo "[INFO] SAMPLE NAME PARENT2: ${SAMPLE_P2}"
echo "[INFO] REFERENCE FILE PATH: ${REFERENCE_FILE_PATH}"
echo "[INFO] CLAIR3 MODEL PATH: ${MODEL_PATH_C3}"
echo "[INFO] CLAIR3-TRIO MODEL PATH: ${MODEL_PATH_C3T}"
echo "[INFO] CLAIR3-TRIO MODEL PREFIX: ${TRIO_PREFIX}"
echo "[INFO] PLATFORM: ${PLATFORM}"
echo "[INFO] THREADS: ${THREADS}"
echo "[INFO] BED FILE PATH: ${BED_FILE_PATH}"
echo "[INFO] VCF FILE PATH: ${VCF_FILE_PATH}"
echo "[INFO] CONTIGS: ${CONTIGS}"
echo "[INFO] CONDA PREFIX: ${CONDA_PREFIX}"
echo "[INFO] SAMTOOLS PATH: ${SAMTOOLS}"
echo "[INFO] PYTHON PATH: ${PYTHON}"
echo "[INFO] PYPY PATH: ${PYPY}"
echo "[INFO] PARALLEL PATH: ${PARALLEL}"
echo "[INFO] WHATSHAP PATH: ${WHATSHAP}"
echo "[INFO] CHUNK SIZE: ${CHUNK_SIZE}"
if [ ${CHUNK_NUM} -gt 0 ]; then echo "[INFO] CHUNK NUM: ${CHUNK_NUM}"; fi
echo "[INFO] FULL ALIGN PROPORTION: ${PRO}"
echo "[INFO] FULL ALIGN REFERENCE PROPORTION: ${REF_PRO}"
echo "[INFO] PHASING PROPORTION: ${PHASING_PCT}"
if [ ${SNP_AF} -gt 0 ]; then echo "[INFO] USER DEFINED SNP THRESHOLD: ${SNP_AF}"; fi
if [ ${INDEL_AF} -gt 0 ]; then echo "[INFO] USER DEFINED INDEL THRESHOLD: ${INDEL_AF}"; fi
echo "[INFO] ENABLE FILEUP ONLY CALLING: ${PILEUP_ONLY}"
echo "[INFO] ENABLE PILEUP CALLING AND PHASING: ${PILEUP_PHASING}"
echo "[INFO] ENABLE FAST MODE CALLING: ${FAST_MODE}"
echo "[INFO] ENABLE CALLING SNP CANDIDATES ONLY: ${SNP_ONLY}"
echo "[INFO] ENABLE PRINTING REFERENCE CALLS: ${SHOW_REF}"
echo "[INFO] ENABLE OUTPUT GVCF: ${GVCF}"
echo "[INFO] ENABLE HAPLOID PRECISE MODE: ${HAP_PRE}"
echo "[INFO] ENABLE HAPLOID SENSITIVE MODE: ${GVCF}"
echo "[INFO] ENABLE INCLUDE ALL CTGS CALLING: ${INCLUDE_ALL_CTGS}"
echo "[INFO] ENABLE NO PHASING FOR FULL ALIGNMENT: ${NO_PHASING}"
echo "[INFO] ENABLE REMOVING INTERMEDIATE FILES: ${RM_TMP_DIR}"
echo "[INFO] ENABLE PHASING VCF OUTPUT: ${ENABLE_PHASING}"
echo "[INFO] ENABLE LONG INDEL CALLING: ${ENABLE_LONG_INDEL}"
echo $''

# file check
if [ ! -f ${BAM_FILE_PATH_C} ]; then echo -e "${ERROR} BAM [C] file ${BAM_FILE_PATH_C} not found${NC}"; exit 1; fi
if [ ! -f ${BAM_FILE_PATH_P1} ]; then echo -e "${ERROR} BAM [P1] file ${BAM_FILE_PATH_P1} not found${NC}"; exit 1; fi
if [ ! -f ${BAM_FILE_PATH_P2} ]; then echo -e "${ERROR} BAM [P2] file ${BAM_FILE_PATH_P2} not found${NC}"; exit 1; fi
if [ ! -f ${BAM_FILE_PATH_C}.bai ] && [ ! -f ${BAM_FILE_PATH_C%.*}.bai ]; then echo -e "${ERROR} BAM [C] index bai file not found, please use 'samtools index \$BAM' first${NC}"; exit 1; fi
if [ ! -f ${BAM_FILE_PATH_P1}.bai ] && [ ! -f ${BAM_FILE_PATH_P1%.*}.bai ]; then echo -e "${ERROR} BAM [P1] index bai file not found, please use 'samtools index \$BAM' first${NC}"; exit 1; fi
if [ ! -f ${BAM_FILE_PATH_P2}.bai ] && [ ! -f ${BAM_FILE_PATH_P2%.*}.bai ]; then echo -e "${ERROR} BAM [P2] index bai file not found, please use 'samtools index \$BAM' first${NC}"; exit 1; fi
if [ ! -f ${REFERENCE_FILE_PATH} ]; then echo -e "${ERROR} Reference file ${REFERENCE_FILE_PATH} not found${NC}"; exit 1; fi
if [ ! -f ${REFERENCE_FILE_PATH}.fai ] && [ ! -f ${REFERENCE_FILE_PATH%.*}.fai ]; then echo -e "${ERROR} Reference index fai file not found, please use 'samtools faidx \$REF' first${NC}"; exit 1; fi

if [ "${BED_FILE_PATH}" != "EMPTY" ] && [ ! -z ${BED_FILE_PATH} ] && [ ! -f ${BED_FILE_PATH} ]; then echo -e "${ERROR} BED file ${BED_FILE_PATH} provides but not found${NC}"; exit 1; fi
if [ "${VCF_FILE_PATH}" != "EMPTY" ] && [ ! -z ${VCF_FILE_PATH} ] && [ ! -f ${VCF_FILE_PATH} ]; then echo -e "${ERROR} VCF file ${VCF_FILE_PATH} provides but not found${NC}"; exit 1; fi
if [ ! -d ${MODEL_PATH_C3} ] && [ -z ${CONDA_PREFIX} ]; then echo -e "${ERROR} Conda prefix not found, please activate clair3 conda environment first, model path: ${MODEL_PATH_C3}${NC}"; exit 1; fi
if [ ! -d ${MODEL_PATH_C3} ]; then echo -e "${ERROR} Model path not found${NC}"; exit 1; fi
if [ ! -d ${MODEL_PATH_C3T} ]; then echo -e "${ERROR} Model path not found${NC}"; exit 1; fi

# max threads detection
MAX_THREADS=$(nproc)
if [[ ! ${THREADS} =~ ^[\-0-9]+$ ]] || (( ${THREADS} <= 0)); then echo -e "${ERROR} Invalid threads input --threads=INT ${NC}"; exit 1; fi
if [[ ${THREADS} -gt ${MAX_THREADS} ]]; then echo -e "${WARNING} Threads setting exceeds maximum available threads ${MAX_THREADS}, set threads=${MAX_THREADS}${NC}"; THREADS=${MAX_THREADS}; fi

# max user ulimit threads detection
MAX_ULIMIT_THREADS=`ulimit -u`
if [ ! -z ${MAX_ULIMIT_THREADS} ]; then PER_ULIMIT_THREADS=$((${MAX_ULIMIT_THREADS}/30)); else MAX_ULIMIT_THREADS="unlimited"; PER_ULIMIT_THREADS=${THREADS}; fi
if [[ ${PER_ULIMIT_THREADS} < 1 ]]; then PER_ULIMIT_THREADS=1; fi
if [ "${MAX_ULIMIT_THREADS}" != "unlimited" ] && [[ ${THREADS} -gt ${PER_ULIMIT_THREADS} ]]; then echo -e "${WARNING} Threads setting exceeds maximum ulimit threads ${THREADS} * 30 > ${MAX_ULIMIT_THREADS} (ulimit -u), set threads=${PER_ULIMIT_THREADS}${NC}"; THREADS=${PER_ULIMIT_THREADS}; fi

# platform check
# if [ ! ${PLATFORM} = "ont" ] && [ ! ${PLATFORM} = "hifi" ] && [ ! ${PLATFORM} = "ilmn" ]; then echo -e "${ERROR} Invalid platform input, optional: {ont, hifi, ilmn}${NC}"; exit 1; fi
if [ ! ${PLATFORM} = "ont" ]; then echo -e "${ERROR} Invalid platform input, optional: {ont}${NC}"; exit 1; fi

# optional parameter detection
if [ -z ${BED_FILE_PATH} ]; then echo -e "${ERROR} Use '--bed_fn=FILE' instead of '--bed_fn FILE' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${VCF_FILE_PATH} ]; then echo -e "${ERROR} Use '--vcf_fn=FILE' instead of '--vcf_fn =FILE' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${CONTIGS} ]; then echo -e "${ERROR} Use '--ctg_name=STR' instead of '--ctg_name STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${SAMPLE_C} ]; then echo -e "${ERROR} Use '--sample_name_c=STR' instead of '--sample_name_c STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${SAMPLE_P1} ]; then echo -e "${ERROR} Use '--sample_name_p1=STR' instead of '--sample_name_p1 STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${SAMPLE_P2} ]; then echo -e "${ERROR} Use '--sample_name_p2=STR' instead of '--sample_name_p2 STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${QUAL} ]; then echo -e "${ERROR} Use '--qual=INT' instead of '--qual INT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${SAMTOOLS} ]; then echo -e "${ERROR} Use '--samtools=STR' instead of '--samtools STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PYTHON} ]; then echo -e "${ERROR} Use '--python=STR' instead of '--python STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PYPY} ]; then echo -e "${ERROR} Use '--pypy=STR' instead of '--pypy STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PARALLEL} ]; then echo -e "${ERROR} Use '--parallel=STR' instead of '--parallel STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${WHATSHAP} ]; then echo -e "${ERROR} Use '--whatshap=STR' instead of '--whatshap STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${CHUNK_SIZE} ]; then echo -e "${ERROR} Use '--chunk_size=INT' instead of '--chunk_size INT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${SNP_AF} ]; then echo -e "${ERROR} Use '--snp_min_af=FLOAT' instead of '--snp_min_af FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${INDEL_AF} ]; then echo -e "${ERROR} Use '--indel_min_af=FLOAT' instead of '--indel_min_af FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PRO} ]; then echo -e "${ERROR} Use '--var_pct_full=FLOAT' instead of '--var_pct_full FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${REF_PRO} ]; then echo -e "${ERROR} Use '--ref_pct_full=FLOAT' instead of '--ref_pct_full FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PHASING_PCT} ]; then echo -e "${ERROR} Use '--var_pct_phasing=FLOAT' instead of '--var_pct_phasing FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PILEUP_PREFIX} ]; then echo -e "${ERROR} Use '--pileup_model_prefix=STR' instead of '--pileup_model_prefix STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${TRIO_PREFIX} ]; then echo -e "${ERROR} Use '--trio_model_prefix=STR' instead of '--trio_model_prefix STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${RESUMN} ]; then echo -e "${ERROR} Use '--resumn=0,1,2,3,4'for optional parameters${NC}"; exit 1 ; fi

# model prefix detection
# if [ ! -f ${MODEL_PATH}/${PILEUP_PREFIX}.index ]; then echo -e "${ERROR} No pileup model found in provided model path and model prefix ${MODEL_PATH}/${PILEUP_PREFIX} ${NC}"; exit 1; fi
if [ ! -f ${MODEL_PATH_C3T}/${TRIO_PREFIX}.index ]; then echo -e "${ERROR} No trio model found in provided model path and model prefix ${MODEL_PATH_C3T}/${TRIO_PREFIX} ${NC}"; exit 1; fi


set -x
${SCRIPT_PATH}/trio/Call_Clair3_Trio.sh \
    --bam_fn_c ${BAM_FILE_PATH_C} \
    --bam_fn_p1 ${BAM_FILE_PATH_P1} \
    --bam_fn_p2 ${BAM_FILE_PATH_P2} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --threads ${THREADS} \
    --model_path_clair3 ${MODEL_PATH_C3} \
    --model_path_clair3_trio ${MODEL_PATH_C3T} \
    --platform ${PLATFORM} \
    --output ${OUTPUT_FOLDER} \
    --bed_fn=${BED_FILE_PATH} \
    --vcf_fn=${VCF_FILE_PATH} \
    --ctg_name=${CONTIGS} \
    --sample_name_c=${SAMPLE_C} \
    --sample_name_p1=${SAMPLE_P1} \
    --sample_name_p2=${SAMPLE_P2} \
    --chunk_num=${CHUNK_NUM} \
    --chunk_size=${CHUNK_SIZE} \
    --samtools=${SAMTOOLS} \
    --python=${PYTHON} \
    --pypy=${PYPY} \
    --parallel=${PARALLEL} \
    --whatshap=${WHATSHAP} \
    --qual=${QUAL} \
    --var_pct_full=${PRO} \
    --ref_pct_full=${REF_PRO} \
    --var_pct_phasing=${PHASING_PCT} \
    --snp_min_af=${SNP_AF} \
    --indel_min_af=${INDEL_AF} \
    --pileup_only=${PILEUP_ONLY} \
    --pileup_phasing=${PILEUP_PHASING} \
    --gvcf=${GVCF} \
    --fast_mode=${FAST_MODE} \
    --call_snp_only=${SNP_ONLY} \
    --print_ref_calls=${SHOW_REF} \
    --haploid_precise=${HAP_PRE} \
    --haploid_sensitive=${HAP_SEN} \
    --include_all_ctgs=${INCLUDE_ALL_CTGS} \
    --no_phasing_for_fa=${NO_PHASING} \
    --pileup_model_prefix=${PILEUP_PREFIX} \
    --trio_model_prefix=${TRIO_PREFIX} \
    --remove_intermediate_dir=${RM_TMP_DIR} \
    --enable_phasing=${ENABLE_PHASING} \
    --enable_long_indel=${ENABLE_LONG_INDEL}


)) |& tee ${OUTPUT_FOLDER}/run_clair3_trio.log





