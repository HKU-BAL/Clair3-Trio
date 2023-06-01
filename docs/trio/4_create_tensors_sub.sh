source /autofs/bal36/jhsu/r10/scripts/data_config_trio.sh

PARALLEL=parallel
PYPY=pypy
SAMTOOLS=samtools
PYTHON3=python3
PLATFORM="ont"
THREADS=8                                            # threads number


# Clair3-Trio's path
CLAIR3_TRIO="/autofs/bal31/jhsu/home/projects/git/clair3-trio-clean-v1/clair3.py"      

# creating working folder
TRAIN_FOLDER_PREFIX="/autofs/bal36/jhsu/r10/output/4_build_tensors"
BUILD_N="SUB_1356"                       # data building data, e.g. "HG002_all"

# Temporary working directories
TRAIN_FOLDER="${TRAIN_FOLDER_PREFIX}"
mkdir -p ${TRAIN_FOLDER}

DATASET_FOLDER_PATH="${TRAIN_FOLDER}/build/${BUILD_N}"
TENSOR_CANDIDATE_FOLDER_PATH="${DATASET_FOLDER_PATH}/tensor_can"
BINS_FOLDER_PATH="${DATASET_FOLDER_PATH}/bins"
READ_FOLDER_PATH="${DATASET_FOLDER_PATH}/read_info"
INDEL_PATH="${DATASET_FOLDER_PATH}/alt_info"
SPLIT_BED_PATH="${DATASET_FOLDER_PATH}/split_beds"
PHASE_VCF_PATH="${DATASET_FOLDER_PATH}/phased_vcf"
PHASE_BAM_PATH="${DATASET_FOLDER_PATH}/phased_bam"
LOG_PATH="${DATASET_FOLDER_PATH}/log"
mkdir -p ${DATASET_FOLDER_PATH}
mkdir -p ${TENSOR_CANDIDATE_FOLDER_PATH}
mkdir -p ${BINS_FOLDER_PATH}
mkdir -p ${READ_FOLDER_PATH}
mkdir -p ${INDEL_PATH}
mkdir -p ${SPLIT_BED_PATH}
mkdir -p ${PHASE_VCF_PATH}
mkdir -p ${PHASE_BAM_PATH}
mkdir -p ${LOG_PATH}
cd ${DATASET_FOLDER_PATH}


# log file suffix name
_LOG_SUF=""                         # log file suffix

# input files and parameters
# training chrosome name, and prefix
CHR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22)
CHR_PREFIX="chr"

# GH (tb update to github)
CHILD_SAMPLE_N="HG002"
P1_SAMPLE_N="HG003"
P2_SAMPLE_N="HG004"

# sample name
ALL_SAMPLE=(
${CHILD_SAMPLE_N}                   # your child sample name
${P1_SAMPLE_N}                      # your parenet-1 sample name
${P2_SAMPLE_N}                      # your parenet-2 sample name
${CHILD_SAMPLE_N}
${P1_SAMPLE_N}
${P2_SAMPLE_N}
${CHILD_SAMPLE_N}
${P1_SAMPLE_N}
${P2_SAMPLE_N}
${CHILD_SAMPLE_N}
${P1_SAMPLE_N}
${P2_SAMPLE_N}
${CHILD_SAMPLE_N}
${P1_SAMPLE_N}
${P2_SAMPLE_N}
${CHILD_SAMPLE_N}
${P1_SAMPLE_N}
${P2_SAMPLE_N}
)

TRIO_N="${CHILD_SAMPLE_N}_TRIO_SUB"     # your trio name, e.g. HG002_TRIO
#echo ${ALL_SAMPLE[@]}

DEPTHS=(                            # data coverage
31
31
31
51
51
51
53
53
53
61
61
61
63
63
63
6550
6550
6550
)

#HOME_DIR="/autofs/bal31/jhsu/home"
#REF_FILE_PATH="${HOME_DIR}/data/reference/grch38_no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
#CHILD_BED_FILE_PATH="${HOME_DIR}/data/giab/${CHILD_SAMPLE_N}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
#P1_BED_FILE_PATH="${HOME_DIR}/data/giab/${P1_SAMPLE_N}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
#P2_BED_FILE_PATH="${HOME_DIR}/data/giab/${P2_SAMPLE_N}_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
#
## GH 
## merge trio's bed file using the gerneate_trio_bed.sh
#_TRIO_BED_PATH="${HOME_DIR}/data/giab/020304.bed"

ALL_REFERENCE_FILE_PATH=(
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
"${REF_FILE_PATH}"
)


ALL_BED_FILE_PATH=(
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
)


# true variants set from Clair3-Trio Representation Unification
# Each line represents one representation-unified path for each input sample
# note the all path have a folder called **var_ru**
# check the representation_unification_trio.md page for more information
# for practical concerns, the representation_unification_trio.md require only run once on the highest depth for each sample, while the low coverage can be sampled from the highest coverage data, i.e. merged.bam in the representation_unification folder

ALL_PILEUP_VCF_FILE_PATH=(
"/autofs/bal36/jhsu/r10/output/3_pileup/HG002_30/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG003_10/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG004_10/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG002_50/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG003_10/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG004_10/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG002_50/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG003_30/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG004_30/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG002_65/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG003_10/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG004_10/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG002_65/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG003_30/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG004_30/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG002_65/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG003_50/pileup.vcf.gz"
"/autofs/bal36/jhsu/r10/output/3_pileup/HG004_50/pileup.vcf.gz"
)

ALL_RU_FILE_PATH=(
"/autofs/bal36/jhsu/r10/output/1_ru/HG002/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG003/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG004/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG002/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG003/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG004/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG002/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG003/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG004/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG002/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG003/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG004/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG002/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG003/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG004/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG002/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG003/"
"/autofs/bal36/jhsu/r10/output/1_ru/HG004/"
)

ALL_PHASED_BAM_FILE_PATH=(
"/autofs/bal36/jhsu/r10/output/2_bam/HG002/HG002_30.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG003/HG003_10.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG004/HG004_10.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG002/HG002_50.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG003/HG003_10.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG004/HG004_10.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG002/HG002_50.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG003/HG003_30.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG004/HG004_30.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG002/HG002_65.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG003/HG003_10.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG004/HG004_10.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG002/HG002_65.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG003/HG003_30.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG004/HG004_30.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG002/HG002_65.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG003/HG003_50.bam"
"/autofs/bal36/jhsu/r10/output/2_bam/HG004/HG004_50.bam"
)


# GH
# set up array for create tensors input
INPUT_PILEUP_VCF_C=()
INPUT_PILEUP_VCF_P1=()
INPUT_PILEUP_VCF_P2=()

TRUE_RU_FILE_C=()
TRUE_RU_FILE_P1=()
TRUE_RU_FILE_P2=()

DEPTH_S=()

# GH
# create trio list for candidates input
for i in $(seq 0 $((${#ALL_SAMPLE[@]}-1)))
do
	
    if [ $(($i % 3)) -eq 0 ]; then
        INPUT_PILEUP_VCF_C+=("${ALL_PILEUP_VCF_FILE_PATH[$(($i))]}")
        INPUT_PILEUP_VCF_P1+=("${ALL_PILEUP_VCF_FILE_PATH[$(($i+1))]}")
        INPUT_PILEUP_VCF_P2+=("${ALL_PILEUP_VCF_FILE_PATH[$(($i+2))]}")

        TRUE_RU_FILE_C+=("${ALL_RU_FILE_PATH[$(($i))]}")
        TRUE_RU_FILE_P1+=("${ALL_RU_FILE_PATH[$(($i+1))]}")
        TRUE_RU_FILE_P2+=("${ALL_RU_FILE_PATH[$(($i+2))]}")
        
        DEPTH_S+=("${DEPTHS[$(($i))]}")
    fi
done

echo ${INPUT_PILEUP_VCF_C[@]}
echo ${INPUT_PILEUP_VCF_P1[@]}
echo ${INPUT_PILEUP_VCF_P2[@]}
echo ${TRUE_RU_FILE_C[@]}
echo ${TRUE_RU_FILE_P1[@]}
echo ${TRUE_RU_FILE_P2[@]}
echo ${DEPTH_S[@]}

# created tensors chunk number for each chr
chunk_num=20
CHUNK_LIST=`seq 1 ${chunk_num}`

# create tensors bin chunk number for each chr
bin_chunk_num=10
BIN_CHUNK_LIST=`seq 1 ${bin_chunk_num}`


# Select trio candidates from pileup candidates using the SelectHetSnp_Trio submodule
time ${PARALLEL} --joblog ${LOG_PATH}/S_fiter_hete_snp_pileup${_LOG_SUF}.log -j${THREADS} \
"${PYPY} ${CLAIR3_TRIO} SelectHetSnp_Trio \
--alt_fn_c {2} \
--alt_fn_p1 {3} \
--alt_fn_p2 {4} \
--var_fn_c {5}/var_ru/var_{1} \
--var_fn_p1 {6}/var_ru/var_{1} \
--var_fn_p2 {7}/var_ru/var_{1} \
--split_folder ${SPLIT_BED_PATH} \
--sampleName ${TRIO_N} \
--depth {8} \
--ref_pct_full 0.2 \
--var_pct_full 1.0 \
--ref_var_max_ratio 1.0 \
--for_train 1 \
--chunk_num ${chunk_num} \
--bed ${_TRIO_BED_PATH} \
--ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${INPUT_PILEUP_VCF_C[@]} :::+ ${INPUT_PILEUP_VCF_P1[@]} :::+ ${INPUT_PILEUP_VCF_P2[@]} :::+ ${TRUE_RU_FILE_C[@]} :::+ ${TRUE_RU_FILE_P1[@]} :::+ ${TRUE_RU_FILE_P2[@]} :::+ ${DEPTH_S[@]} |& tee ${LOG_PATH}/FHSP${_LOG_SUF}.log

echo "[INFO] Create Tensors"
time ${PARALLEL} --joblog ${LOG_PATH}/S_create_tensor${_LOG_SUF}.log -j${THREADS} \
"${PYPY} ${CLAIR3_TRIO} CreateTensorFullAlignment \
--bam_fn {4} \
--ref_fn {5} \
--tensor_can_fn ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_{2}_{3}_{1}_{6} \
--indel_fn ${INDEL_PATH}/{2}_{3}_{1}_{6} \
--ctgName ${CHR_PREFIX}{1} \
--samtools ${SAMTOOLS} \
--platform ${PLATFORM} \
--full_aln_regions ${SPLIT_BED_PATH}/${TRIO_N}_{3}_{1}_{6} \
--add_no_phasing_data_training \
--phasing_info_in_bam" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_PHASED_BAM_FILE_PATH[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} ::: ${CHUNK_LIST[@]} |& tee ${LOG_PATH}/CT${_LOG_SUF}.log

# merge trio tensor, noted that merged no depth info
time ${PARALLEL} --joblog ${LOG_PATH}/S_merge_tensors${_LOG_SUF}.log -j${THREADS} \
"${PYTHON3} ${CLAIR3_TRIO} Merge_Tensors_Trio \
--tensor_fn_c ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${ALL_SAMPLE[0]}_{3}_{1}_{2} \
--tensor_fn_p1 ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${ALL_SAMPLE[1]}_{3}_{1}_{2} \
--tensor_fn_p2 ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${ALL_SAMPLE[2]}_{3}_{1}_{2} \
--candidate_fn_c ${INDEL_PATH}/${ALL_SAMPLE[0]}_{3}_{1}_{2} \
--candidate_fn_p1 ${INDEL_PATH}/${ALL_SAMPLE[1]}_{3}_{1}_{2} \
--candidate_fn_p2 ${INDEL_PATH}/${ALL_SAMPLE[2]}_{3}_{1}_{2} \
--tensor_fn ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${TRIO_N}_{3}_{1}_{2} \
--candidate_fn ${INDEL_PATH}/${TRIO_N}_{3}_{1}_{2} \
" ::: ${CHR[@]} ::: ${CHUNK_LIST[@]} ::: ${DEPTH_S[@]} |& tee ${LOG_PATH}/MT${_LOG_SUF}.log


IF_CHECK_MCV=1  # whether filter MCV in training data

time ${PARALLEL} --joblog ${LOG_PATH}/S_tensor2Bin${_LOG_SUF}.log -j${THREADS} \
"${PYTHON3} ${CLAIR3_TRIO} Tensor2Bin_Trio \
--tensor_fn ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${TRIO_N}_{3}_{1} \
--var_fn_c {4}/var_ru/var_{1} \
--var_fn_p1 {5}/var_ru/var_{1} \
--var_fn_p2 {6}/var_ru/var_{1} \
--bin_fn ${BINS_FOLDER_PATH}/${TRIO_N}_{3}_{1}_{2} \
--chunk_id {2} \
--chunk_num ${bin_chunk_num} \
--platform ${PLATFORM} \
--allow_duplicate_chr_pos \
--maximum_non_variant_ratio 1.0 \
--check_mcv ${IF_CHECK_MCV} \
--shuffle" ::: ${CHR[@]} ::: ${BIN_CHUNK_LIST[@]} ::: ${DEPTH_S[@]} :::+ ${TRUE_RU_FILE_C[@]} :::+ ${TRUE_RU_FILE_P1[@]} :::+ ${TRUE_RU_FILE_P2[@]} |& tee ${LOG_PATH}/T2B${_LOG_SUF}.log

