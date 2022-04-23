PARALLEL='parallel'
PYPY="pypy3"
SAMTOOLS="samtools"
PLATFORM="ont"                      
THREADS=8
PYTHON3='python3'

# Clair3 folder
_ORI_CLAIR3=/autofs/bal31/jhsu/home/projects/clair3_t/Clair3_gh/
_MODEL_DIR=/autofs/bal31/jhsu/home/projects/clair3_t/Clair3_gh/models/ont
C3_THREADS=8


# Clair3-Trio folder
CLAIR3_TRIO="${HOME_DIR}/projects/clair3_t/clair3_trio/clair3.py"

# creating working folder
TRAIN_FOLDER_PREFIX="/autofs/bal31/jhsu/home/data/clair3/train/ont_trio_build"
BUILD_N="ALL_1"


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
PILEUP_OUTPUT_PATH="${DATASET_FOLDER_PATH}/pileup"
mkdir -p ${DATASET_FOLDER_PATH}
mkdir -p ${TENSOR_CANDIDATE_FOLDER_PATH}
mkdir -p ${BINS_FOLDER_PATH}
mkdir -p ${READ_FOLDER_PATH}
mkdir -p ${INDEL_PATH}
mkdir -p ${SPLIT_BED_PATH}
mkdir -p ${PHASE_VCF_PATH}
mkdir -p ${PHASE_BAM_PATH}
mkdir -p ${LOG_PATH}
mkdir -p ${PILEUP_OUTPUT_PATH}
cd ${DATASET_FOLDER_PATH}

# log file suffix name
_LOG_SUF=""


# input files

# sample name
ALL_SAMPLE=(
${CHILD_SAMPLE_N}
${P1_SAMPLE_N}
${P2_SAMPLE_N}
)

TRIO_N="${CHILD_SAMPLE_N}_TRIO"

DEPTHS=(
10
10
10
)



# true variants set from Clair3-Trio Representation Unification
# Each line represents one representation-unified path for each input sample
# note the all path have a folder called var_ru
ALL_RU_FILE_PATH=(
"/autofs/bal31/jhsu/home/data/clair3/train/ont_ru_new/HG002/"
"/autofs/bal31/jhsu/home/data/clair3/train/ont_ru_new/HG003/"
"/autofs/bal31/jhsu/home/data/clair3/train/ont_ru_new/HG004/"
)


ALL_PHASED_BAM_FILE_PATH=(
"/autofs/bal31/jhsu/home/data/ont_downsample/HG002/HG002_30.bam"
"/autofs/bal31/jhsu/home/data/ont_downsample/HG003/HG003_30.bam"
"/autofs/bal31/jhsu/home/data/ont_downsample/HG004/HG004_30.bam"
)


ALL_REFERENCE_FILE_PATH=(
${REF_FILE_PATH}
${REF_FILE_PATH}
${REF_FILE_PATH}
)

ALL_ORI_BED_FILE_PATH=(
${HG002_GIAB_BED}
${HG003_GIAB_BED}
${HG004_GIAB_BED}
)

CHR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22)
CHR_PREFIX="chr"

# merge trio's bed file
BED2=${ALL_ORI_BED_FILE_PATH[0]}
BED3=${ALL_ORI_BED_FILE_PATH[1]}
BED4=${ALL_ORI_BED_FILE_PATH[2]}
_TRIO_BED_PATH=${DATASET_FOLDER_PATH}/020304.bed

docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools intersect -a ${BED2} -b ${BED3} > ${_INPUT_DIR}/tmp_out
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools sort -i ${_INPUT_DIR}/tmp_out > ${_INPUT_DIR}/0203.bed
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools intersect -a ${_INPUT_DIR}/0203.bed -b ${BED4} > ${_INPUT_DIR}/tmp_out
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools sort -i ${_INPUT_DIR}/tmp_out > ${_TRIO_BED_PATH}

ALL_BED_FILE_PATH=(
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
${_TRIO_BED_PATH}
)

# Run Clair3 pileup model
time ${PARALLEL} -j ${C3_THREADS} --joblog  ${LOG_PATH}/input_pileup${_LOG_SUF}.log ${_ORI_CLAIR3}/run_clair3.sh \
  --bam_fn={5} \
  --ref_fn={2} \
  --threads=${C3_THREADS} \
  --platform="ont" \
  --model_path="${_MODEL_DIR}" \
  --output=${PILEUP_OUTPUT_PATH}/{1}_{4} \
  --bed_fn={3} \
  --pileup_only ::: ${ALL_SAMPLE[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} :::+ ${ALL_BED_FILE_PATH[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_PHASED_BAM_FILE_PATH[@]}



INPUT_PILEUP_VCF_C=()
INPUT_PILEUP_VCF_P1=()
INPUT_PILEUP_VCF_P2=()

TRUE_RU_FILE_C=()
TRUE_RU_FILE_P1=()
TRUE_RU_FILE_P2=()

DEPTH_S=()

# create trio list for candidates input
for i in $(seq 0 $((${#ALL_SAMPLE[@]}-1)))
do
    if [ $(($i % 3)) -eq 0]; then
        INPUT_PILEUP_VCF_C+=("${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($i))]}_${DEPTHS[$(($i))]}")
        INPUT_PILEUP_VCF_P1+=("${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($i+1))]}_${DEPTHS[$(($i+1))]}")
        INPUT_PILEUP_VCF_P2+=("${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($i+2))]}_${DEPTHS[$(($i+2))]}")

        TRUE_RU_FILE_C+=("${ALL_RU_FILE_PATH[$(($i))]}")
        TRUE_RU_FILE_P1+=("${ALL_RU_FILE_PATH[$(($i+1))]}")
        TRUE_RU_FILE_P2+=("${ALL_RU_FILE_PATH[$(($i+2))]}")
        
        DEPTH_S+=("${DEPTHS[$(($i))]}")
    fi
done


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
--candidate_fn ${INDEL_PATH}/${TRIO_N}_{1}_{2} \
" ::: ${CHR[@]} ::: ${CHUNK_LIST[@]} ::: ${DEPTH_S[@]} |& tee ${LOG_PATH}/MT${_LOG_SUF}.log

IF_CHECK_MCV=1

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


# Merge compressed binaries
ALL_BINS_FOLDER_PATH="${BINS_FOLDER_PATH}/../all_bins"
mkdir -p ${ALL_BINS_FOLDER_PATH}
${PARALLEL} --joblog ${LOG_PATH}/S_mergeBin${_LOG_SUF}.log -j${THREADS} \
echo "${PYTHON3} ${CLAIR3_TRIO} MergeBin_Trio \
    ${BINS_FOLDER_PATH}/${TRIO_N}_{2}_{1}_* \
    --out_fn ${ALL_BINS_FOLDER_PATH}/bin_{2}_{1}" ::: ${CHR[@]} ::: ${DEPTH_S[@]}


# Training trio model
TRAIN_N="test_train"
MODEL_FOLDER_PATH="${TRAIN_FOLDER_PREFIX}/train/{TRAIN_N}"					      
mkdir -p ${MODEL_FOLDER_PATH}
cd ${MODEL_FOLDER_PATH}

# training setting
add_indel_length=1
MODEL_ARC=NN
MODEL_ALS="Clair3_Trio_Out3"
BATCH_SIZE=800
IF_ADD_MCV_LOSS=0
MCVLOSS_ALPHA=0

# A single GPU is used for model training
export CUDA_VISIBLE_DEVICES="0"

echo "[INFO] Model training"
time ${PYTHON3} ${CLAIR3_TRIO} Train_Trio \
--bin_fn ${ALL_BINS_FOLDER_PATH} \
--ochk_prefix ${MODEL_FOLDER_PATH} \
--add_indel_length ${add_indel_length} \
--platform ${PLATFORM} \
--validation_dataset \
--learning_rate 0.001 \
--maxEpoch 30 \
--model_arc ${MODEL_ARC} \
--model_cls ${MODEL_ALS} \
--batch_size ${BATCH_SIZE} \
--add_mcv_loss ${IF_ADD_MCV_LOSS} \
--mcv_alpha ${MCVLOSS_ALPHA} \
 |& tee ${MODEL_FOLDER_PATH}/train_log


# finetune with MCVLoss
add_indel_length=1
MODEL_ARC=NN
MODEL_ALS="Clair3_Trio_Out3"
BATCH_SIZE=800
IF_ADD_MCV_LOSS=1
MCVLOSS_ALPHA=0.1


# set pretrained models
PRETRAINED_N="test_train"
PRETRAINED_M_EPCH='10'    #01, 10, 30
PRETRAINED_MODEL="${TRAIN_FOLDER_PREFIX}/train/{PRETAINED_N}/.${PRETRAINED_M_EPCH}"					     

# Training trio model
TRAIN_N="test_train_finetune"
MODEL_FOLDER_PATH="${TRAIN_FOLDER_PREFIX}/train/{TRAIN_N}"					      
mkdir -p ${MODEL_FOLDER_PATH}
cd ${MODEL_FOLDER_PATH}

echo "[INFO] Model training"
time ${PYTHON3} ${CLAIR3_TRIO} Train_Trio \
--bin_fn ${ALL_BINS_FOLDER_PATH} \
--ochk_prefix ${MODEL_FOLDER_PATH} \
--add_indel_length ${add_indel_length} \
--platform ${PLATFORM} \
--validation_dataset \
--learning_rate 1e-5 \
--maxEpoch 10 \
--model_arc ${MODEL_ARC} \
--batch_size ${BATCH_SIZE} \
--chkpnt_fn ${PRETRAINED_MODEL} \
--add_mcv_loss ${IF_ADD_MCV_LOSS} \
--mcv_alpha ${MCVLOSS_ALPHA} \
 |& tee ${MODEL_FOLDER_PATH}/train_log
