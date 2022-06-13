# Train a model for Clair3-trio trio calling (revision 1)

This document shows how to train and fine-tune a deep learning model for Clair3-trio trio calling. For training a model for pileup calling, please check [here](../pileup_training.md). The training materials are grouped according to sample, coverage, and chromosome. The groups are converted into tensor binaries. The binaries are much space-efficient and easier to process. As required, multiples tensor binaries can be used together for model training and fine-tuning. 

We provided a step by step scripts to train a Clair3-Trio with the example of data from HG002 trio, please check it if needed.
- [0_generate_trio_bed.sh](0_generate_trio_bed.sh)
- [1_run_RU.sh](1_run_RU.sh)
- [2_genreate_downsample_phased_bam.sh](2_genreate_downsample_phased_bam.sh)
- [3_generate_downsample_pileup.sh](3_generate_downsample_pileup.sh)
- [4_create_tensors.sh](4_create_tensors.sh)
- [5_train.sh](5_train.sh)


---

## Prerequisites

- Clair3 installed
- Clair3-Trio installed
- GNU Parallel installed
- Sufficient hard-disk space
- Truth VCF file after representation unification (check [here](representation_unification_trio.md) on how to generate unified VCF)
- A high-end GPU (have tested in RTX Titan, and RTX 2080Ti)

---

## Contents

  - [1. Setup variables](#1-setup-variables)
  - [2. Run Clair3 pileup model](#2-run-clair3-pileup-model)
  - [3. Create and merge trio tensors](#3-create-and-merge-trio-tensors)
  - [4. Train a Clair3-Trio model](#4-train-a-clair3-trio-model)
  - [5. Finetune a Clair3-Trio model](#5-finetune-a-clair3-trio-model)

---



## 1. Setup variables


trio model utilizes phased alignment, phased alignments is required for training a trio model. 

please using this [page](representation_unification_trio.md) to run phasing alignments in all your family samples.

Note that all samples in family should be phased before training.



```bash
PARALLEL="[WHATSHAP_PATH]"                           # e.g. "parallel"
PYPY="[PYPY_PATH]"                                   # e.g. "pypy3"
SAMTOOLS="[SAMTOOLS_PATH]"                           # e.g. "samtools"
PYTHON3="[PYTHON3_PATH]"                             # e.g. "python3"
THREADS=8                                            # threads number
PLATFORM="ont"                      

# Clair3 folder
_ORI_CLAIR3="[CLAIR3_PATH]"
_MODEL_DIR="[CLAIR3_MODEL_PATH]"
C3_THREADS=8                                         # Clair3 threads number


# Clair3-Trio's path
CLAIR3_TRIO="[CLAIR3-TRIO_PATH]/clair3.py"      

# creating working folder
TRAIN_FOLDER_PREFIX="[YOU_TRAINING_FOLDER]"
BUILD_N="[DATA_BUILDING_NAME]"                       # data building data, e.g. "HG002_all"

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
_LOG_SUF=""                         # log file suffix

# input files and parameters

# sample name
ALL_SAMPLE=(
${CHILD_SAMPLE_N}                   # your child sample name
${P1_SAMPLE_N}                      # your parenet-1 sample name
${P2_SAMPLE_N}                      # your parenet-2 sample name
)

TRIO_N="${CHILD_SAMPLE_N}_TRIO"     # your trio name, e.g. HG002_TRIO

DEPTHS=(                            # data coverage
10
10
10
)

# true variants set from Clair3-Trio Representation Unification
# Each line represents one representation-unified path for each input sample
# note the all path have a folder called **var_ru**
# check the representation_unification_trio.md page for more information
# for practical concerns, the representation_unification_trio.md require only run once on the highest depth for each sample, while the low coverage can be sampled from the highest coverage data, i.e. merged.bam in the representation_unification folder

ALL_RU_FILE_PATH=(
"[child representation unificated folder]"
"[parent-1 representation unificated folder]"
"[parent-2 representation unificated folder]"
)

ALL_PHASED_BAM_FILE_PATH=(
"[child representation unificated folder]/merged.bam"          
"[parent-1 representation unificated folder]/merged.bam"
"[parent-2 representation unificated folder]/merged.bam"
)

ALL_REFERENCE_FILE_PATH=(
"[YOUR_REF_FILE]"
"[YOUR_REF_FILE]"
"[YOUR_REF_FILE]"
)

ALL_ORI_BED_FILE_PATH=(
"[YOUR_BED_FILE_CHILD]"
"[YOUR_BED_FILE_PARENET1]"
"[YOUR_BED_FILE_PARENET2]"
)

# training chrosome name, and prefix
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

```
---

## 2. Run Clair3 pileup model

```
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

```

## 3. Create and merge trio tensors

```

# Get output pileup vcf path
ALL_PILEUP_VCF_FILE_PATH=(
"${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($0))]}_${DEPTHS[$(($0))]}"
"${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($1))]}_${DEPTHS[$(($1))]}"
"${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($2))]}_${DEPTHS[$(($2))]}"
)

# set up input for trio
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
        INPUT_PILEUP_VCF_C+=("${ALL_PILEUP_VCF_FILE_PATH[$(($i))]}")
        INPUT_PILEUP_VCF_P1+=("${ALL_PILEUP_VCF_FILE_PATH[$(($i+1))]}")
        INPUT_PILEUP_VCF_P2+=("${ALL_PILEUP_VCF_FILE_PATH[$(($i+2))]}")

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
--candidate_fn ${INDEL_PATH}/${TRIO_N}_{3}_{1}_{2} \
" ::: ${CHR[@]} ::: ${CHUNK_LIST[@]} ::: ${DEPTH_S[@]} |& tee ${LOG_PATH}/MT${_LOG_SUF}.log

IF_CHECK_MCV=0  # whether filter MCV in training data

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
```


## 4. Train a Clair3-Trio model 

```
# Training trio model
TRAIN_N="[YOUR_CLAIR3-TRIO_MODEL_NAME]"
MODEL_FOLDER_PATH="${TRAIN_FOLDER_PREFIX}/train/{TRAIN_N}"					      
mkdir -p ${MODEL_FOLDER_PATH}
cd ${MODEL_FOLDER_PATH}

# training setting
BATCH_SIZE="[YOUR_BATCH_SIZE]"  #training batch size, e.g. 800
add_indel_length=1
MODEL_ARC=NN
MODEL_ALS="Clair3_Trio_Out3"
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

```

## 5. Finetune a Clair3-Trio model 

```
# finetune with MCVLoss
BATCH_SIZE="[YOUR_BATCH_SIZE]"  #training batch size, e.g. 800
add_indel_length=1
MODEL_ARC=NN
MODEL_ALS="Clair3_Trio_Out3"
IF_ADD_MCV_LOSS=1
MCVLOSS_ALPHA=0.1


# set pretrained models
PRETRAINED_N="[YOUR_PRETRAINED_CALIR3-TRIO_MODEL_NAME]"
PRETRAINED_M_EPCH="[BEGIN_EPOCH_FROM_PRETRAINED_MODEL]"    #10
PRETRAINED_MODEL="${TRAIN_FOLDER_PREFIX}/train/{PRETAINED_N}/.${PRETRAINED_M_EPCH}"					     

# Training trio model
TRAIN_N="[YOUR_NEW_CALIR3-TRIO_MODEL_NAME]"
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
```
