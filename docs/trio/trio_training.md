# Train a model for Clair3-trio trio calling (v3)

2023/06/01

This document shows how to train and fine-tune a deep learning model for Clair3-trio trio calling. For training a model for pileup calling, please check [here](../pileup_training.md). The training materials are grouped according to sample, coverage, and chromosome. The groups are converted into tensor binaries. The binaries are much space-efficient and easier to process. As required, multiples tensor binaries can be used together for model training and fine-tuning. 

We provided a list of step-by-step scripts to train a Clair3-Trio with the example of data from the HG002 trio, please check it if needed.

- [0_generate_trio_bed.sh](0_generate_trio_bed.sh)
- [1 representation_unification_trio.md](representation_unification_trio.md)
- [2_genreate_downsample_phased_bam.sh](2_genreate_downsample_phased_bam.sh)
- [3_generate_downsample_pileup.sh](3_generate_downsample_pileup.sh)
- [4.1_create even coverage_tensors.sh](4_create_tensors.sh)
- [4.2_create_uneven coverage tensors.sh](4_create_tensors_sub.sh)
- [4.3_gather even and unenven coverage data](4_1_downsample_bin_1.sh)
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
  - [3-1 Create and merge trio tensors even coverage](#3-1-create-and-merge-trio-tensors-even-coverage)
  - [3-2 Create and merge trio tensors uneven coverage](#3-2-create-and-merge-trio-tensors-uneven-coverage)
  - [3-3 Gather all bins](#3-3-gather-all-bins)
  - [4. Train a Clair3-Trio model](#4-train-a-clair3-trio-model)
  - [5. Finetune a Clair3-Trio model](#5-finetune-a-clair3-trio-model)
---

## 0. Preparation of needed files

The original input files for Clair3-Trio includes:

- reference file: [R]
- child/parents bam files: [C BAM], [P1 BAM], [P2 BAM]
- child/parents truth vcf files: [C TRUTH VCF], [P1 TRUTH VCF], [P2 TRUTH VCF]
- child/parents bed files for truth variants: [C TRUTH BED], [P1 TRUTH BED], [P2 TRUTH BED]

After gathering all above files, a few files need to be generated for Clair3-Trio training in the first place.
- phased alignment for **each individual** by running the [representation_unifivation_trio](https://github.com/HKU-BAL/Clair3-Trio/blob/trio/docs/trio/representation_unification_trio.md) 


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
${CHILD_SAMPLE_N}                   # your child sample name
${P1_SAMPLE_N}                      # your parenet-1 sample name
${P2_SAMPLE_N}                      # your parenet-2 sample name
)

TRIO_N="${CHILD_SAMPLE_N}_TRIO"     # your trio name, e.g. HG002_TRIO

# Note: here we set it to 10x and 30x for example, for using a maximum of 120x data, here, all following paths of (ALL_RU_FILE_PATH, ALL_PHASED_BAM_FILE_PATH, ALL_REFERENCE_FILE_PATH, etc.),  need to be set to 10x, 30x, 50x, 70x, 90x, 120x accordingly.DEPTHS=(                            # data coverage
10
10
10
30
30
30
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
"[child representation unificated folder]"
"[parent-1 representation unificated folder]"
"[parent-2 representation unificated folder]"
)

# downsampled RU bam
ALL_PHASED_BAM_FILE_PATH=(
"[child representation unificated folder]/${CHILD_SAMPLE_N}_10.bam"          
"[parent-1 representation unificated folder]/${P1_SAMPLE_N}_10.bam"
"[parent-2 representation unificated folder]/${P2_SAMPLE_N}_10.bam"
"[child representation unificated folder]/${CHILD_SAMPLE_N}_30.bam"          
"[parent-1 representation unificated folder]/${P1_SAMPLE_N}_30.bam"
"[parent-2 representation unificated folder]/${P2_SAMPLE_N}_30.bam"
)

ALL_REFERENCE_FILE_PATH=(
"[YOUR_REF_FILE]"
"[YOUR_REF_FILE]"
"[YOUR_REF_FILE]"
"[YOUR_REF_FILE]"
"[YOUR_REF_FILE]"
"[YOUR_REF_FILE]"
)

ALL_ORI_BED_FILE_PATH=(
"[YOUR_BED_FILE_CHILD]"
"[YOUR_BED_FILE_PARENET1]"
"[YOUR_BED_FILE_PARENET2]"
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

## 3-1 Create and merge trio tensors even coverage

generate the even coverage for Clair-Trio input, for example (child 10x, parent 1 10x, parent2 10x) + (child 30x, parent 1 30x, parent2 30)

```

# Get output pileup vcf path
ALL_PILEUP_VCF_FILE_PATH=(
"${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($0))]}_${DEPTHS[$(($0))]}/pileup.vcf.gz"
"${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($1))]}_${DEPTHS[$(($1))]}/pileup.vcf.gz"
"${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($2))]}_${DEPTHS[$(($2))]}/pileup.vcf.gz"
"${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($3))]}_${DEPTHS[$(($3))]}/pileup.vcf.gz"
"${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($4))]}_${DEPTHS[$(($4))]}/pileup.vcf.gz"
"${PILEUP_OUTPUT_PATH}/${ALL_SAMPLE[$(($5))]}_${DEPTHS[$(($5))]}/pileup.vcf.gz"
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

```
## 3-2 Create and merge trio tensors uneven coverage

generate the uneven coverage for Clair-Trio input, for example (child 30x, parent 1 10x, parent2 10x) + (child 50x, parent 1 10x, parent2 10)
rerun 1, 2, and 3-1 to generate uneven coverage.
check [4.2_create_uneven coverage tensors.sh](4_create_tensors_sub.sh) for example

for a trio data from [r10](https://labs.epi2me.io/askenazi-kit14-2022-12/) with the following data coverage,

| sample | coverage |
| ------ | -------- |
| HG002  | 70.13    |
| HG003  | 79.66    |
| HG004  | 64.72    |


we recommend generating the below trio data coverage config:

| even coverage   | HG002 | HG003 | HG004 |
| --------------- | ----- | ----- | ----- |
|                 | 10    | 10    | 10    |
|                 | 30    | 30    | 30    |
|                 | 50    | 50    | 50    |
|                 | 65    | 65    | 65    |
|                 |       |       |       |
| uneven coverage | HG002 | HG003 | HG004 |
|                 | 30    | 10    | 10    |
|                 | 50    | 10    | 10    |
|                 | 50    | 30    | 30    |
|                 | 65    | 10    | 10    |
|                 | 65    | 30    | 30    |
|                 | 65    | 50    | 50    |


## 3-3 Gather all bins

gather all bin files for model training.
by copying all bins files from even/uneven coverage `${BINS_FOLDER_PATH}` into `$ALL_BINS_FOLDER_PATH`.

We recommend to downsample bins files, as an example in [4.3_gather even and uneven coverage data](4_1_downsample_bin_1.sh).


## 4. Train a Clair3-Trio model 

please set the `${ALL_BINS_FOLDER_PATH}` which contains all bins files.

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
--ochk_prefix ${MODEL_FOLDER_PATH}/trio \
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
--ochk_prefix ${MODEL_FOLDER_PATH}/trio \
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
