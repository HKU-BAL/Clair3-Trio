# using Clair3 to generating pileup calling on each downsampled sample

PARALLEL=parallel
PYPY=pypy
SAMTOOLS=samtools
PYTHON3=python3
PLATFORM="ont"

# Clair3 folder
_ORI_CLAIR3="XXX"

# note the use right models for your training
# check https://github.com/HKU-BAL/Clair3#pre-trained-models
_MODEL_DIR="XXXX/r941_prom_sup_g5014/"
C3_THREADS=36                                         # Clair3 threads number

# Clair3-Trio's clair3.py path
CLAIR3_TRIO="XXX/clair3.py"      

# output folder for pileup files
DATASET_FOLDER_PATH="XXX"

# creating working folder
PILEUP_OUTPUT_PATH="${DATASET_FOLDER_PATH}/3_pileup"
LOG_PATH="${PILEUP_OUTPUT_PATH}"
mkdir -p ${PILEUP_OUTPUT_PATH}
cd ${PILEUP_OUTPUT_PATH}


# input files and parameters
# training chrosome name, and prefix
# excluded the chr20 for testing
CHR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22) # please update this part for your dataset
CHR_PREFIX="chr"

CHILD_SAMPLE_N="HG002" #please update this part for your dataset
P1_SAMPLE_N="HG003"
P2_SAMPLE_N="HG004"

# bam file from 2_generate_downsample_phased_bam
# please update this part for your dataset
# sample name
# if you have four downsampled dataset, you need to secific 4 * 3 sample name in this setting.
ALL_SAMPLE=(
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


# list number same as in $ALL_SAMPLE
DEPTHS=(                            # data coverage
10
10
10
30
30
30
60
60
60
80
80
80
)

# list number same as in $ALL_SAMPLE
# phased bam files, from the 2_generate_downsample_phased_bam.sh
ALL_PHASED_BAM_FILE_PATH=(
"XXX/bam/phased/HG002/HG002_10.bam"
"XXX/bam/phased/HG003/HG003_10.bam"
"XXX/bam/phased/HG004/HG004_10.bam"
"XXX/bam/phased/HG002/HG002_30.bam"
"XXX/bam/phased/HG003/HG003_30.bam"
"XXX/bam/phased/HG004/HG004_30.bam"
"XXX/bam/phased/HG002/HG002_60.bam"
"XXX/bam/phased/HG003/HG003_60.bam"
"XXX/bam/phased/HG004/HG004_60.bam"
"XXX/bam/phased/HG002/HG002_80.bam"
"XXX/bam/phased/HG003/HG003_80.bam"
"XXX/bam/phased/HG004/HG004_80.bam"
)

# please update this part for your dataset
# HOME_DIR="/autofs/bal31/jhsu/home"
# REF_FILE_PATH="${HOME_DIR}/data/reference/grch38_no_alt_plus_hs38d1/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
# merged trio's bed file using the 0_gerneate_trio_bed.sh
# _TRIO_BED_PATH="${HOME_DIR}/data/giab/020304.bed" 

# your reference file
REF_FILE_PATH="XXXX/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna"
# your trio bed file, from 0_generate_trio_bed.sh
_TRIO_BED_PATH="XXXX/data/giab/020304.bed" 

# list number same as in $ALL_SAMPLE
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
)


# list number same as in $ALL_SAMPLE
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
)


# log file suffix name
_LOG_SUF=""                         # log file suffix

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

