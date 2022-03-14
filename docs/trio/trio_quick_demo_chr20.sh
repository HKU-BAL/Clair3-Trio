#!/bin/bash
_INPUT_DIR=/autofs/bal31/jhsu/home/projects/git/c3t_server_date/demo_ori
_OUTPUT_DIR=/autofs/bal31/jhsu/home/projects/git/c3t_server_date/demo_ori_output
_MODEL_DIR_C3="/autofs/bal31/jhsu/home/projects/git/clair3-trio/models/ont"
_MODEL_DIR_C3T="/autofs/bal31/jhsu/home/projects/git/c3t_server_date/clair3_trio_models/c3t_hg002_g422"

_BAM_C=${_INPUT_DIR}/HG002_60.bam        
_BAM_P1=${_INPUT_DIR}/HG003_60.bam          
_BAM_P2=${_INPUT_DIR}/HG004_60.bam          
_SAMPLE_C="HG002"
_SAMPLE_P1="HG003"
_SAMPLE_P2="HG004"
_REF=${_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna  
_THREADS="36"                
_CONTIGS="chr20"


# run Clair3-Trio
./run_clair3_trio.sh \
  --bam_fn_c=${_BAM_C} \
  --bam_fn_p1=${_BAM_P1} \
  --bam_fn_p2=${_BAM_P2} \
  --output=${_OUTPUT_DIR} \
  --ref_fn=${_REF} \
  --threads=${_THREADS} \
  --model_path_clair3="${_MODEL_DIR_C3}" \
  --model_path_clair3_trio="${_MODEL_DIR_C3T}" \
  --sample_name_c=${_SAMPLE_C} \
  --sample_name_p1=${_SAMPLE_P1} \
  --sample_name_p2=${_SAMPLE_P2}


# benchmarking
BASELINE_VCF_FILE_PATH_C="HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BASELINE_BED_FILE_PATH_C="HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
BASELINE_VCF_FILE_PATH_P1="HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BASELINE_BED_FILE_PATH_P1="HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
BASELINE_VCF_FILE_PATH_P2="HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
BASELINE_BED_FILE_PATH_P2="HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"
OUTPUT_VCF_FILE_C="HG002.vcf.gz"
OUTPUT_VCF_FILE_P1="HG003.vcf.gz"
OUTPUT_VCF_FILE_P2="HG004.vcf.gz"



# Benchmark variant calls against truth set with hap.py
mkdir -p ${_OUTPUT_DIR}/hap

BASELINE_VCF_FILE_PATH=${BASELINE_VCF_FILE_PATH_C}
BASELINE_BED_FILE_PATH=${BASELINE_BED_FILE_PATH_C}
OUTPUT_VCF_FILE_PATH=${OUTPUT_VCF_FILE_C}
HAPPY_PATH=hap/${_SAMPLE_C}_happy

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
${_OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
-f "${_INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
-r "${_REF}" \
-o "${_OUTPUT_DIR}/${HAPPY_PATH}" \
-l ${_CONTIGS} \
--engine=vcfeval \
--threads="8" \
--pass-only

BASELINE_VCF_FILE_PATH=${BASELINE_VCF_FILE_PATH_P1}
BASELINE_BED_FILE_PATH=${BASELINE_BED_FILE_PATH_P1}
OUTPUT_VCF_FILE_PATH=${OUTPUT_VCF_FILE_P1}
HAPPY_PATH=hap/${_SAMPLE_P1}_happy

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
${_OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
-f "${_INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
-r "${_REF}" \
-o "${_OUTPUT_DIR}/${HAPPY_PATH}" \
-l ${_CONTIGS} \
--engine=vcfeval \
--threads="8" \
--pass-only


BASELINE_VCF_FILE_PATH=${BASELINE_VCF_FILE_PATH_P2}
BASELINE_BED_FILE_PATH=${BASELINE_BED_FILE_PATH_P2}
OUTPUT_VCF_FILE_PATH=${OUTPUT_VCF_FILE_P2}
HAPPY_PATH=hap/${_SAMPLE_P2}_happy

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
${_OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
-f "${_INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
-r "${_REF}" \
-o "${_OUTPUT_DIR}/${HAPPY_PATH}" \
-l ${_CONTIGS} \
--engine=vcfeval \
--threads="8" \
--pass-only


# Calculate number of Mendelian violation and de novo variants

mkdir -p ${_OUTPUT_DIR}/trio
M_VCF=${_OUTPUT_DIR}/trio/${_SAMPLE_C}_TRIO.vcf.gz
M_VCF_annotated=${_OUTPUT_DIR}/trio/${_SAMPLE_C}_TRIO_ann.vcf.gz

# merge predicted VCFs
docker run \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
staphb/bcftools:1.12 bcftools merge \
${_OUTPUT_DIR}/${_SAMPLE_C}.vcf.gz \
${_OUTPUT_DIR}/${_SAMPLE_P1}.vcf.gz \
${_OUTPUT_DIR}/${_SAMPLE_P2}.vcf.gz \
--threads 8 -f PASS -0 -m all -O z -o ${M_VCF}

docker run \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
staphb/bcftools:1.12 bcftools index ${M_VCF}

# data preparing
docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
realtimegenomics/rtg-tools:3.12.1 format \
-o ${_INPUT_DIR}/GRCh38_no_alt_analysis_set.sdf "${_REF}"

# generarte BED file
FILE="${_OUTPUT_DIR}/trio.ped"
cat <<EOM >$FILE
#PED format pedigree
#fam-id ind-id pat-id mat-id sex phen
1 HG002 HG003 HG004 1 0
1 HG003 0 0 1 0
1 HG004 0 0 2 0
EOM

# get Mendelian vilations
_TRIO_PED=${_OUTPUT_DIR}/trio.ped
REF_SDF_FILE_PATH=${_INPUT_DIR}/GRCh38_no_alt_analysis_set.sdf 
docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
-v "${_OUTPUT_DIR}":"${_OUTPUT_DIR}" \
realtimegenomics/rtg-tools:3.12.1 mendelian \
-i ${M_VCF} -o ${M_VCF_annotated} --pedigree ${_TRIO_PED} -t ${REF_SDF_FILE_PATH} |& tee ${_OUTPUT_DIR}/trio/MDL.log


# benchmark de novo variants
# note that checking de novo variants require a region specific in ${_CONTIGS}

# merge trio's bed file
BED2=${_INPUT_DIR}/${BASELINE_BED_FILE_PATH_C}
BED3=${_INPUT_DIR}/${BASELINE_BED_FILE_PATH_P1}
BED4=${_INPUT_DIR}/${BASELINE_BED_FILE_PATH_P2}
_TRIO_BED_PATH=${_INPUT_DIR}/020304.bed

docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools intersect -a ${BED2} -b ${BED3} > ${_INPUT_DIR}/tmp_out
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools sort -i ${_INPUT_DIR}/tmp_out > ${_INPUT_DIR}/0203.bed
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools intersect -a ${_INPUT_DIR}/0203.bed -b ${BED4} > ${_INPUT_DIR}/tmp_out
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools sort -i ${_INPUT_DIR}/tmp_out > ${_TRIO_BED_PATH}

# merge trio's true set
_TRIO_GIAB_MERGED=${_INPUT_DIR}/${_SAMPLE_C}_TRIO.vcf.gz

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
staphb/bcftools:1.12 bcftools merge \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH_C} \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH_P1} \
${_INPUT_DIR}/${BASELINE_VCF_FILE_PATH_P2} \
--threads 8 -f PASS -0 -m all -O z -o ${_TRIO_GIAB_MERGED}

docker run \
-v "${_INPUT_DIR}":"${_INPUT_DIR}" \
staphb/bcftools:1.12 bcftools index ${_TRIO_GIAB_MERGED}

# get de nove variants
python clair3.py Check_de_novo --call_vcf ${M_VCF} --ctgName ${_CONTIGS} --bed_fn ${_TRIO_BED_PATH} --true_vcf ${_TRIO_GIAB_MERGED} |& tee ${_OUTPUT_DIR}/trio/denovo_rst


