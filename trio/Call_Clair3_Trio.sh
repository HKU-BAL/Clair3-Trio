#!/bin/bash
SCRIPT_NAME=$(basename "$0")
Usage="Usage: ./${SCRIPT_NAME} --bam_fn_c=BAM --bam_fn_p1=BAM --bam_fn_p2=BAM --ref_fn=REF --output=OUTPUT_DIR --threads=THREADS --model_path_clair3=MODEL_PREFIX --model_path_clair3_trio=MODEL_PREFIX [--bed_fn=BED] [options]"
# INFO: whole calling workflow of clair3

set -e
ARGS=`getopt -o f:t:p:o:r::c::s::h::g \
-l bam_fn_c:,bam_fn_p1:,bam_fn_p2:,ref_fn:,threads:,model_path_clair3:,model_path_clair3_trio:,platform:,output:,\
bed_fn::,vcf_fn::,ctg_name::,sample_name_c::,sample_name_p1::,sample_name_p2::,help::,qual::,samtools::,python::,pypy::,parallel::,whatshap::,chunk_num::,chunk_size::,var_pct_full::,var_pct_phasing::,\
snp_min_af::,indel_min_af::,ref_pct_full::,pileup_only::,pileup_phasing::,fast_mode::,gvcf::,print_ref_calls::,haploid_precise::,haploid_sensitive::,include_all_ctgs::,\
no_phasing_for_fa::,pileup_model_prefix::,trio_model_prefix::,call_snp_only::,remove_intermediate_dir::,enable_phasing::,enable_output_haplotagging::,enable_long_indel:: -n 'run_clair3_trio.sh' -- "$@"`

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
PILEUP_PREFIX="pileup"
TRIO_PREFIX="trio"
RM_TMP_DIR=False
ENABLE_PHASING=False
ENABLE_LONG_INDEL=False
ENABLE_OUTPUT_HAPLOTAGGING=False


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
    --pileup_only ) PILEUP_ONLY="$2"; shift 2 ;;
    --pileup_phasing ) PILEUP_PHASING="$2"; shift 2 ;;
    --fast_mode ) FAST_MODE="$2"; shift 2 ;;
    --call_snp_only ) SNP_ONLY="$2"; shift 2 ;;
    --print_ref_calls ) SHOW_REF="$2"; shift 2 ;;
    --gvcf ) GVCF="$2"; shift 2 ;;
    --snp_min_af ) SNP_AF="$2"; shift 2 ;;
    --indel_min_af ) INDEL_AF="$2"; shift 2 ;;
    --pileup_model_prefix ) PILEUP_PREFIX="$2"; shift 2 ;;
    --trio_model_prefix ) TRIO_PREFIX="$2"; shift 2 ;;
    --haploid_precise ) HAP_PRE="$2"; shift 2 ;;
    --haploid_sensitive ) HAP_SEN="$2"; shift 2 ;;
    --include_all_ctgs ) INCLUDE_ALL_CTGS="$2"; shift 2 ;;
    --no_phasing_for_fa ) NO_PHASING="$2"; shift 2 ;;
    --remove_intermediate_dir ) RM_TMP_DIR="$2"; shift 2 ;;
    --enable_long_indel ) ENABLE_LONG_INDEL="$2"; shift 2 ;;
    --enable_phasing ) ENABLE_PHASING="$2"; shift 2 ;;
    --enable_output_haplotagging ) ENABLE_OUTPUT_HAPLOTAGGING="$2"; shift 2 ;;

    -- ) shift; break; ;;
    -h|--help ) print_help_messages; break ;;
    * ) print_help_messages; exit 0 ;;
   esac
done

if [ ${GVCF} == True ]
then
SHOW_REF=True
fi

SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
CLAIR3_TRIO="${SHELL_FOLDER}/../clair3.py"


# if [ ${BED_FILE_PATH} = "EMPTY" ] ; then BED_FILE_PATH= ; fi
RETRIES=4

PILEUP_CHECKPOINT_PATH="${MODEL_PATH}/${PILEUP_PREFIX}"
FULL_ALIGNMENT_CHECKPOINT_PATH="${MODEL_PATH}/${FA_PREFIX}"
LOG_PATH="${OUTPUT_FOLDER}/log"
TMP_FILE_PATH="${OUTPUT_FOLDER}/tmp"
GVCF_TMP_PATH="${TMP_FILE_PATH}/gvcf_tmp_output"

PILEUP_VCF_PATH="${TMP_FILE_PATH}/pileup_output"
# TENSOR_CANDIDATE_FOLDER_PATH="${TMP_FILE_PATH}/tensor_can"
# SPLIT_BED_PATH="${TMP_FILE_PATH}/split_beds"
# INDEL_PATH="${TMP_FILE_PATH}/alt_info"

CALL_PATH="${TMP_FILE_PATH}/trio_output"
CANDIDATE_BED_PATH="${CALL_PATH}/candidate_bed"
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1

echo "[TRIO INFO] * Clir3-Trio pipeline start"
echo "[TRIO INFO] * 0 Check environment variables"

${PYTHON} ${CLAIR3_TRIO} CheckEnvs_Trio \
   --bam_fn_c ${BAM_FILE_PATH_C} \
   --bam_fn_p1 ${BAM_FILE_PATH_P1} \
   --bam_fn_p2 ${BAM_FILE_PATH_P2} \
   --output_fn_prefix ${OUTPUT_FOLDER} \
   --ctg_name ${CONTIGS} \
   --bed_fn ${BED_FILE_PATH} \
   --ref_fn ${REFERENCE_FILE_PATH} \
   --vcf_fn ${VCF_FILE_PATH} \
   --chunk_num ${CHUNK_NUM} \
   --chunk_size ${CHUNK_SIZE} \
   --include_all_ctgs ${INCLUDE_ALL_CTGS} \
   --threads ${THREADS} \
   --python ${PYTHON} \
   --pypy ${PYPY} \
   --samtools ${SAMTOOLS} \
   --whatshap ${WHATSHAP} \
   --parallel ${PARALLEL} \
   --qual ${QUAL} \
   --sampleName_c ${SAMPLE_C} \
   --sampleName_p1 ${SAMPLE_P1} \
   --sampleName_p2 ${SAMPLE_P2} \
   --var_pct_full ${PRO} \
   --ref_pct_full ${REF_PRO} \
   --snp_min_af ${SNP_AF} \
   --indel_min_af ${INDEL_AF}


readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
if [ ${#CHR[@]} -eq 0 ]; then echo "[INFO] Exit in environment checking"; exit 0; fi
THREADS_LOW=$((${THREADS}*3/4))
if [[ ${THREADS_LOW} < 1 ]]; then THREADS_LOW=1; fi

cd ${OUTPUT_FOLDER}

export CUDA_VISIBLE_DEVICES=""
echo "[TRIO INFO] * 1 call variants using pileup model"


ALL_SAMPLE=(
${SAMPLE_C}
${SAMPLE_P1}
${SAMPLE_P2}
)

ALL_UNPHASED_BAM_FILE_PATH=(
${BAM_FILE_PATH_C}
${BAM_FILE_PATH_P1}
${BAM_FILE_PATH_P2}
)

CLAIR3_THREADS=$((${THREADS}/3))
if [[ ${CLAIR3_THREADS} < 1 ]]; then CLAIR3_THREADS=1; fi

echo "trio sample" ${ALL_SAMPLE[@]}
echo "trio bam" ${ALL_UNPHASED_BAM_FILE_PATH[@]}
echo "pileup threads" ${CLAIR3_THREADS}

echo "running pileup model"


${PARALLEL} --retries ${RETRIES} -j3 --joblog  ${LOG_PATH}/parallel_1_clair3_pileup.log \
"${SHELL_FOLDER}/..//scripts/clair3.sh \
   --bam_fn={2} \
   --ref_fn=${REFERENCE_FILE_PATH} \
   --threads=${CLAIR3_THREADS} \
   --model_path=${MODEL_PATH_C3} \
   --output=${PILEUP_VCF_PATH}/{1} \
   --platform="ont" \
   --bed_fn=${BED_FILE_PATH} \
   --vcf_fn=${VCF_FILE_PATH} \
   --ctg_name=${CONTIGS} \
   --sample_name={1} \
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
   --snp_min_af=${SNP_AF} \
   --indel_min_af=${INDEL_AF} \
   --var_pct_phasing=${PHASING_PCT} \
   --pileup_only=False \
   --pileup_phasing=True \
   --gvcf=${GVCF} \
   --fast_mode=${FAST_MODE} \
   --call_snp_only=${SNP_ONLY} \
   --print_ref_calls=${SHOW_REF} \
   --haploid_precise=${HAP_PRE} \
   --haploid_sensitive=${HAP_SEN} \
   --include_all_ctgs=${INCLUDE_ALL_CTGS} \
   --no_phasing_for_fa=${NO_PHASING} \
   --remove_intermediate_dir=${RM_TMP_DIR} \
   --enable_phasing=${ENABLE_PHASING} \
   --enable_long_indel=${ENABLE_LONG_INDEL} \
   --pileup_model_prefix=${PILEUP_PREFIX} \
   --fa_model_prefix=full_alignment" ::: ${ALL_SAMPLE[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]} |& tee ${LOG_PATH}/1_call_var_bam_pileup.log


echo "[TRIO INFO] Finish Pileup model"


# phased bam organized in different chr
ALL_PHASED_BAM_FILE_DIR=(
${PILEUP_VCF_PATH}/${SAMPLE_C}/tmp/phase_output/phase_bam/
${PILEUP_VCF_PATH}/${SAMPLE_P1}/tmp/phase_output/phase_bam/
${PILEUP_VCF_PATH}/${SAMPLE_P2}/tmp/phase_output/phase_bam/
)



INPUT_PILEUP_VCF=(
${PILEUP_VCF_PATH}/${SAMPLE_C}/pileup.vcf.gz
${PILEUP_VCF_PATH}/${SAMPLE_P1}/pileup.vcf.gz
${PILEUP_VCF_PATH}/${SAMPLE_P2}/pileup.vcf.gz
)


TRIO_N="${ALL_SAMPLE}_TRIO"

# note that the phased bam stored in separate files between chromosome
#echo ${CHR[@]}

# select candidate from trio input
echo "[TRIO INFO] 2 Select Trio Candidates"
time ${PARALLEL} --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_2_fiter_hete_snp_pileup.log -j${THREADS} \
"${PYPY} ${CLAIR3_TRIO}  SelectCandidates_Trio \
--alt_fn_c ${INPUT_PILEUP_VCF[0]} \
--alt_fn_p1 ${INPUT_PILEUP_VCF[1]} \
--alt_fn_p2 ${INPUT_PILEUP_VCF[2]} \
--candidate_bed ${CANDIDATE_BED_PATH} \
--sampleName ${TRIO_N} \
--ref_pct_full 0.03 \
--var_pct_full 1.0 \
--ref_var_max_ratio 1.2 \
--ctgName {1}" ::: ${CHR[@]} |& tee ${LOG_PATH}/2_FHSP.log

# using gpu are preserved for debugging only
use_gpu=0
# export CUDA_VISIBLE_DEVICES="0"

echo "[TRIO INFO] * 3 Call Clair3-Trio model"
cat ${CANDIDATE_BED_PATH}/FULL_ALN_FILE_* > ${CANDIDATE_BED_PATH}/FULL_ALN_FILES

time ${PARALLEL} --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_3_callvarbam.log -j${THREADS_LOW} \
"${PYTHON} ${CLAIR3_TRIO} CallVarBam_Trio \
--chkpnt_fn ${MODEL_PATH_C3T}/${TRIO_PREFIX} \
--bam_fn_c ${ALL_PHASED_BAM_FILE_DIR[0]}/{1/.}.bam \
--bam_fn_p1 ${ALL_PHASED_BAM_FILE_DIR[1]}/{1/.}.bam \
--bam_fn_p2 ${ALL_PHASED_BAM_FILE_DIR[2]}/{1/.}.bam \
--sampleName_c ${SAMPLE_C} \
--sampleName_p1 ${SAMPLE_P1} \
--sampleName_p2 ${SAMPLE_P2} \
--call_fn_c ${CALL_PATH}/${SAMPLE_C}/trio_${SAMPLE_C}_{1/}.vcf \
--call_fn_p1 ${CALL_PATH}/${SAMPLE_P1}/trio_${SAMPLE_P1}_{1/}.vcf \
--call_fn_p2 ${CALL_PATH}/${SAMPLE_P2}/trio_${SAMPLE_P2}_{1/}.vcf \
--use_gpu ${use_gpu} \
--ref_fn ${REFERENCE_FILE_PATH} \
--ctgName {1/.} \
--samtools ${SAMTOOLS} \
--platform ${PLATFORM} \
--full_aln_regions {1} \
--gvcf=${GVCF} \
--showRef ${SHOW_REF} \
--phasing_info_in_bam" :::: ${CANDIDATE_BED_PATH}/FULL_ALN_FILES |& tee ${LOG_PATH}/3_CV.log

${PARALLEL}  -j${THREADS} \
"${PYPY} ${CLAIR3_TRIO} SortVcf_Trio \
   --input_dir ${CALL_PATH}/{1} \
   --vcf_fn_prefix "trio_{1}" \
   --output_fn ${OUTPUT_FOLDER}/{1}_c3t.vcf \
   --sampleName {1} \
   --ref_fn ${REFERENCE_FILE_PATH} \
   --contigs_fn ${TMP_FILE_PATH}/CONTIGS" ::: ${ALL_SAMPLE[@]}

if [ "$( gzip -fdc ${OUTPUT_FOLDER}/${SAMPLE_C}_c3t.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then echo "[INFO] Exit in trio variant calling"; exit 0; fi


INPUT_PILEUP_GVCF_PATH=(
${PILEUP_VCF_PATH}/${SAMPLE_C}/tmp/gvcf_tmp_output/
${PILEUP_VCF_PATH}/${SAMPLE_P1}/tmp/gvcf_tmp_output/
${PILEUP_VCF_PATH}/${SAMPLE_P2}/tmp/gvcf_tmp_output/
)

TRIO_M_OUTPUT_VCF=(
${OUTPUT_FOLDER}/${SAMPLE_C}_c3t.vcf.gz
${OUTPUT_FOLDER}/${SAMPLE_P1}_c3t.vcf.gz
${OUTPUT_FOLDER}/${SAMPLE_P2}_c3t.vcf.gz
)

mkdir -p ${GVCF_TMP_PATH}/${SAMPLE_C}
mkdir -p ${GVCF_TMP_PATH}/${SAMPLE_P1}
mkdir -p ${GVCF_TMP_PATH}/${SAMPLE_P2}


echo $''
echo "[TRIO INFO] * 4 Merge trio VCF"
time ${PARALLEL} --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_4_merge_vcf.log -j${THREADS} \
"${PYPY} ${CLAIR3_TRIO} MergeVcf_Trio \
   --pileup_vcf_fn {3} \
   --bed_fn_prefix ${CANDIDATE_BED_PATH} \
   --trio_vcf_fn {4} \
   --output_fn ${GVCF_TMP_PATH}/{2}/merge_{1}.vcf \
   --gvcf_fn ${GVCF_TMP_PATH}/{2}/merge_{1}.gvcf \
   --platform ${PLATFORM} \
   --print_ref_calls ${SHOW_REF} \
   --gvcf ${GVCF} \
   --haploid_precise ${HAP_PRE} \
   --haploid_sensitive ${HAP_SEN} \
   --non_var_gvcf_fn {5}/non_var.gvcf \
   --ref_fn ${REFERENCE_FILE_PATH} \
   --ctgName {1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${INPUT_PILEUP_VCF[@]} :::+ ${TRIO_M_OUTPUT_VCF[@]} :::+ ${INPUT_PILEUP_GVCF_PATH[@]} |& tee ${LOG_PATH}/4_MV.log

${PARALLEL}  -j${THREADS} \
"${PYPY} ${CLAIR3_TRIO} SortVcf_Trio \
   --input_dir ${GVCF_TMP_PATH}/{1} \
   --vcf_fn_prefix "merge" \
   --output_fn ${OUTPUT_FOLDER}/{1}.vcf \
   --sampleName {1} \
   --ref_fn ${REFERENCE_FILE_PATH} \
   --contigs_fn ${TMP_FILE_PATH}/CONTIGS" ::: ${ALL_SAMPLE[@]}

if [ "$( gzip -fdc ${OUTPUT_FOLDER}/${SAMPLE_C}.vcf.gz | grep -v '#' | wc -l )" -eq 0 ]; then echo "[INFO] Exit in trio variant calling"; exit 0; fi

if [ ${GVCF} == True ]
then
${PARALLEL}  -j${THREADS} \
${PYPY} ${CLAIR3_TRIO} SortVcf_Trio \
       --input_dir ${GVCF_TMP_PATH}/{1} \
       --vcf_fn_prefix "merge" \
       --vcf_fn_suffix ".gvcf" \
       --output_fn ${OUTPUT_FOLDER}/{1}.gvcf \
       --sampleName {1} \
       --ref_fn ${REFERENCE_FILE_PATH} \
       --contigs_fn ${TMP_FILE_PATH}/CONTIGS ::: ${ALL_SAMPLE[@]}
fi



if [ ${ENABLE_PHASING} == True ]
then
    echo "[INFO] Phasing VCF output in parallel using WhatsHap"
    time ${PARALLEL} --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_5_phase_vcf_output.log -j${THREADS} \
    "${WHATSHAP} phase \
        --output ${GVCF_TMP_PATH}/{2}/phased_merge_{1}.vcf \
        --reference ${REFERENCE_FILE_PATH} \
        --ignore-read-groups \
        ${GVCF_TMP_PATH}/{2}/merge_{1}.vcf \
        {3}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]} |& tee ${LOG_PATH}/5_phase_vcf_output.log


    ${PARALLEL}  -j${THREADS} \
    "${PYPY} ${CLAIR3_TRIO} SortVcf_Trio \
        --input_dir ${GVCF_TMP_PATH}/{1} \
        --vcf_fn_prefix "phased_merge" \
        --output_fn ${OUTPUT_FOLDER}/phased_{1}.vcf \
        --sampleName {1} \
        --ref_fn ${REFERENCE_FILE_PATH} \
        --contigs_fn ${TMP_FILE_PATH}/CONTIGS" ::: ${ALL_SAMPLE[@]}
    echo "[INFO] Finish Phasing calling, output phased VCF file: ${OUTPUT_FOLDER}/phased_*.vcf.gz"


    # if [ ${GVCF} == True ]
    # then
    # echo "[INFO] merge phased gvcf"
    # time ${PARALLEL} --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_4_merge_vcf.log -j${THREADS} \
    # "${PYPY} ${CLAIR3_TRIO} MergeVcf_Trio \
    #     --output_fn ${GVCF_TMP_PATH}/{2}/phased_merge_{1}.vcf \
    #     --gvcf_fn ${GVCF_TMP_PATH}/{2}/phased_merge_{1}.gvcf \
    #     --platform ${PLATFORM} \
    #     --print_ref_calls ${SHOW_REF} \
    #     --gvcf ${GVCF} \
    #     --haploid_precise ${HAP_PRE} \
    #     --haploid_sensitive ${HAP_SEN} \
    #     --non_var_gvcf_fn {5}/non_var.gvcf \
    #     --ref_fn ${REFERENCE_FILE_PATH} \
    #     --merge_gvcf_only True \
    #     --ctgName {1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${INPUT_PILEUP_VCF[@]} :::+ ${TRIO_M_OUTPUT_VCF[@]} :::+ ${INPUT_PILEUP_GVCF_PATH[@]} |& tee ${LOG_PATH}/4_MV.log

    # ${PARALLEL}  -j${THREADS} \
    # ${PYPY} ${CLAIR3_TRIO} SortVcf_Trio \
    #         --input_dir ${GVCF_TMP_PATH}/{1} \
    #         --vcf_fn_prefix "phased_merge" \
    #         --vcf_fn_suffix ".gvcf" \
    #         --output_fn ${OUTPUT_FOLDER}/phased_{1}.gvcf \
    #         --sampleName {1} \
    #         --ref_fn ${REFERENCE_FILE_PATH} \
    #         --contigs_fn ${TMP_FILE_PATH}/CONTIGS ::: ${ALL_SAMPLE[@]}

    # fi

fi



if [ ${ENABLE_OUTPUT_HAPLOTAGGING} == True ]
then
    echo "[INFO] Haplotagging input BAM file using Whatshap, need some time to finish!"

    time ${PARALLEL} --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_5_haptag.log -j${THREADS} \
    "${WHATSHAP} haplotag \
        --output ${OUTPUT_FOLDER}/phased_{1}.bam \
        --reference ${REFERENCE_FILE_PATH} \
        --ignore-read-groups \
        ${OUTPUT_FOLDER}/phased_{1}.vcf.gz \
        {2}"  ::: ${ALL_SAMPLE[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]} |& tee ${LOG_PATH}/5_haplotag.log

    time ${PARALLEL} --retries ${RETRIES} --joblog ${LOG_PATH}/parallel_5_haptag.log -j${THREADS} \
    "${SAMTOOLS} index -@12 \
        ${OUTPUT_FOLDER}/phased_{1}.bam" ::: ${ALL_SAMPLE[@]}

    echo $''
    echo "[INFO] Finish Haplotagging, output haplotagged file [Child]: ${OUTPUT_FOLDER}/phased_${SAMPLE_C}.bam"
    echo "[INFO] Finish Haplotagging, output haplotagged file [Parent 1]: ${OUTPUT_FOLDER}/phased_${SAMPLE_P1}.bam"
    echo "[INFO] Finish Haplotagging, output haplotagged file [Parent 2]: ${OUTPUT_FOLDER}/phased_${SAMPLE_P2}.bam"
fi


echo $''
echo "[INFO] Finish calling, output VCF file [Child]: ${OUTPUT_FOLDER}/${SAMPLE_C}.vcf.gz"
echo "[INFO] Finish calling, output VCF file [Parent 1]: ${OUTPUT_FOLDER}/${SAMPLE_P1}.vcf.gz"
echo "[INFO] Finish calling, output VCF file [Parent 2]: ${OUTPUT_FOLDER}/${SAMPLE_P2}.vcf.gz"

if [ ${GVCF} == True ]
then
echo $''
echo "[INFO] Finish calling, output gVCF file [Child]: ${OUTPUT_FOLDER}/${SAMPLE_C}.gvcf.gz"
echo "[INFO] Finish calling, output gVCF file [Parent 1]: ${OUTPUT_FOLDER}/${SAMPLE_P1}.gvcf.gz"
echo "[INFO] Finish calling, output gVCF file [Parent 2]: ${OUTPUT_FOLDER}/${SAMPLE_P2}.gvcf.gz"
fi
