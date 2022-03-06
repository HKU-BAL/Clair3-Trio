#!/bin/bash
SCRIPT_NAME=$(basename "$0")
Usage="Usage: ./${SCRIPT_NAME} --bam_fn_c=BAM --bam_fn_p1=BAM --bam_fn_p2=BAM --ref_fn=REF --output=OUTPUT_DIR --threads=THREADS --model_path_clair3=MODEL_PREFIX --model_path_clair3_trio=MODEL_PREFIX [--bed_fn=BED] [options]"
# INFO: whole calling workflow of clair3

set -e
ARGS=`getopt -o f:t:p:o:r::c::s::h::g \
-l bam_fn_c:,bam_fn_p1:,bam_fn_p2:,ref_fn:,threads:,model_path_clair3:,model_path_clair3_trio:,platform:,output:,\
bed_fn::,vcf_fn::,ctg_name::,sample_name_c::,sample_name_p1::,sample_name_p2::,help::,qual::,samtools::,python::,pypy::,parallel::,whatshap::,chunk_num::,chunk_size::,var_pct_full::,\
snp_min_af::,indel_min_af::,ref_pct_full::,pileup_only::,pileup_phasing::,fast_mode::,gvcf::,print_ref_calls::,haploid_precise::,haploid_sensitive::,include_all_ctgs::,\
no_phasing_for_fa::,pileup_model_prefix::,trio_model_prefix::,call_snp_only:: -n 'run_clair3_trio.sh' -- "$@"`

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

    -- ) shift; break; ;;
    -h|--help ) print_help_messages; break ;;
    * ) print_help_messages; exit 0 ;;
   esac
done



SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
CLAIR3_TRIO="${SHELL_FOLDER}/../clair3.py"



# if [ ${BED_FILE_PATH} = "EMPTY" ] ; then BED_FILE_PATH= ; fi
RETRIES=4

PILEUP_CHECKPOINT_PATH="${MODEL_PATH}/${PILEUP_PREFIX}"
FULL_ALIGNMENT_CHECKPOINT_PATH="${MODEL_PATH}/${FA_PREFIX}"
LOG_PATH="${OUTPUT_FOLDER}/log"
TMP_FILE_PATH="${OUTPUT_FOLDER}/tmp"
SPLIT_BED_PATH="${TMP_FILE_PATH}/trio_input/split_beds"
TENSOR_CANDIDATE_FOLDER_PATH="${TMP_FILE_PATH}/trio_input/tensor_can"
INDEL_PATH="${TMP_FILE_PATH}/trio_input/alt_info"
PILEUP_VCF_PATH="${TMP_FILE_PATH}/pileup_output"
GVCF_TMP_PATH="${TMP_FILE_PATH}/gvcf_tmp_output"
PHASE_OUTPUT_PATH="${TMP_FILE_PATH}/phase_output"
TRIO_ALIGNMENT_OUTPUT_PATH="${TMP_FILE_PATH}/trio_alignment_output"
PHASE_VCF_PATH="${PHASE_OUTPUT_PATH}/phase_vcf"
PHASE_BAM_PATH="${PHASE_OUTPUT_PATH}/phase_bam"
CANDIDATE_BED_PATH="${FULL_ALIGNMENT_OUTPUT_PATH}/candidate_bed"
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1


mkdir -p ${SPLIT_BED_PATH}
mkdir -p ${TENSOR_CANDIDATE_FOLDER_PATH}
mkdir -p ${INDEL_PATH}
mkdir -p ${LOG_PATH}



echo "[INFO] * Clir3-Trio pipeline start"
echo "[INFO] * 0/7 Check environment variables"

# ${PYTHON} ${CLAIR3_TRIO} CheckEnvs_Trio \
#     --bam_fn_c ${BAM_FILE_PATH_C} \
#     --bam_fn_p1 ${BAM_FILE_PATH_P1} \
#     --bam_fn_p2 ${BAM_FILE_PATH_P2} \
#     --bed_fn ${BED_FILE_PATH} \
#     --output_fn_prefix ${OUTPUT_FOLDER} \
#     --ref_fn ${REFERENCE_FILE_PATH} \
#     --vcf_fn ${VCF_FILE_PATH} \
#     --ctg_name ${CONTIGS} \
#     --chunk_num ${CHUNK_NUM} \
#     --chunk_size ${CHUNK_SIZE} \
#     --include_all_ctgs ${INCLUDE_ALL_CTGS} \
#     --threads ${THREADS} \
#     --python ${PYTHON} \
#     --pypy ${PYPY} \
#     --samtools ${SAMTOOLS} \
#     --whatshap ${WHATSHAP} \
#     --parallel ${PARALLEL} \
#     --qual ${QUAL} \
#     --sampleName_c ${SAMPLE_C} \
#     --sampleName_p1 ${SAMPLE_P1} \
#     --sampleName_p2 ${SAMPLE_P2} \
#     --var_pct_full ${PRO} \
#     --ref_pct_full ${REF_PRO} \
#     --snp_min_af ${SNP_AF} \
#     --indel_min_af ${INDEL_AF}


# readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
# if [ ${#CHR[@]} -eq 0 ]; then echo "[INFO] Exit in environment checking"; exit 0; fi
# THREADS_LOW=$((${THREADS}*3/4))
# if [[ ${THREADS_LOW} < 1 ]]; then THREADS_LOW=1; fi

# cd ${OUTPUT_FOLDER}

# export CUDA_VISIBLE_DEVICES=""
# echo "[INFO] * 1/7 Call variants using pileup model"



ALL_SAMPLE=(
${SAMPLE_C}
${SAMPLE_P1}
${SAMPLE_P2}
)

# ALL_UNPHASED_BAM_FILE_PATH=(
# ${BAM_FILE_PATH_C}
# ${BAM_FILE_PATH_P1}
# ${BAM_FILE_PATH_P2}
# )

# CLAIR3_THREADS=$((${THREADS}/3))
# if [[ ${CLAIR3_THREADS} < 1 ]]; then CLAIR3_THREADS=1; fi

# echo "trio sample" ${ALL_SAMPLE[@]}
# echo "trio bam" ${ALL_UNPHASED_BAM_FILE_PATH[@]}
# echo "pileup threads" ${CLAIR3_THREADS}

# echo "running pileup"


# ${PARALLEL} -j3 --joblog  ${LOG_PATH}/parallel_1_clair3_pileup.log \
# "${SHELL_FOLDER}/..//scripts/clair3.sh \
#     --bam_fn={2} \
#     --ref_fn=${REFERENCE_FILE_PATH} \
#     --threads=${CLAIR3_THREADS} \
#     --model_path=${MODEL_PATH_C3} \
#     --output=${PILEUP_VCF_PATH}/{1} \
#     --platform="ont" \
#     --bed_fn=${BED_FILE_PATH} \
#     --vcf_fn=${VCF_FILE_PATH} \
#     --ctg_name=${CONTIGS} \
#     --sample_name={1} \
#     --chunk_num=${CHUNK_NUM} \
#     --chunk_size=${CHUNK_SIZE} \
#     --samtools=${SAMTOOLS} \
#     --python=${PYTHON} \
#     --pypy=${PYPY} \
#     --parallel=${PARALLEL} \
#     --whatshap=${WHATSHAP} \
#     --qual=${QUAL} \
#     --var_pct_full=${PRO} \
#     --ref_pct_full=${REF_PRO} \
#     --snp_min_af=${SNP_AF} \
#     --indel_min_af=${INDEL_AF} \
#     --pileup_only=True \
#     --gvcf=${GVCF} \
#     --fast_mode=${FAST_MODE} \
#     --call_snp_only=${SNP_ONLY} \
#     --print_ref_calls=${SHOW_REF} \
#     --haploid_precise=${HAP_PRE} \
#     --haploid_sensitive=${HAP_SEN} \
#     --include_all_ctgs=${INCLUDE_ALL_CTGS} \
#     --no_phasing_for_fa=${NO_PHASING} \
#     --pileup_model_prefix=${PILEUP_PREFIX} \
#     --fa_model_prefix=full_alignment" ::: ${ALL_SAMPLE[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]} |& tee ${LOG_PATH}/1_call_var_bam_pileup.log


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
echo ${CHR[@]}

chunk_num=20
CHUNK_LIST=`seq 1 ${chunk_num}`

# # select candidate from trio input
# echo "[INFO] * 2/7 Select Trio Candidates"
# time ${PARALLEL} --joblog ${LOG_PATH}/parallel_2_fiter_hete_snp_pileup.log -j${THREADS} \
# "${PYPY} ${CLAIR3_TRIO} SelectHetSnp_Trio \
# --alt_fn_c ${INPUT_PILEUP_VCF[0]} \
# --alt_fn_p1 ${INPUT_PILEUP_VCF[1]} \
# --alt_fn_p2 ${INPUT_PILEUP_VCF[2]} \
# --split_folder ${SPLIT_BED_PATH} \
# --sampleName ${TRIO_N} \
# --ref_pct_full 0.2 \
# --var_pct_full 1.0 \
# --chunk_num ${chunk_num} \
# --chr_prefix '' \
# --ctgName {1}" ::: ${CHR[@]} |& tee ${LOG_PATH}/2_FHSP.log



# # note that in training included the add_no_phasing_data_training
# # no using the giab bed files
# echo "[INFO] * 3/7 Creating Tensors"
# time ${PARALLEL}  --joblog ${LOG_PATH}/parallel_3_create_tensor.log -j${THREADS} \
# "${PYPY} ${CLAIR3_TRIO} CreateTensorFullAlignment \
# --bam_fn {3}/{1}.bam \
# --ref_fn ${REFERENCE_FILE_PATH} \
# --tensor_can_fn ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_{2}_{1}_{4} \
# --indel_fn ${INDEL_PATH}/{2}_{1}_{4} \
# --ctgName {1} \
# --samtools ${SAMTOOLS} \
# --platform ${PLATFORM} \
# --full_aln_regions ${SPLIT_BED_PATH}/${TRIO_N}_1000_{1}_{4} \
# --phasing_info_in_bam" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${ALL_PHASED_BAM_FILE_DIR[@]} ::: ${CHUNK_LIST[@]} |& tee ${LOG_PATH}/3_CT.log


# echo "[INFO] * 4/7 Merging Trio Tensors"
# time ${PARALLEL} --joblog ${LOG_PATH}/parallel_4_merge_tensors.log -j${THREADS} \
# "${PYTHON} ${CLAIR3_TRIO} Merge_Tensors_Trio \
# --tensor_fn_c ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${ALL_SAMPLE[0]}_{1}_{2} \
# --tensor_fn_p1 ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${ALL_SAMPLE[1]}_{1}_{2} \
# --tensor_fn_p2 ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${ALL_SAMPLE[2]}_{1}_{2} \
# --candidate_fn_c ${INDEL_PATH}/${ALL_SAMPLE[0]}_{1}_{2} \
# --candidate_fn_p1 ${INDEL_PATH}/${ALL_SAMPLE[1]}_{1}_{2} \
# --candidate_fn_p2 ${INDEL_PATH}/${ALL_SAMPLE[2]}_{1}_{2} \
# --tensor_fn ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${TRIO_N}_{1}_{2} \
# --candidate_fn ${INDEL_PATH}/${TRIO_N}_{1}_{2} \
# " ::: ${CHR[@]} ::: ${CHUNK_LIST[@]} |& tee ${LOG_PATH}/4_MT.log





# BINS_FOLDER_PATH="${TMP_FILE_PATH}/trio_input/bins"
# mkdir -p ${BINS_FOLDER_PATH}

# cp ${BAM_FILE_PATH_C}/CONTIGS ${TMP_FILE_PATH}
readarray -t CHR < "${BAM_FILE_PATH_C}/CONTIGS"
BINS_FOLDER_PATH="${BAM_FILE_PATH_C}/trio_input/bins"

echo $BINS_FOLDER_PATH
# exit 0;




# create tensors bin chunk number for each chr
bin_chunk_num=1
BIN_CHUNK_LIST=`seq 1 ${bin_chunk_num}`


# echo "[INFO] * 5/7 Coverting Tensors to Bins"
# time ${PARALLEL} --joblog ${LOG_PATH}/parallel_5_tensor2Bin.log -j${THREADS} \
# "${PYTHON} ${CLAIR3_TRIO} Tensor2Bin_Trio \
# --tensor_fn ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${TRIO_N}_{1} \
# --bin_fn ${BINS_FOLDER_PATH}/${TRIO_N}_{1}_{2} \
# --chunk_id {2} \
# --chunk_num ${bin_chunk_num} \
# --platform ${PLATFORM} \
# " ::: ${CHR[@]} ::: ${BIN_CHUNK_LIST[@]} |& tee ${LOG_PATH}/5_T2B.log

call_chunk_num=6
CALL_CHUNK_LIST=`seq 1 ${call_chunk_num}`
PREDICT_THREADS=6
use_gpu=1
add_indel_length=1
MODEL_ALS="Clair3_Trio_Out3"
MODEL_ARC=NN
export CUDA_VISIBLE_DEVICES="0"


CALL_PATH=${TMP_FILE_PATH}/predict_tensors
mkdir -p ${CALL_PATH}

# PYTHON=/autofs/bal31/jhsu/home/_env/anaconda3/envs/py38_3090/bin/python


echo "[INFO] * 6/7 Calling Trio Variants"
time ${PARALLEL} --joblog ${LOG_PATH}/parallel_6_predict.log -j${PREDICT_THREADS} \
"${PYTHON} ${CLAIR3_TRIO} CallVariants_Trio \
--tensor_fn ${BINS_FOLDER_PATH}/${TRIO_N}_{1}_{2} \
--chkpnt_fn ${MODEL_PATH_C3T}/${TRIO_PREFIX} \
--predict_fn ${CALL_PATH}/predict_${TRIO_N}_{1}_{2}_{3} \
--sampleName ${TRIO_N} \
--chunk_id {3} \
--chunk_num ${call_chunk_num} \
--is_from_tables \
--use_gpu ${use_gpu} \
--add_indel_length ${add_indel_length} \
--model_arc ${MODEL_ARC} \
--model_cls ${MODEL_ALS} \
--platform ${PLATFORM} \
--output_probabilities" ::: ${CHR[@]} ::: ${BIN_CHUNK_LIST[@]} ::: ${CALL_CHUNK_LIST[@]} |& tee ${LOG_PATH}/6_PREDICT.log

ALL_SAMPLE_IDX=(
0
1
2
)

time ${PARALLEL}  --joblog ${LOG_PATH}/parallel_7_call.log  -j ${THREADS} \
"${PYTHON} ${CLAIR3_TRIO} CallVariants_Trio \
--tensor_fn ${CALL_PATH}/predict_${TRIO_N}_{1}_{4}_{5} \
--chkpnt_fn 0 \
--call_fn ${CALL_PATH}/{2}_{1}_{4}_{5}.vcf \
--sampleName {2} \
--ref_fn ${REFERENCE_FILE_PATH} \
--add_indel_length ${add_indel_length} \
--platform ${PLATFORM} \
--model_arc ${MODEL_ARC} \
--trio_n_id {3} \
--showRef \
--input_probabilities" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${ALL_SAMPLE_IDX[@]} ::: ${BIN_CHUNK_LIST[@]} ::: ${CALL_CHUNK_LIST[@]} |& tee ${LOG_PATH}/7_CV.log

time ${PARALLEL}  --joblog ${LOG_PATH}/parallel_8_sort.log  -j ${THREADS} \
"${PYTHON} ${CLAIR3_TRIO} SortVcf \
--input_dir ${CALL_PATH} \
--vcf_fn_prefix {1} \
--output_fn ${OUTPUT_FOLDER}/{1}.vcf \
--sampleName {1} \
--ref_fn ${REFERENCE_FILE_PATH} \
--contigs_fn ${TMP_FILE_PATH}/CONTIGS" ::: ${ALL_SAMPLE[@]}  |& tee ${LOG_PATH}/8_SORT.log




echo $''
echo "[INFO] Finish calling, output folder: ${OUTPUT_FOLDER}"



# call_chunk_num=1
# CALL_CHUNK_LIST=`seq 1 ${call_chunk_num}`





# CALL_PATH="${TMP_FILE_PATH}/trio_input/call"
# mkdir -p ${CALL_PATH}

# cd ${CALL_PATH}
# echo "[INFO] Call var parallel for testing"



# MODEL_ALS="Clair3_Trio_Out3"
# MODEL_ARC=NN
# use_gpu=0
# add_indel_length=1
# export CUDA_VISIBLE_DEVICES="1"


# # 10G GPU: maximum 6 threads 
# PREDICT_THREADS=6

# time ${PARALLEL} --joblog ${LOG_PATH}/parallel_5_call_v.txt -j${PREDICT_THREADS} \
# echo "${PYTHON} ${CLAIR3_TRIO} CallVariants_Trio \
# --tensor_fn ${TENSOR_CANDIDATE_FOLDER_PATH}/tensor_can_${TRIO_N}_{1}_1 \
# --chkpnt_fn ${MODEL_PATH_C3T}/${TRIO_PREFIX} \
# --call_fn ${CALL_PATH}/predict_tensors/predict_${TRIO_N}.vcf \
# --predict_fn ${CALL_PATH}/predict_tensors/predict_${TRIO_N}_{1}_{3}_{4} \
# --ref_fn ${REFERENCE_FILE_PATH} \
# --sampleName ${TRIO_N} \
# --chunk_id {4} \
# --ctgName {1} \
# --chunk_num ${call_chunk_num} \
# --use_gpu ${use_gpu} \
# --add_indel_length ${add_indel_length} \
# --model_arc ${MODEL_ARC} \
# --model_cls ${MODEL_ALS} \
# --platform ${PLATFORM} \
# " ::: ${CHR[@]} ::: ${EPOCH[@]} ::: ${BIN_CHUNK_LIST[@]} ::: ${CALL_CHUNK_LIST[@]} |& tee ${CALL_PATH}/PREDICT.log


















# THREADS
# CLAIR3
# CLAIR3_THREADS
# CLAIR3_MODEL
# REFERENCE_FILE_PATH
# BED_FILE_PATH

# CLAIR3_THREADS=$((${THREADS}/3))
# if [[ ${CLAIR3_THREADS} < 1 ]]; then CLAIR3_THREADS=1; fi


# PARALLEL
# LOG_PATH
# PILEUP_OUTPUT_PATH
# # Step 1 run Clair3 pileup



# ALL_SAMPLE=(
# )

# ALL_UNPHASED_BAM_FILE_PATH=(
# )


# ${LOG_PATH}/parallel_1_pileup_c.log

# # Call variants using Clair3‘s pileup model with the --pileup_only option
# # Only select the candidates in the high-confident BED regions for model training (with --bed_fn)
# ${PARALLEL} -j3 --joblog  ${LOG_PATH}/parallel_1_clair3_pileup.log \
# "${CLAIR3}/run_clair3.sh \
# --bam_fn={3} \
# --ref_fn=${REFERENCE_FILE_PATH} \
# --threads=${CLAIR3_THREADS} \
# --platform="ont" \
# --model_path="${CLAIR3_MODEL}" \
# --output=${PILEUP_OUTPUT_PATH}/{1} \
# --bed_fn=${ALL_BED_FILE_PATH} \
# --bed_fn ${BED_FILE_PATH} \
# --ctg_name ${CONTIGS}  \
# --include_all_ctgs ${INCLUDE_ALL_CTGS} \
# --python ${PYTHON} \
# --pypy ${PYPY} \
# --samtools ${SAMTOOLS} \
# --whatshap ${WHATSHAP} \
# --parallel ${PARALLEL} \
# --qual ${QUAL} \
# --sampleName ${SAMPLE} \
# --var_pct_full ${PRO} \
# --ref_pct_full ${REF_PRO} \
# --snp_min_af ${SNP_AF} \
# --indel_min_af ${INDEL_AF}
# --pileup_phasing" ::: ${ALL_SAMPLE[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]}

