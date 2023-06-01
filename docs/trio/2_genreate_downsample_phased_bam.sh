# please using this script to generate the downsampled of Representation Unification bam for each sample for training

# input unifold
tar_dir="XXX" # Representation Unification folder from https://github.com/HKU-BAL/Clair3-Trio/blob/trio/docs/trio/representation_unification_trio.md
_SAMPLE_N='HG002' # your sample name
BAM_FILE_PATH=${tar_dir}/merged.bam

SUBSAMPLED_BAMS_FOLDER_PATH="XXX" #output downsampling folder
mkdir -p ${SUBSAMPLED_BAMS_FOLDER_PATH}

# setting the downsamping fraction code for samtools -s
# for example, when the original data is in 120x coverage, you can using the code following the generate 10x to 80x data
# (you can check your data coverage via $mosdepth -t 5 -n -x --quantize 0:15:150: ${_SAMPLE_N}_phased ${_SAMPLE_N}_phased.bam)
# 120.33, downsample to 10x, 20x, etc.
DEPTHS_N=(10 20 30 40 50 60 70 80)
# FRAC code for samtools -s option, check samtools manpage
DEPTHS=(083 166 249 332 416 499 582 665)


# Other parameters
THREADS=8
PARALLEL='parallel'
SAMTOOLS='samtools'

# Subsample BAM to specific depths in parallel
${PARALLEL} -j${THREADS} "${SAMTOOLS} view -@12 -s 42.{1} -b -o ${SUBSAMPLED_BAMS_FOLDER_PATH}/${_SAMPLE_N}_{2}.bam ${BAM_FILE_PATH}"  ::: ${DEPTHS[@]}  :::+ ${DEPTHS_N[@]}
${PARALLEL} -j${THREADS} "${SAMTOOLS} index ${SUBSAMPLED_BAMS_FOLDER_PATH}/${_SAMPLE_N}_{1}.bam" ::: ${DEPTHS_N[@]}

# comfirm the downsampled coverage for each files
cd ${SUBSAMPLED_BAMS_FOLDER_PATH}
${PARALLEL} -j ${THREADS} "mosdepth -t 5 -n -x --quantize 0:15:150: ${_SAMPLE_N}_{1} ${SUBSAMPLED_BAMS_FOLDER_PATH}/${_SAMPLE_N}_{1}.bam " ::: ${DEPTHS_N[@]}

