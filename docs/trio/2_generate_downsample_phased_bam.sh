# please using this script to generate the downsampled of Representation Unification bam for each sample for training
# Note that for training a Clair-Trio model, this process should **separately apply to All child, parent_1, and parent_2 samples**. 

# input unifold
tar_dir="XXX" # Representation Unification folder from https://github.com/HKU-BAL/Clair3-Trio/blob/trio/docs/trio/representation_unification_trio.md
BAM_FILE_PATH=${tar_dir}/merged.bam

_SAMPLE_N='HG002' # your sample name

SUBSAMPLED_BAMS_FOLDER_PATH="XXX" #output downsampling folder
mkdir -p ${SUBSAMPLED_BAMS_FOLDER_PATH}

# setting the downsamping fraction code for samtools -s
# for example, when the original data is in 120x coverage, you can using the code following the generate 10x to 80x data
# (you can check your data coverage via $mosdepth -t 5 -n -x --quantize 0:15:150: ${_SAMPLE_N}_phased ${_SAMPLE_N}_phased.bam)
# 120.33, downsample to 10x, 20x, etc.
# we recommand to downsample your date to 10x, 30x, 50x or 60x, and high coverage (or full-depth) for each of your sample for Clair3-Trio model training
# For high coverage dataset, like your data have 120x coverage, you can downsample your data to 10x, 30x, 60x, 80x
# For median coverage dataset, like your data have 65x coverage, you can downsample your data to 10x, 30x, 50x, 65x (full-depth)
DEPTHS_N=(10 30 50 80)
# FRAC code for samtools -s option, check samtools manpage
DEPTHS=(083 249 416 665)


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

