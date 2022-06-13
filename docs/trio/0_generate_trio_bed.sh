# generate trio bed from individual bed file
# e.g. input HG002, HG003, and HG004 bed, output HG002 trio bed

# please provide your BED file here
ALL_ORI_BED_FILE_PATH=(
"${CHILD_BED_FILE_PATH}"
"${P1_BED_FILE_PATH}"
"${P2_BED_FILE_PATH}"
)
# target output folder
DATASET_FOLDER_PATH="output"


# merge trio's bed file
BED2=${ALL_ORI_BED_FILE_PATH[0]}
BED3=${ALL_ORI_BED_FILE_PATH[1]}
BED4=${ALL_ORI_BED_FILE_PATH[2]}

_INPUT_DIR=${DATASET_FOLDER_PATH}/bed
mkdir -p ${_INPUT_DIR}
cp ${CHILD_BED_FILE_PATH} ${_INPUT_DIR}
cp ${P1_BED_FILE_PATH} ${_INPUT_DIR}
cp ${P2_BED_FILE_PATH} ${_INPUT_DIR}

_TRIO_BED_PATH=${DATASET_FOLDER_PATH}/bed/trio.bed

docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools intersect -a ${BED2} -b ${BED3} > ${_INPUT_DIR}/tmp_out
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools sort -i ${_INPUT_DIR}/tmp_out > ${_INPUT_DIR}/0203.bed
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools intersect -a ${_INPUT_DIR}/0203.bed -b ${BED4} > ${_INPUT_DIR}/tmp_out
docker run -v "${_INPUT_DIR}":"${_INPUT_DIR}" biocontainers/bedtools:v2.26.0dfsg-3-deb_cv1 bedtools sort -i ${_INPUT_DIR}/tmp_out > ${_TRIO_BED_PATH}
