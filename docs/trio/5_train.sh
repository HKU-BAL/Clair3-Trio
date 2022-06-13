# Training trio model
PYTHON3=python3
PLATFORM="ont"
# Clair3-Trio's path
CLAIR3_TRIO="XXX/clair3.py"      

# bins folder for training, obtained from 4_create_tensors.sh
ALL_BINS_FOLDER_PATH="XXX/build/ALL_1368/bins_for_train"

# training output folder
TRAIN_FOLDER_PREFIX="XXX"
# training name
TRAIN_N="ALL"
MODEL_FOLDER_PATH="${TRAIN_FOLDER_PREFIX}/train/${TRAIN_N}"                        
mkdir -p ${MODEL_FOLDER_PATH}
cd ${MODEL_FOLDER_PATH}

# training setting
BATCH_SIZE="800"  #training batch size, e.g. 800
add_indel_length=1
MODEL_ARC=NN
MODEL_ALS="Clair3_Trio_Out3"
IF_ADD_MCV_LOSS=0
MCVLOSS_ALPHA=0

# A single GPU is used for model training
export CUDA_VISIBLE_DEVICES="1"

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
