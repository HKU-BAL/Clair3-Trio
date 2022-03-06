# Clair3 full alignment parameters
REPO_NAME = "Clair3"
from itertools import accumulate

zstd='zstd'
default_optimizer = "Radam"
default_loss_function = "FocalLoss"
matrix_depth_dict = {'ont': 89, 'hifi': 55, 'ilmn': 55}
support_platform = {'ont'}
# min_af_dict = {'ont':0.15}
min_af_dict = {'ont':0.8}

# Full alignment input feature list
channel = (
'reference_base', 'alternative_base', 'mapping_quality', 'base_quality', 'strand_info', 'variant_type', 'insert_base',
'phasing_info')  # phasing info if add_phasing
channel_size = len(channel)
flankingBaseNum = 16
no_of_positions = 2 * flankingBaseNum + 1
input_shape = [matrix_depth_dict['hifi'], no_of_positions, channel_size]
ont_input_shape = [matrix_depth_dict['ont'], no_of_positions, channel_size]
label_shape = [21, 3, no_of_positions, no_of_positions]
label_size = sum(label_shape)
label_shape_cum = list(accumulate(label_shape))
expandReferenceRegion = 1000000
SAMTOOLS_VIEW_FILTER_FLAG = 2316
NORMALIZE_NUM = 100

# Realignment parameters
partition_size = 500000
realign_chunk_size = 5000
phasing_window_size = 100000
illumina_phasing_window_size = 10000
max_phasing_depth = 15
min_phasing_read_coverage = 2
split_region_size = 1000
extend_bp = 10

# Training hyperparameters
chunk_size = 200
trainBatchSize = 1000 # require mod(chunk_size) == 0
predictBatchSize = 200
initialLearningRate = 1e-3
l2RegularizationLambda = 1e-7
trainingDatasetPercentage = 0.9
maxEpoch = 30
OPERATION_SEED = None
RANDOM_SEED = 42



ont_input_shape_trio = [matrix_depth_dict['ont'] * 3, no_of_positions, channel_size]
label_shape_trio = label_shape * 3
label_size_trio = label_size * 3
label_shape_cum_trio = list(accumulate(label_shape_trio))

# make sure each alt info. length <= 5000
max_alt_info_length = 5000 * 3


padding_channel = (
'reference_base', 'alternative_base', 'mapping_quality', 'base_quality', 'strand_info', 'variant_type', 'insert_base',
'phasing_info', 'is_no_padding') 
padding_channel_size = len(padding_channel)
p_ont_input_shape = [matrix_depth_dict['ont'], no_of_positions, padding_channel_size]
p_ont_input_shape_trio = [matrix_depth_dict['ont'] * 3, no_of_positions, padding_channel_size]


# iinfo(min=-128, max=127, dtype=int8)
# padding_value_c = "100"
# padding_value_p1 = "100"
# padding_value_p2 = "100"

padding_value_c = "30"
padding_value_p1 = "60"
padding_value_p2 = "90"