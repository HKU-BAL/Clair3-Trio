import logging
import random
import numpy as np
from argparse import ArgumentParser, SUPPRESS
import tensorflow_addons as tfa
import tensorflow as tf
import tables
import os
import sys
from itertools import accumulate
from time import time

import trio.model as model_path
from shared.utils import str2bool



logging.basicConfig(format='%(message)s', level=logging.INFO)
tables.set_blosc_max_threads(512)
os.environ['NUMEXPR_MAX_THREADS'] = '1028'
os.environ['NUMEXPR_NUM_THREADS'] = '16'


def get_label_task(label, label_shape_cum, task):
    if task == 0:
        return label[:label_shape_cum[task]]
    elif task == len(label_shape_cum) - 1:
        return label[label_shape_cum[task - 1]:]
    else:
        return label[label_shape_cum[task - 1]:label_shape_cum[task]]


def cal_class_weight(samples_per_cls, no_of_classes, beta=0.999):
    effective_num = 1.0 - np.power(beta, samples_per_cls)
    cls_weights = (1.0 - beta) / np.array(effective_num)
    cls_weights = cls_weights / np.sum(cls_weights) * no_of_classes
    return cls_weights



class FocalLoss(tf.keras.losses.Loss):
    """
    updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
    more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
    """

    def __init__(self, label_shape_cum, task, effective_label_num=None, gamma=2):
        super(FocalLoss, self).__init__()
        self.gamma = gamma
        self.cls_weights = None
        if effective_label_num is not None:
            task_label_num = get_label_task(effective_label_num, label_shape_cum, task)
            cls_weights = cal_class_weight(task_label_num, len(task_label_num))
            cls_weights = tf.constant(cls_weights, dtype=tf.float32)
            cls_weights = tf.expand_dims(cls_weights, axis=0)
            self.cls_weights = cls_weights

    def call(self, y_true, y_pred):
        y_pred = tf.clip_by_value(y_pred, clip_value_min=1e-9, clip_value_max=1 - 1e-9)
        cross_entropy = -y_true * tf.math.log(y_pred)
        # legacy
        # weight = ((1 - y_pred) ** self.gamma) * y_true
        weight = ((1 - y_pred) ** self.gamma)
        FCLoss = cross_entropy * weight
        if self.cls_weights is not None:
            FCLoss = FCLoss * self.cls_weights
        reduce_fl = tf.reduce_sum(FCLoss, axis=-1)
        # print(reduce_fl.shape)
        # print(int(sum(reduce_fl)))
        # if int(sum(reduce_fl)) > 92536321156295163904:
        #     import pdb; pdb.set_trace()
        return reduce_fl



GT21_LABELS = [
    'AA',
    'AC',
    'AG',
    'AT',
    'CC',
    'CG',
    'CT',
    'GG',
    'GT',
    'TT',
    'DelDel',
    'ADel',
    'CDel',
    'GDel',
    'TDel',
    'InsIns',
    'AIns',
    'CIns',
    'GIns',
    'TIns',
    'InsDel'
]

def get_allel_from_gt(gt):
    if len(gt) <= 4:
        return gt[0], gt[1]
    return gt[0], gt[3]

def check_if_MC(gt1, gt2, gt3):
    all_c_gt = [''.join(sorted([i, j])) for i in get_allel_from_gt(gt2) for j in get_allel_from_gt(gt3)]
    tar_c_gt = ''.join(sorted(get_allel_from_gt(gt1)))
    # print(gt1, gt2, gt3, all_c_gt, tar_c_gt)
    if tar_c_gt in all_c_gt:
        return 2
    if tar_c_gt[0] in ''.join(all_c_gt):
        return 1
    if tar_c_gt[1] in ''.join(all_c_gt):
        return 1
    return 0



class MCVLoss(tf.keras.losses.Loss):
    """
    updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
    more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
    """

    def __init__(self, alpha=1, gt_error_rate=1e-8, single_tensor_size=21):
        super(MCVLoss, self).__init__()
        self.alpha = alpha
        self.gt_error_rate = gt_error_rate
        self.single_tensor_size = single_tensor_size

        # print("MCVLoss, alpha is %.2f" % (self.alpha))
        tmp_df = np.zeros(21 ** 3)
        for i in range(21):
            for j in range(21):
                for k in range(21):
                    gt_1, gt_2, gt_3 = GT21_LABELS[i], GT21_LABELS[j], GT21_LABELS[k]
                    tar_mc_type = check_if_MC(gt_1, gt_2, gt_3)
                    # print(tar_mc_type)
                    if tar_mc_type == 2:
                        tmp_df[k + j * 21 + i * 21 * 21] = 1
                    if tar_mc_type == 1:
                        tmp_df[k + j * 21 + i * 21 * 21] = gt_error_rate
                        # import pdb; pdb.set_trace()
        self.cls_weights = tf.constant(tmp_df , dtype=tf.float32)
        # import pdb; pdb.set_trace()

        
    def call(self, y_true, y_pred):
        y_pred_1 = y_pred[:, :self.single_tensor_size][:, :21]
        y_pred_2 = y_pred[:, self.single_tensor_size : (self.single_tensor_size * 2)][:, :21]
        y_pred_3 = y_pred[:, (self.single_tensor_size * 2): ][:, :21]
        # import pdb; pdb.set_trace()

        y_pred_trio = tf.matmul(tf.expand_dims(y_pred_1, axis=-1), tf.expand_dims(y_pred_2, axis=-2))
        # y_pred_trio = tf.reshape(y_pred_trio, (y_pred_trio.shape[0], -1))
        y_pred_trio = tf.reshape(y_pred_trio, (tf.shape(y_pred_trio)[0], 21 * 21))
        y_pred_trio = tf.matmul(tf.expand_dims(y_pred_trio, axis=-1), tf.expand_dims(y_pred_3, axis=-2))
        # y_pred_trio = tf.reshape(y_pred_trio, (y_pred_trio.shape[0], -1))
        # specific tensor shape to increase GPU utilize
        y_pred_trio = tf.reshape(y_pred_trio, (tf.shape(y_pred_trio)[0], 21 * 21 * 21))

        y_pred_trio_w  = y_pred_trio * self.cls_weights
        reduce_ls = tf.reduce_sum(y_pred_trio_w, axis=-1)
        reduce_ls = tf.clip_by_value(reduce_ls, clip_value_min=1e-9, clip_value_max=1 - 1e-9)
        reduce_ls = -self.alpha * tf.math.log(reduce_ls)
        # reduce_ls = (1 - reduce_ls) ** 2

        # import pdb; pdb.set_trace()
        return reduce_ls


def get_chunk_list(chunk_offset, train_data_size, chunk_size):
    """
    get chunk list for training and validation data. we will randomly split training and validation dataset,
    all training data is directly acquired from various tensor bin files.

    """
    all_shuffle_chunk_list = []
    total_size = 0
    offset_idx = 0
    for bin_idx, chunk_num in enumerate(chunk_offset):
        all_shuffle_chunk_list += [(bin_idx, chunk_idx) for chunk_idx in range(chunk_num)]
    np.random.seed(0)
    np.random.shuffle(all_shuffle_chunk_list)  # keep the same random validate dataset
    for bin_idx, chunk_num in enumerate(chunk_offset):
        if chunk_num * chunk_size + total_size >= train_data_size:
            chunk_num = (train_data_size - total_size) // chunk_size
            offset_idx += chunk_num
            return np.array(all_shuffle_chunk_list[:offset_idx]), np.array(all_shuffle_chunk_list[offset_idx + 1:])
        else:
            total_size += chunk_num * chunk_size
            offset_idx += chunk_num


def exist_file_prefix(exclude_training_samples, f):
    for prefix in exclude_training_samples:
        if prefix in f:
            return True
    return False

# basic 3 to 1 model (N to 1)
def train_model_N1(args):
    platform = args.platform
    pileup = args.pileup
    add_indel_length = args.add_indel_length
    exclude_training_samples = args.exclude_training_samples
    exclude_training_samples = set(exclude_training_samples.split(',')) if exclude_training_samples else set()
    add_validation_dataset = args.validation_dataset
    ochk_prefix = args.ochk_prefix if args.ochk_prefix is not None else ""
    
    training_start_time = time()
    import trio.param_t as param
    model = model_path.Clair3_Trio_Basic(add_indel_length=add_indel_length, is_padding=args.add_padding)

    tar_model_type = args.tar_model_type
    _is_reverse_23 = args.add_reverse_23
    tensor_shape = param.ont_input_shape_trio

    if args.add_padding:
        tensor_shape = param.p_ont_input_shape_trio
    label_size, label_shape, label_shape_cum = param.label_size, param.label_shape, param.label_shape_cum
    label_size_trio = param.label_size_trio
    # 2000, 200
    batch_size, chunk_size = param.trainBatchSize, param.chunk_size
    if args.batch_size != 0:
        batch_size = args.batch_size
    random.seed(param.RANDOM_SEED)
    np.random.seed(param.RANDOM_SEED)
    learning_rate = args.learning_rate if args.learning_rate else param.initialLearningRate
    max_epoch = args.maxEpoch if args.maxEpoch else param.maxEpoch

    task_num = 4 if add_indel_length else 2
    TensorShape = (tf.TensorShape([None] + tensor_shape), \
        tuple(tf.TensorShape([None, label_shape[task]]) for task in range(task_num)))
    TensorDtype = (tf.int32, tuple(tf.float32 for _ in range(task_num)))

    bin_list = os.listdir(args.bin_fn)
    if args.is_debuging:
        bin_list = bin_list[:1]
        max_epoch = 2
    # default we exclude sample hg003 and all chr20 for training
    bin_list = [f for f in bin_list if '_20_' not in f and not exist_file_prefix(exclude_training_samples, f)]



    logging.info("[INFO] total {} training bin files: {}".format(len(bin_list), ','.join(bin_list)))
    total_data_size = 0
    table_dataset_list = []
    validate_table_dataset_list = []
    chunk_offset = np.zeros(len(bin_list), dtype=int)

    # import pdb; pdb.set_trace()
    effective_label_num = None
    for bin_idx, bin_file in enumerate(bin_list):
        table_dataset = tables.open_file(os.path.join(args.bin_fn, bin_file), 'r')
        validate_table_dataset = tables.open_file(os.path.join(args.bin_fn, bin_file), 'r')
        table_dataset_list.append(table_dataset)
        validate_table_dataset_list.append(validate_table_dataset)

        chunk_num = (len(table_dataset.root.label) - batch_size) // chunk_size
        data_size = int(chunk_num * chunk_size)
        chunk_offset[bin_idx] = chunk_num
        total_data_size += data_size

    if args.is_debuging:
        total_data_size = int(total_data_size / 10)

    train_data_size = total_data_size * param.trainingDatasetPercentage
    validate_data_size = int((total_data_size - train_data_size) // chunk_size) * chunk_size
    train_data_size = int(train_data_size // chunk_size) * chunk_size
    train_shuffle_chunk_list, validate_shuffle_chunk_list = get_chunk_list(chunk_offset, train_data_size, chunk_size)
    # import pdb; pdb.set_trace()

    # # for _1, _t = 10000, 143614
    # _t = 10000
    # _t = 33 # _2
    # _t = 33 # _2
    # table_dataset_list[0].root.position_matrix[_t].shape
    # print(table_dataset_list[0].root.position[_t], table_dataset_list[0].root.alt_info[_t], table_dataset_list[0].root.label[_t])
    # table_dataset_list[0].root.position[-1]
    # _s = 1
    # for _d in range(0, 89):
    #     print(_d, table_dataset_list[0].root.position_matrix[_t][_d + _s * 89,:,:].reshape(-1)[:20])
    # len(table_dataset_list[0].root.label)
    # len(table_dataset_list[0].root.position)

    # print(table_dataset_list[0].root.position[_t], table_dataset_list[0].root.alt_info[_t], table_dataset_list[0].root.label[_t])

    def DataGenerator(x, data_size, shuffle_chunk_list, train_flag=True):

        """
        data generator for pileup or full alignment data processing, pytables with blosc:lz4hc are used for extreme fast
        compression and decompression. random chunk shuffling and random start position to increase training model robustness.

        """
        _idx = 0
        chunk_iters = batch_size // chunk_size
        batch_num = data_size // batch_size
        position_matrix = np.empty([batch_size] + tensor_shape, np.int32)
        label = np.empty((batch_size, label_size_trio), np.float32)

        random_start_position = np.random.randint(0, batch_size) if train_flag else 0
        if train_flag:
            np.random.shuffle(shuffle_chunk_list)

        _tar_lebel_id = 0 if tar_model_type == 'child' else 1
        # import pdb; pdb.set_trace()

        # if _is_reverse_23, loop two time
        _tmp_loop_i = 2 if _is_reverse_23 else 1
        for _ in range(_tmp_loop_i):
            for batch_idx in range(batch_num):
                for chunk_idx in range(chunk_iters):
                    offset_chunk_id = shuffle_chunk_list[batch_idx * chunk_iters + chunk_idx]
                    bin_id, chunk_id = offset_chunk_id
                    position_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.position_matrix[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    label[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.label[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    # for _ii, _t in enumerate(range(random_start_position + chunk_id * chunk_size, random_start_position + (chunk_id + 1) * chunk_size)):
                    #     print(_ii, _t)
                    #     print(x[bin_id].root.position[_t], x[bin_id].root.alt_info[_t], x[bin_id].root.label[_t])
                    #     for _i in range(3):
                    #         print(', '.join([str(i) for i in x[bin_id].root.label[_t].reshape(3, -1)[_i, :]]))

                    # import pdb; pdb.set_trace()

                _tar_label_idx = _tar_lebel_id * label_size
                # import pdb; pdb.set_trace()
                # print(len(position_matrix), len(label))

                # # adhoc rescure data
                # tmp_l = position_matrix[:,:,:,-1:]
                # tmp_l[tmp_l == 100] = 90
                # tmp_l[tmp_l == -56] = 60
                # tmp_l[tmp_l == 44] = 30

                # # tmp_l[tmp_l == 100] = 100
                # # tmp_l[tmp_l == -56] = 100
                # # tmp_l[tmp_l == 44] = 100
                # position_matrix[:,:,:,-1:] = tmp_l


                # import pdb; pdb.set_trace()
                yield position_matrix, (
                            label[:, _tar_label_idx + 0                 : _tar_label_idx + label_shape_cum[0]],
                            label[:, _tar_label_idx + label_shape_cum[0]: _tar_label_idx + label_shape_cum[1]],
                            label[:, _tar_label_idx + label_shape_cum[1]: _tar_label_idx + label_shape_cum[2]],
                            label[:, _tar_label_idx + label_shape_cum[2]: _tar_label_idx + label_size]
                        )
                # if _idx > 10:
                #     break
                # _idx += 1

        # two rounds for generating reversed data
        # _is_reverse_23 = False
        if _is_reverse_23:
            one_tensor_shape = int(tensor_shape[0] / 3)
            for batch_idx in range(batch_num):
                for chunk_idx in range(chunk_iters):
                    offset_chunk_id = shuffle_chunk_list[batch_idx * chunk_iters + chunk_idx]
                    bin_id, chunk_id = offset_chunk_id
                    position_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.position_matrix[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    label[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.label[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    # for _ii, _t in enumerate(range(random_start_position + chunk_id * chunk_size, random_start_position + (chunk_id + 1) * chunk_size)):
                    #     print(_ii, _t)
                    #     print(x[bin_id].root.position[_t], x[bin_id].root.alt_info[_t], x[bin_id].root.label[_t])
                    #     for _i in range(3):
                    #         print(', '.join([str(i) for i in x[bin_id].root.label[_t].reshape(3, -1)[_i, :]]))

                if _is_reverse_23:
                    position_matrix = np.concatenate((
                    position_matrix[:, one_tensor_shape * 0 : one_tensor_shape * 1,:,:], \
                    position_matrix[:, one_tensor_shape * 2 : one_tensor_shape * 3:,:], \
                    position_matrix[:, one_tensor_shape * 1 : one_tensor_shape * 2,:,:]), axis=1)

                    label = np.concatenate((
                    label[:, label_size * 0 : label_size * 1], \
                    label[:, label_size * 2 : label_size * 3], \
                    label[:, label_size * 1 : label_size * 2]), axis=1)

                _tar_label_idx = _tar_lebel_id * label_size

                yield position_matrix, (
                            label[:, _tar_label_idx + 0                 : _tar_label_idx + label_shape_cum[0]],
                            label[:, _tar_label_idx + label_shape_cum[0]: _tar_label_idx + label_shape_cum[1]],
                            label[:, _tar_label_idx + label_shape_cum[1]: _tar_label_idx + label_shape_cum[2]],
                            label[:, _tar_label_idx + label_shape_cum[2]: _tar_label_idx + label_size]
                        )
                


    train_dataset = tf.data.Dataset.from_generator(
        lambda: DataGenerator(table_dataset_list, train_data_size, train_shuffle_chunk_list, True), TensorDtype,
        TensorShape).prefetch(buffer_size=tf.data.experimental.AUTOTUNE)
    validate_dataset = tf.data.Dataset.from_generator(
        lambda: DataGenerator(validate_table_dataset_list, validate_data_size, validate_shuffle_chunk_list, False), TensorDtype,
        TensorShape).prefetch(buffer_size=tf.data.experimental.AUTOTUNE)

    total_steps = max_epoch * train_data_size // batch_size

    # import pdb; pdb.set_trace()
    # for example in train_dataset:
    #     print(example)
    # import pdb; pdb.set_trace()

    #RectifiedAdam with warmup start
    optimizer = tfa.optimizers.Lookahead(tfa.optimizers.RectifiedAdam(
            lr=learning_rate,
            total_steps=total_steps,
            warmup_proportion=0.1,
            min_lr=learning_rate*0.75,
        ))

    # policy = tf.keras.mixed_precision.experimental.Policy('mixed_float16')
    # tf.keras.mixed_precision.experimental.set_policy(policy)
    # print('Compute dtype: %s' % policy.compute_dtype)
    # print('Variable dtype: %s' % policy.variable_dtype)
    # optimizer = tf.keras.mixed_precision.experimental.LossScaleOptimizer(optimizer, "dynamic")

    # import pdb; pdb.set_trace()

    loss_func = [FocalLoss(label_shape_cum, task, effective_label_num) for task in range(task_num)]
    loss_task = {"output_{}".format(task + 1): loss_func[task] for task in range(task_num)}
    metrics = {"output_{}".format(task + 1): tfa.metrics.F1Score(num_classes=label_shape[task], average='micro') for
               task in range(task_num)}

    model.compile(
        loss=loss_task,
        metrics=metrics,
        optimizer=optimizer,
        # run_eagerly=True
    )
    early_stop_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10, mode="min")
    model_save_callbakck = tf.keras.callbacks.ModelCheckpoint(ochk_prefix + ".{epoch:02d}", period=1, save_weights_only=False)

    # Use first 20 element to initialize tensorflow model using graph mode
    # output = model(np.array(table_dataset_list[0].root.position_matrix[:20]))
    # print(np.array(table_dataset_list[0].root.position_matrix[:20]).shape)

    model.build(tf.TensorShape([None] + tensor_shape))
    logging.info(model.summary(print_fn=logging.info))


    #import pdb; pdb.set_trace()
    # from tensorflow.keras.utils import plot_model
    # plot_model(model, to_file=ochk_prefix + 'model.png', show_layer_names=True, show_shapes=True)

    logging.info("[INFO] The size of dataset: {}".format(total_data_size))
    logging.info("[INFO] The training batch size: {}".format(batch_size))
    logging.info("[INFO] The training learning_rate: {}".format(learning_rate))
    logging.info("[INFO] Total training steps: {}".format(total_steps))
    logging.info("[INFO] Maximum training epoch: {}".format(max_epoch))
    logging.info("[INFO] Start training...")

    validate_dataset = validate_dataset if add_validation_dataset else None
    if args.chkpnt_fn is not None:
        model.load_weights(args.chkpnt_fn)
        logging.info("pretrained model at: %s" % (args.chkpnt_fn))

    # import wandb
    # from wandb.keras import WandbCallback
    # wandb.login(key='43ff292c29b48569c32e92bbf82356581cc71cae')

    # _train_name = "trio-" + ochk_prefix.split('/')[-2]
    # #import datetime
    # #log_dir = ochk_prefix  + "N" + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    # #wandb.tensorboard.patch(root_logdir=log_dir)
    #             #sync_tensorboard=False,
    # wandb.init(project='clair3-trio', entity='junhao', name=_train_name,
    #             config={
    #               "learning_rate": learning_rate,
    #               "epochs": max_epoch,
    #               "batch_size": batch_size,
    #               "model_s_type": "N1",
    #               "output_tensor_shape": label_size,
    #               "tar_model_type": tar_model_type,
    #               "add_reverse_23": args.add_reverse_23,
    #            })

    
    # # import pdb; pdb.set_trace()
    # #tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1, write_graph=True, write_images=True)


    

    train_history = model.fit(x=train_dataset,
                              epochs=max_epoch,
                              validation_data=validate_dataset,
                              callbacks=[early_stop_callback, model_save_callbakck],
                              verbose=1,
                              shuffle=False)

    # import datetime
    # log_dir = ochk_prefix  + datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    # tensorboard_callback = tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1, write_graph=True, write_images=True)
    # train_history = model.fit(x=train_dataset,
    #                           epochs=max_epoch,
    #                           validation_data=validate_dataset,
    #                           callbacks=[tensorboard_callback, early_stop_callback, model_save_callbakck],
    #                           verbose=1,
    #                           shuffle=False)


    for table_dataset in table_dataset_list:
        table_dataset.close()

    for table_dataset in validate_table_dataset_list:
        table_dataset.close()

    # show the parameter set with the smallest validation loss
    if 'val_loss' in train_history.history:
        best_validation_epoch = np.argmin(np.array(train_history.history["val_loss"])) + 1
        logging.info("[INFO] Best validation loss at epoch: %d" % best_validation_epoch)
    else:
        best_train_epoch = np.argmin(np.array(train_history.history["loss"])) + 1
        logging.info("[INFO] Best train loss at epoch: %d" % best_train_epoch)
    logging.info("Total time elapsed: %.2f s" % (time() - training_start_time))



# input N, output N model
def train_model_NN(args):
    platform = args.platform

    #legacy
    pileup = args.pileup
    tar_model_type = args.tar_model_type
    add_indel_length = args.add_indel_length

    exclude_training_samples = args.exclude_training_samples
    exclude_training_samples = set(exclude_training_samples.split(',')) if exclude_training_samples else set()
    add_validation_dataset = args.validation_dataset
    ochk_prefix = args.ochk_prefix if args.ochk_prefix is not None else ""
    

    if abs(args.mcv_alpha) > 1e-10:
        args.add_mcv_loss = True


    logging.info(args.__dict__)
    # import pdb; pdb.set_trace()
    

    training_start_time = time()
    import trio.param_t as param

    if args.model_cls == "Clair3_Trio_Out3":
        model = model_path.Clair3_Trio_Out3(add_mcv_loss=args.add_mcv_loss)
    elif args.model_cls == "Clair3_Trio_V_res":
        model = model_path.Clair3_Trio_V_res()
    elif args.model_cls == "Clair3_Trio_V_rres":
        model = model_path.Clair3_Trio_V_rres(is_padding=args.add_padding)
    elif args.model_cls == "Clair3_Trio_V_o1":
        model = model_path.Clair3_Trio_V_o1()
    elif args.model_cls == "Clair3_Trio_V_o2":
        model = model_path.Clair3_Trio_V_o2()
    elif args.model_cls == "Clair3_Trio_V_o1_rres":
        model = model_path.Clair3_Trio_V_o1_rres()
    else:
        raise ValueError('Unsupported model_cls name: ',args.model_cls)
    logging.info("Model class name: %s" % (args.model_cls))

    # import pdb; pdb.set_trace()


    _is_reverse_23 = args.add_reverse_23
    tensor_shape = param.ont_input_shape_trio

    if args.add_padding:
        tensor_shape = param.p_ont_input_shape_trio
    # import pdb; pdb.set_trace()

    label_size, label_shape, label_shape_cum = param.label_size_trio, param.label_shape_trio, param.label_shape_cum_trio
    
    tensor_size_one = int(tensor_shape[0] / 3)
    label_size_one = param.label_size


    # 2000, 200
    batch_size, chunk_size = param.trainBatchSize, param.chunk_size
    if args.batch_size != 0:
        batch_size = args.batch_size
    random.seed(param.RANDOM_SEED)
    np.random.seed(param.RANDOM_SEED)
    learning_rate = args.learning_rate if args.learning_rate else param.initialLearningRate
    max_epoch = args.maxEpoch if args.maxEpoch else param.maxEpoch
    task_num = len(label_shape)
    TensorShape = (tf.TensorShape([None] + tensor_shape),
         tuple(tf.TensorShape([None, label_shape[task]]) for task in range(task_num)))
    TensorDtype = (tf.int32, tuple(tf.float32 for _ in range(task_num)))

    # trio loss
    if args.add_mcv_loss:
        TensorShape = (tf.TensorShape([None] + tensor_shape),
             tuple([tf.TensorShape([None, label_shape[task]]) for task in range(task_num)] + [tf.TensorShape([None, param.label_shape[0] * 3])]))
        TensorDtype = (tf.int32, tuple([tf.float32 for _ in range(task_num+1)]))

    # import pdb; pdb.set_trace()


    bin_list = os.listdir(args.bin_fn)

    # default we exclude sample hg003 and all chr20 for training
    bin_list = [f for f in bin_list if '_20_' not in f and not exist_file_prefix(exclude_training_samples, f)]

    logging.info("[INFO] total {} training bin files: {}".format(len(bin_list), ','.join(bin_list)))
    total_data_size = 0
    table_dataset_list = []
    validate_table_dataset_list = []
    chunk_offset = np.zeros(len(bin_list), dtype=int)

    effective_label_num = None
    for bin_idx, bin_file in enumerate(bin_list):
        table_dataset = tables.open_file(os.path.join(args.bin_fn, bin_file), 'r')
        validate_table_dataset = tables.open_file(os.path.join(args.bin_fn, bin_file), 'r')
        table_dataset_list.append(table_dataset)
        validate_table_dataset_list.append(validate_table_dataset)
        chunk_num = (len(table_dataset.root.label) - batch_size) // chunk_size
        data_size = int(chunk_num * chunk_size)
        chunk_offset[bin_idx] = chunk_num
        total_data_size += data_size

    train_data_size = total_data_size * param.trainingDatasetPercentage
    validate_data_size = int((total_data_size - train_data_size) // chunk_size) * chunk_size
    train_data_size = int(train_data_size // chunk_size) * chunk_size
    train_shuffle_chunk_list, validate_shuffle_chunk_list = get_chunk_list(chunk_offset, train_data_size, chunk_size)
    # import pdb; pdb.set_trace()

    def DataGenerator(x, data_size, shuffle_chunk_list, train_flag=True):

        """
        data generator for pileup or full alignment data processing, pytables with blosc:lz4hc are used for extreme fast
        compression and decompression. random chunk shuffling and random start position to increase training model robustness.

        """

        chunk_iters = batch_size // chunk_size
        batch_num = data_size // batch_size
        position_matrix = np.empty([batch_size] + tensor_shape, np.int32)
        label = np.empty((batch_size, label_size), np.float32)

        random_start_position = np.random.randint(0, batch_size) if train_flag else 0
        if train_flag:
            np.random.shuffle(shuffle_chunk_list)

        _tar_lebel_id = 0 if tar_model_type == 'child' else 1
        # import pdb; pdb.set_trace()

        _idx = 0
        for batch_idx in range(batch_num):
            for chunk_idx in range(chunk_iters):
                offset_chunk_id = shuffle_chunk_list[batch_idx * chunk_iters + chunk_idx]
                bin_id, chunk_id = offset_chunk_id
                position_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.position_matrix[
                        random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                label[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.label[
                        random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                # for _ii, _t in enumerate(range(random_start_position + chunk_id * chunk_size, random_start_position + (chunk_id + 1) * chunk_size)):
                #     print(_ii, _t)
                #     print(x[bin_id].root.position[_t], x[bin_id].root.alt_info[_t], x[bin_id].root.label[_t])
                #     for _i in range(3):
                #         print(', '.join([str(i) for i in x[bin_id].root.label[_t].reshape(3, -1)[_i, :]]))

                # import pdb; pdb.set_trace()

            # legacy
            # _tar_label_idx = _tar_lebel_id * label_size_one
            # yield position_matrix, (
            #             label[:, _tar_label_idx + 0                 : _tar_label_idx + label_shape_cum[0]],
            #             label[:, _tar_label_idx + label_shape_cum[0]: _tar_label_idx + label_shape_cum[1]],
            #             label[:, _tar_label_idx + label_shape_cum[1]: _tar_label_idx + label_shape_cum[2]],
            #             label[:, _tar_label_idx + label_shape_cum[2]: _tar_label_idx + label_size_one]
            #         )
            # if _idx > 10:
            #     break
            # _idx += 1
            # print(len(position_matrix), len(label), label_shape_cum, label[-6])
            # for i in range(1, 600):
            #     # print(i)
            #     if position_matrix[-i].sum() != 0:
            #         print(i)
            #         print(position_matrix[-i].sum())
            #         break

            # for i in range(1, 600):
            #     # print(i)
            #     if label[-i].sum() != 0:
            #         print(i)
            #         print(label[-i].sum())
            #         break

            # import pdb; pdb.set_trace()
            _t_cum = [0] + label_shape_cum
            if args.add_mcv_loss:
                trio_pred_tar = [label[:, _t_cum[_i]: _t_cum[_i+1]] for _i in range(len(label_shape_cum)) if _i % 4 == 0]
                trio_pred_tar = np.concatenate(trio_pred_tar, axis=1)
                yield position_matrix, tuple([label[:, _t_cum[_i]: _t_cum[_i+1]] for _i in range(len(label_shape_cum))] + [trio_pred_tar])
            else:
                yield position_matrix, tuple([label[:, _t_cum[_i]: _t_cum[_i+1]] for _i in range(len(label_shape_cum))]) 

            # trio loss
            # yield position_matrix, label

            
            # trio_pred_tar = np.concatenate((label, np.multiply(label[:, :1], 10)), axis=1) 
            # # trio_pred_tar = np.concatenate((label, np.ones((len(label),1), dtype=np.float32)), axis=1) 
            # yield position_matrix, tuple([label[:, _t_cum[_i]: _t_cum[_i+1]] for _i in range(len(label_shape_cum))] + [trio_pred_tar]) 


        # two rounds for generating reversed data
        # _is_reverse_23 = False
        if _is_reverse_23:
            for batch_idx in range(batch_num):
                for chunk_idx in range(chunk_iters):
                    offset_chunk_id = shuffle_chunk_list[batch_idx * chunk_iters + chunk_idx]
                    bin_id, chunk_id = offset_chunk_id
                    position_matrix[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.position_matrix[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    label[chunk_idx * chunk_size:(chunk_idx + 1) * chunk_size] = x[bin_id].root.label[
                            random_start_position + chunk_id * chunk_size:random_start_position + (chunk_id + 1) * chunk_size]
                    # for _ii, _t in enumerate(range(random_start_position + chunk_id * chunk_size, random_start_position + (chunk_id + 1) * chunk_size)):
                    #     print(_ii, _t)
                    #     print(x[bin_id].root.position[_t], x[bin_id].root.alt_info[_t], x[bin_id].root.label[_t])
                    #     for _i in range(3):
                    #         print(', '.join([str(i) for i in x[bin_id].root.label[_t].reshape(3, -1)[_i, :]]))

                if _is_reverse_23:
                    position_matrix = np.concatenate((
                    position_matrix[:, tensor_size_one * 0 : tensor_size_one * 1,:,:], \
                    position_matrix[:, tensor_size_one * 2 : tensor_size_one * 3:,:], \
                    position_matrix[:, tensor_size_one * 1 : tensor_size_one * 2,:,:]), axis=1)

                    label = np.concatenate((
                    label[:, label_size_one * 0 : label_size_one * 1], \
                    label[:, label_size_one * 2 : label_size_one * 3], \
                    label[:, label_size_one * 1 : label_size_one * 2]), axis=1)


                # _t_cum = [0] + label_shape_cum
                # yield position_matrix, tuple([label[:, _t_cum[_i]: _t_cum[_i+1]] for _i in range(len(label_shape_cum))]) 

                _t_cum = [0] + label_shape_cum
                if args.add_mcv_loss:
                    trio_pred_tar = [label[:, _t_cum[_i]: _t_cum[_i+1]] for _i in range(len(label_shape_cum)) if _i % 4 == 0]
                    trio_pred_tar = np.concatenate(trio_pred_tar, axis=1)
                    yield position_matrix, tuple([label[:, _t_cum[_i]: _t_cum[_i+1]] for _i in range(len(label_shape_cum))] + [trio_pred_tar])
                else:
                    yield position_matrix, tuple([label[:, _t_cum[_i]: _t_cum[_i+1]] for _i in range(len(label_shape_cum))]) 


    train_dataset = tf.data.Dataset.from_generator(
        lambda: DataGenerator(table_dataset_list, train_data_size, train_shuffle_chunk_list, True), TensorDtype,
        TensorShape).prefetch(buffer_size=tf.data.experimental.AUTOTUNE)
    validate_dataset = tf.data.Dataset.from_generator(
        lambda: DataGenerator(validate_table_dataset_list, validate_data_size, validate_shuffle_chunk_list, False), TensorDtype,
        TensorShape).prefetch(buffer_size=tf.data.experimental.AUTOTUNE)

    total_steps = max_epoch * train_data_size // batch_size

    # for example in train_dataset:
    #     print(example)
    # import pdb; pdb.set_trace()

    # RectifiedAdam with warmup start
    optimizer = tfa.optimizers.Lookahead(tfa.optimizers.RectifiedAdam(
            lr=learning_rate,
            total_steps=total_steps,
            warmup_proportion=0.1,
            min_lr=learning_rate*0.75,
        ))
    logging.info("finetune optimizer with RectifiedAdam")
    

    # epochs = max_epoch
    # steps_per_epoch = train_data_size // batch_size
    # num_train_steps = steps_per_epoch * epochs
    # num_warmup_steps = int(0.1*num_train_steps)

    # init_lr = 3e-5
    # optimizer = model_path.create_BERT_optimizer(init_lr=init_lr,
    #                                           num_train_steps=num_train_steps,
    #                                           num_warmup_steps=num_warmup_steps)
    # logging.info("finetune optimizer with BERT's AdamW")
    # import pdb; pdb.set_trace()


    # import pdb; pdb.set_trace()

    # trio loss
    # task_num = 1
    loss_func = [FocalLoss(label_shape_cum, task, effective_label_num) for task in range(task_num)]
    loss_task = {"output_{}".format(task + 1): loss_func[task] for task in range(task_num)}
    metrics = {"output_{}".format(task + 1): tfa.metrics.F1Score(num_classes=label_shape[task], average='micro') for
               task in range(task_num)}

    if args.add_mcv_loss:
        if args.mcv_alpha == 0:
            args.mcv_alpha = 1
        loss_task["output_{}".format(task_num + 1)] = MCVLoss(alpha=args.mcv_alpha)

    model.compile(
        loss=loss_task,
        metrics=metrics,
        optimizer=optimizer,
        # run_eagerly=True
    )


    # model.compile(
    #     loss=MCVLoss,
    #     metrics=[test_accuracy],
    #     optimizer=optimizer,
    # )


        # run_eagerly=True
    # import pdb; pdb.set_trace()
    early_stop_callback = tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10, mode="min")
    model_save_callback = tf.keras.callbacks.ModelCheckpoint(ochk_prefix + ".{epoch:02d}", period=1, save_weights_only=False)
    model_best_callback = tf.keras.callbacks.ModelCheckpoint("best_val_loss", monitor='val_loss', save_best_only=True, mode="min")
    train_log_callback = tf.keras.callbacks.CSVLogger("training.log", separator='\t')

    # Use first 20 element to initialize tensorflow model using graph mode
    # output = model(np.array(table_dataset_list[0].root.position_matrix[:20]))

    model.build(tf.TensorShape([None] + tensor_shape))
    logging.info(model.summary(print_fn=logging.info))

    logging.info("[INFO] The size of dataset: {}".format(total_data_size))
    logging.info("[INFO] The training batch size: {}".format(batch_size))
    logging.info("[INFO] The training learning_rate: {}".format(learning_rate))
    logging.info("[INFO] Total training steps: {}".format(total_steps))
    logging.info("[INFO] Maximum training epoch: {}".format(max_epoch))
    logging.info("[INFO] MCVLoss alpha: {}".format(args.mcv_alpha))
    logging.info("[INFO] Start training...")

    validate_dataset = validate_dataset if add_validation_dataset else None
    # import pdb; pdb.set_trace()
    if args.chkpnt_fn is not None:
        model.load_weights(args.chkpnt_fn)
        logging.info("[INFO] Starting from model {}".format(args.chkpnt_fn))


    # import wandb
    # from wandb.keras import WandbCallback
    # wandb.login(key='43ff292c29b48569c32e92bbf82356581cc71cae')

    # _train_name = "trio-" + ochk_prefix.split('/')[-2]
    # wandb.init(project='clair3-trio', entity='junhao', name=_train_name,
    #        config={
    #           "learning_rate": learning_rate,
    #           "epochs": max_epoch,
    #           "batch_size": batch_size,
    #           "model_s_type": "NN",
    #           "output_tensor_shape": label_size,
    #           "tar_model_type": tar_model_type,
    #           "add_reverse_23": args.add_reverse_23,
    #        })
    # train_history = model.fit(x=train_dataset,
    #                           epochs=max_epoch,
    #                           validation_data=validate_dataset,
    #                           callbacks=[WandbCallback(save_model=False, log_weights=True), early_stop_callback, model_save_callbakck],
    #                           verbose=1,
    #                           shuffle=False)
    


    # # import pdb; pdb.set_trace()
    train_history = model.fit(x=train_dataset,
                              epochs=max_epoch,
                              validation_data=validate_dataset,
                              callbacks=[early_stop_callback, model_save_callback, model_best_callback, train_log_callback],
                              verbose=1,
                              shuffle=False)

    for table_dataset in table_dataset_list:
        table_dataset.close()

    for table_dataset in validate_table_dataset_list:
        table_dataset.close()

    # show the parameter set with the smallest validation loss
    if 'val_loss' in train_history.history:
        best_validation_epoch = np.argmin(np.array(train_history.history["val_loss"])) + 1
        logging.info("[INFO] Best validation loss at epoch: %d" % best_validation_epoch)
    else:
        best_train_epoch = np.argmin(np.array(train_history.history["loss"])) + 1
        logging.info("[INFO] Best train loss at epoch: %d" % best_train_epoch)
    logging.info("Total time elapsed: %.2f s" % (time() - training_start_time))


def main():
    parser = ArgumentParser(description="Train a Clair3 model")

    parser.add_argument('--platform', type=str, default="ont",
                        help="Sequencing platform of the input. Options: 'ont,hifi,ilmn', default: %(default)s")

    parser.add_argument('--bin_fn', type=str, default="", required=True,
                        help="Binary tensor input generated by Tensor2Bin.py, support multiple bin readers using pytables")

    parser.add_argument('--chkpnt_fn', type=str, default=None,
                        help="Input a model to resume training or for fine-tuning")

    parser.add_argument('--ochk_prefix', type=str, default=None, required=True,
                        help="Prefix for model output after each epoch")

    parser.add_argument('--tar_model_type', type=str, default='child',
                        help="target model type, child or parante, trio, default: %(default)s")

    parser.add_argument('--add_reverse_23', type=str2bool, default=False,
                        help=SUPPRESS)

    parser.add_argument('--model_arc', type=str, default="N1",
                        help="model architecture, N1 for N to 1, NN for N to N model, default: %(default)s")

    parser.add_argument('--model_cls', type=str, default="Clair3_Trio_Out3",
                        help="model class name, %(default)s")

    # options for advanced users
    parser.add_argument('--maxEpoch', type=int, default=None,
                        help="Maximum number of training epochs")

    parser.add_argument('--batch_size', type=int, default=0,
                        help="traing batch size, %(default)s")

    parser.add_argument('--learning_rate', type=float, default=1e-3,
                        help="Set the initial learning rate, default: %(default)s")

    parser.add_argument('--validation_dataset', action='store_true',
                        help="Use validation dataset when training, default: %(default)s")

    parser.add_argument('--exclude_training_samples', type=str, default=None,
                        help="Define training samples to be excluded")
    
    parser.add_argument('--add_padding', action='store_true',
                        help=SUPPRESS)


    parser.add_argument('--mcv_alpha', type=float, default=0,
                        help="Set MCVLoss rate, default: %(default)s")

    parser.add_argument('--add_mcv_loss', type=str2bool, default=False,
                        help=SUPPRESS)
    # Internal process control
    ## In pileup training mode or not
    parser.add_argument('--pileup', action='store_true', help=SUPPRESS)

    ## Add indel length for training and calling, default true for full alignment
    parser.add_argument('--add_indel_length', type=str2bool, default=False,
                        help=SUPPRESS)


    parser.add_argument('--is_debuging', type=str2bool, default=False,
                        help=SUPPRESS)



    args = parser.parse_args()

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        sys.exit(1)

    if args.model_arc == 'N1':
        train_model_N1(args)
    else:
        train_model_NN(args)



if __name__ == "__main__":
    main()
