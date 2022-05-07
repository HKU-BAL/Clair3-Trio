import warnings
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', category=DeprecationWarning)
    warnings.filterwarnings("ignore", category=FutureWarning)
    from tensorflow.python.util import deprecation
    deprecation._PRINT_DEPRECATION_WARNINGS = False
    import tensorflow as tf
import logging
import numpy as np
logging.basicConfig(format='%(message)s', level=logging.INFO)
tf.get_logger().setLevel(logging.ERROR)

from clair3.task.main import GT21, GENOTYPE, VARIANT_LENGTH_1, VARIANT_LENGTH_2
import re
import trio.param_t as param
params = dict(
            float_type=tf.float32,
            task_loss_weights=[
                1,                       # gt21
                1,                       # genotype
                1,                       # variant/indel length 0
                1,                       # variant/indel length 1
                1                        # l2 loss
            ],
            output_shape=GT21.output_label_count + \
                         GENOTYPE.output_label_count + \
                         VARIANT_LENGTH_1.output_label_count + \
                         VARIANT_LENGTH_2.output_label_count,
            output_gt21_shape=GT21.output_label_count,
            output_genotype_shape=GENOTYPE.output_label_count,
            output_indel_length_shape_1=VARIANT_LENGTH_1.output_label_count,
            output_indel_length_shape_2=VARIANT_LENGTH_2.output_label_count,
            output_gt21_entropy_weights=[1] * GT21.output_label_count,
            output_genotype_entropy_weights=[1] * GENOTYPE.output_label_count,
            output_indel_length_entropy_weights_1=[1] * VARIANT_LENGTH_1.output_label_count,
            output_indel_length_entropy_weights_2=[1] * VARIANT_LENGTH_2.output_label_count,
            L3_dropout_rate=0.2,
            L4_num_units=256,
            L4_pileup_num_units=128,
            L4_dropout_rate=0.5,
            L5_1_num_units=128,
            L5_1_dropout_rate=0.2,
            L5_2_num_units=128,
            L5_2_dropout_rate=0.2,
            L5_3_num_units=128,
            L5_3_dropout_rate=0.2,
            L5_4_num_units=128,
            L5_4_dropout_rate=0.2,
            LSTM1_num_units=128,
            LSTM2_num_units=160,
            LSTM1_dropout_rate=0,
            LSTM2_dropout_rate=0.5,
            l2_regularization_lambda=param.l2RegularizationLambda,
        )

add_l2_regulation = True
L2_regularizers = tf.keras.regularizers.l2(params['l2_regularization_lambda']) if add_l2_regulation else None

class Clair3_P(tf.keras.Model):
    # Bi-lstm model for clair3 pileup input
    def __init__(self, add_indel_length=False, predict=False):
        super(Clair3_P, self).__init__()

        # output
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_pileup_num_units = params['L4_pileup_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']
        self.LSTM1_num_units = params['LSTM1_num_units']
        self.LSTM2_num_units = params['LSTM2_num_units']
        self.LSTM1_dropout_rate = params['LSTM1_dropout_rate']
        self.LSTM2_dropout_rate = params['LSTM2_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        self.add_indel_length = add_indel_length
        self.predict = predict

        self.LSTM1 = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(
            units=self.LSTM1_num_units,
            return_sequences=True,
            kernel_regularizer=L2_regularizers
        ))

        self.LSTM2 = tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(
            units=self.LSTM2_num_units,
            return_sequences=True,
            kernel_regularizer=L2_regularizers
        ))

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.L3_dropout_flatten = tf.keras.layers.Flatten()

        self.L4 = tf.keras.layers.Dense(units=self.L4_pileup_num_units, activation='selu',kernel_regularizer=L2_regularizers)

        self.L4_dropout = tf.keras.layers.Dropout(rate=self.LSTM2_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_2_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)

        self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)

        if self.add_indel_length:

            self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

            self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

            self.Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu', kernel_regularizer=L2_regularizers)

            self.Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu', kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()


    def call(self, x,):

        x = tf.cast(x, tf.float32)

        x = self.LSTM1(x)  # (batch_size, inp_seq_len, d_model)

        x = self.LSTM2(x)

        x = self.L3_dropout(x)

        x = self.L3_dropout_flatten(x)

        x = self.L4(x)

        x = self.L4_dropout(x)

        l5_1_dropout = self.L5_1_dropout(self.L5_1(x))

        l5_2_dropout = self.L5_2_dropout(self.L5_2(x))

        y_gt21_logits = self.softmax(self.Y_gt21_logits(l5_1_dropout))

        y_genotype_logits = self.softmax(self.Y_genotype_logits(l5_2_dropout))

        if self.add_indel_length:
            l5_3_dropout = self.L5_3_dropout(self.L5_3(x))

            l5_4_dropout = self.L5_4_dropout(self.L5_4(x))

            y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(l5_3_dropout))

            y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(l5_4_dropout))

            if self.predict:
                return tf.concat([y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2], axis=1)

            return [y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2]

        if self.predict:
            return tf.concat([y_gt21_logits, y_genotype_logits],axis=1)

        return [y_gt21_logits, y_genotype_logits]


class BasicConv2D(tf.keras.layers.Layer):
    def __init__(self, filters, kernel_size, strides, padding, SeparableConv=False):
        super(BasicConv2D, self).__init__()
        conv = tf.keras.layers.SeparableConv2D if SeparableConv else tf.keras.layers.Conv2D
        self.conv = conv(filters=filters,
                                           kernel_size=kernel_size,
                                           strides=strides,
                                           padding=padding,
                                           kernel_regularizer=L2_regularizers)
        self.bn = tf.keras.layers.BatchNormalization()
        self.relu = tf.keras.layers.ReLU()

    def call(self, inputs):
        output = self.conv(inputs)
        output = self.bn(output)
        output = self.relu(output)

        return output

class BasicBlock(tf.keras.layers.Layer):

    def __init__(self, filter_num, stride=1,SeparableConv=False):
        super(BasicBlock, self).__init__()
        conv = tf.keras.layers.SeparableConv2D if SeparableConv else tf.keras.layers.Conv2D

        self.conv1 = conv(filters=filter_num,
                                            kernel_size=(3, 3),
                                            strides=stride,
                                            padding="same",
                                            kernel_regularizer=L2_regularizers)
        self.bn1 = tf.keras.layers.BatchNormalization()
        self.conv2 = conv(filters=filter_num,
                                            kernel_size=(3, 3),
                                            strides=1,
                                            padding="same",
                                            kernel_regularizer=L2_regularizers)
        self.bn2 = tf.keras.layers.BatchNormalization()
        if stride != 1:
            self.downsample = tf.keras.Sequential()
            self.downsample.add(tf.keras.layers.Conv2D(filters=filter_num,
                                                       kernel_size=(1, 1),
                                                       strides=stride,
                                                       kernel_regularizer=L2_regularizers))
            self.downsample.add(tf.keras.layers.BatchNormalization())
        else:
            self.downsample = lambda x: x

    def call(self, inputs):
        residual = self.downsample(inputs)

        x = self.conv1(inputs)
        x = self.bn1(x, )
        x = tf.nn.relu(x)
        x = self.conv2(x)
        x = self.bn2(x, )

        output = tf.nn.relu(tf.keras.layers.add([residual, x]))

        return output

def make_basic_block_layer(filter_num, blocks, stride=1, SeparableConv=False):

    res_block = tf.keras.Sequential()

    res_block.add(BasicBlock(filter_num, stride=stride, SeparableConv=SeparableConv))

    for _ in range(1, blocks):
        res_block.add(BasicBlock(filter_num, stride=1,SeparableConv=SeparableConv))

    return res_block

class PyramidPolling(tf.keras.layers.Layer):
    def __init__(self, spatial_pool_size=(3, 2, 1)):
        super(PyramidPolling, self).__init__()

        self.spatial_pool_size = spatial_pool_size

        self.flatten = tf.keras.layers.Flatten()

    def call(self, x):

        height = int(x.get_shape()[1])
        width = int(x.get_shape()[2])
        # print(height, width)

        for i in range(len(self.spatial_pool_size)):

            window_h = stride_h = int(np.ceil(height / self.spatial_pool_size[i]))

            window_w = stride_w = int(np.ceil(width / self.spatial_pool_size[i]))
            # print(i, window_h, window_w)

            max_pool = tf.nn.max_pool(x, ksize=[1, window_h, window_w, 1], strides=[1, stride_h, stride_w, 1],
                                      padding='SAME')

            # print(max_pool.shape)
            if i == 0:
                pp = self.flatten(max_pool)

            else:
                pp = tf.concat([pp, self.flatten(max_pool)], axis=-1)

        return pp

class Clair3_F(tf.keras.Model):
    # Residual CNN model for clair3 full alignment input
    def __init__(self, add_indel_length=False, predict=False):
        super(Clair3_F, self).__init__()
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        self.add_indel_length = add_indel_length
        self.predict = predict

        self.conv1 = BasicConv2D(filters=64,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same",)

        self.res_block1 = make_basic_block_layer(filter_num=64,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv3 = BasicConv2D(filters=128,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block2 = make_basic_block_layer(filter_num=128,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv5 = BasicConv2D(filters=256,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block3 = make_basic_block_layer(filter_num=256,
                                            blocks=1, stride=1)

        self.pyramidpolling = PyramidPolling()

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.flatten = tf.keras.layers.Flatten()

        self.L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)

        self.L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)

        self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)

        if self.add_indel_length:
            self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

            self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

            self.Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)

            self.Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()


    def call(self, inputs):

        x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM

        x = self.conv1(x)
        x = self.res_block1(x)
        x = self.conv3(x)
        x = self.res_block2(x)
        x = self.conv5(x)
        x = self.res_block3(x)
        x = self.pyramidpolling(x)
        x = self.flatten(self.L3_dropout(x))

        x = self.L4(x)
        x = self.L4_dropout(x)

        #x = tf.keras.layers.Activation('linear', dtype='float32')(x)

        l5_1_dropout = self.L5_1_dropout(self.L5_1(x))

        l5_2_dropout = self.L5_2_dropout(self.L5_2(x))

        y_gt21_logits = self.softmax(self.Y_gt21_logits(l5_1_dropout))

        y_genotype_logits = self.softmax(self.Y_genotype_logits(l5_2_dropout))

        if self.add_indel_length:

            l5_3_dropout = self.L5_3_dropout(self.L5_3(x))

            l5_4_dropout = self.L5_4_dropout(self.L5_4(x))

            y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(l5_3_dropout))

            y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(l5_4_dropout))

            if self.predict:

                return tf.concat([y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2], axis=1)

            return [y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2]

        if self.predict:

            return tf.concat([y_gt21_logits, y_genotype_logits],axis=1)

        return [y_gt21_logits, y_genotype_logits]

    




class Clair3_Trio_Basic(tf.keras.Model):
    # Residual CNN model for clair3 full alignment input
    def __init__(self, add_indel_length=False, predict=False, is_padding=False):
        super(Clair3_Trio_Basic, self).__init__()
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        self.add_indel_length = add_indel_length
        self.is_padding = is_padding
        self.predict = predict

        self.conv1 = BasicConv2D(filters=64,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same",)

        self.res_block1 = make_basic_block_layer(filter_num=64,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv3 = BasicConv2D(filters=128,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block2 = make_basic_block_layer(filter_num=128,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv5 = BasicConv2D(filters=256,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block3 = make_basic_block_layer(filter_num=256,
                                            blocks=1, stride=1)

        self.pyramidpolling = PyramidPolling()

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.flatten = tf.keras.layers.Flatten()

        self.L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)

        self.L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)

        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)

        self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)

        if self.add_indel_length:
            self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

            self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)

            self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

            self.Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)

            self.Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()

    def call(self, inputs):

        x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM

        x = self.conv1(x)
        x = self.res_block1(x)
        x = self.conv3(x)
        x = self.res_block2(x)
        x = self.conv5(x)
        x = self.res_block3(x)
        x = self.pyramidpolling(x)
        x = self.flatten(self.L3_dropout(x))

        x = self.L4(x)
        x = self.L4_dropout(x)

        l5_1_dropout = self.L5_1_dropout(self.L5_1(x))

        l5_2_dropout = self.L5_2_dropout(self.L5_2(x))

        y_gt21_logits = self.softmax(self.Y_gt21_logits(l5_1_dropout))

        y_genotype_logits = self.softmax(self.Y_genotype_logits(l5_2_dropout))

        if self.add_indel_length:

            l5_3_dropout = self.L5_3_dropout(self.L5_3(x))

            l5_4_dropout = self.L5_4_dropout(self.L5_4(x))

            y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(l5_3_dropout))

            y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(l5_4_dropout))

            if self.predict:

                return tf.concat([y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2], axis=1)

            return [y_gt21_logits, y_genotype_logits, y_indel_length_logits_1, y_indel_length_logits_2]

        if self.predict:

            return tf.concat([y_gt21_logits, y_genotype_logits],axis=1)

        return [y_gt21_logits, y_genotype_logits]


from tensorflow.keras.layers import Layer

class MyActivityRegularizer(Layer):
  """Layer that creates an activity sparsity regularization loss."""

  def __init__(self, rate=1e-2):
    super(MyActivityRegularizer, self).__init__()
    self.rate = rate

  def call(self, inputs):
    # We use `add_loss` to create a regularization loss
    # that depends on the inputs.
    self.add_loss(self.rate * tf.reduce_sum(tf.square(inputs)))
    return inputs



class MCVLoss(tf.keras.losses.Loss):
    """
    updated version of focal loss function, for multi class classification, we remove alpha parameter, which the loss
    more stable, and add gradient clipping to avoid gradient explosion and precision overflow.
    """

    def __init__(self, gamma=0.01, alpha=1, gt_error_rate=1e-8, single_tensor_size=90):
        super(MCVLoss, self).__init__()
        self.gamma = gamma
        self.alpha = alpha
        self.gt_error_rate = gt_error_rate
        self.single_tensor_size = single_tensor_size

        self.gamma = 1
        print(self.gamma, self.alpha)
        tmp_df = np.zeros(9261)
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

        
    def call(self, y_pred_1, y_pred_2, y_pred_3):
        # # y_pred = tf.clip_by_value(y_pred, clip_value_min=1e-9, clip_value_max=1 - 1e-9)
        # y_pred_1 = y_pred[:, :self.single_tensor_size][:, :21]
        # y_pred_2 = y_pred[:, self.single_tensor_size : (self.single_tensor_size * 2)][:, :21]
        # y_pred_3 = y_pred[:, (self.single_tensor_size * 2): ][:, :21]


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



class Clair3_Trio_Out3(tf.keras.Model):
    # Residual CNN model for clair3 full alignment input
    def __init__(self, add_indel_length=False, predict=False, is_padding=False, add_mcv_loss=False):
        super(Clair3_Trio_Out3, self).__init__()
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        # legacy
        self.add_indel_length = add_indel_length
        self.is_padding = is_padding
        self.add_mcv_loss = add_mcv_loss
        self.predict = predict


        self.conv1 = BasicConv2D(filters=64,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same",)

        self.res_block1 = make_basic_block_layer(filter_num=64,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv3 = BasicConv2D(filters=128,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block2 = make_basic_block_layer(filter_num=128,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv5 = BasicConv2D(filters=256,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block3 = make_basic_block_layer(filter_num=256,
                                            blocks=1, stride=1)

        self.pyramidpolling = PyramidPolling()

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.flatten = tf.keras.layers.Flatten()

        self.L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_2_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

        self.c_Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.c_Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.c_Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.c_Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.p1_Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p1_Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p1_Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.p1_Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.p2_Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p2_Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p2_Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.p2_Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()
        # self.regularization = MyActivityRegularizer(1e-2)

    def call(self, inputs):

        # if self.is_padding:
        #     i1 = inputs[:,:,:,:-1]
        #     i2 = inputs[:,:,:,-1:]
        #     i1 = tf.cast(i1, tf.float32) / param.NORMALIZE_NUM
        #     i2 = tf.cast(i2, tf.float32)
        #     x = tf.concat([i1, i2], -1)
        # else:
        #     x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM

        x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM
    
        # import pdb; pdb.set_trace()

        x = self.conv1(x)
        x = self.res_block1(x)
        x = self.conv3(x)
        x = self.res_block2(x)
        x = self.conv5(x)
        x = self.res_block3(x)
        x = self.pyramidpolling(x)
        x = self.flatten(self.L3_dropout(x))

        x = self.L4(x)
        x = self.L4_dropout(x)
        l5_1_dropout = self.L5_1_dropout(self.L5_1(x))

        l5_2_dropout = self.L5_2_dropout(self.L5_2(x))

        l5_3_dropout = self.L5_3_dropout(self.L5_3(x))

        l5_4_dropout = self.L5_4_dropout(self.L5_4(x))


        c_y_gt21_logits = self.softmax(self.c_Y_gt21_logits(l5_1_dropout))
        c_y_genotype_logits = self.softmax(self.c_Y_genotype_logits(l5_2_dropout))
        c_y_indel_length_logits_1 = self.softmax(self.c_Y_indel_length_logits_1(l5_3_dropout))
        c_y_indel_length_logits_2 = self.softmax(self.c_Y_indel_length_logits_2(l5_4_dropout))


        p1_y_gt21_logits = self.softmax(self.p1_Y_gt21_logits(l5_1_dropout))
        p1_y_genotype_logits = self.softmax(self.p1_Y_genotype_logits(l5_2_dropout))
        p1_y_indel_length_logits_1 = self.softmax(self.p1_Y_indel_length_logits_1(l5_3_dropout))
        p1_y_indel_length_logits_2 = self.softmax(self.p1_Y_indel_length_logits_2(l5_4_dropout))

        p2_y_gt21_logits = self.softmax(self.p2_Y_gt21_logits(l5_1_dropout))
        p2_y_genotype_logits = self.softmax(self.p2_Y_genotype_logits(l5_2_dropout))
        p2_y_indel_length_logits_1 = self.softmax(self.p2_Y_indel_length_logits_1(l5_3_dropout))
        p2_y_indel_length_logits_2 = self.softmax(self.p2_Y_indel_length_logits_2(l5_4_dropout))

        if self.predict:
            return tf.concat([c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
                              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
                              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2], axis=1)

        if self.add_mcv_loss:
            return [c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
                   p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
                   p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2, \
                    tf.concat([c_y_gt21_logits, p1_y_gt21_logits, p2_y_gt21_logits], axis=1)]
        return [c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2]



class Clair3_Trio_V_o1(tf.keras.Model):
    # Residual CNN model for clair3 full alignment input
    def __init__(self, add_indel_length=False, predict=False):
        super(Clair3_Trio_V_o1, self).__init__()
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        # legacy
        self.add_indel_length = add_indel_length

        self.predict = predict

        self.conv1 = BasicConv2D(filters=64,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same",)

        self.res_block1 = make_basic_block_layer(filter_num=64,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv3 = BasicConv2D(filters=128,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block2 = make_basic_block_layer(filter_num=128,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv5 = BasicConv2D(filters=256,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block3 = make_basic_block_layer(filter_num=256,
                                            blocks=1, stride=1)

        self.pyramidpolling = PyramidPolling()

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.flatten = tf.keras.layers.Flatten()


        self.L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.c_L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.c_L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.p1_L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.p1_L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.p2_L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.p2_L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_2_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

        self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()

    def call(self, inputs):

        x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM
    
        # import pdb; pdb.set_trace()

        x = self.conv1(x)
        x = self.res_block1(x)
        x = self.conv3(x)
        x = self.res_block2(x)
        x = self.conv5(x)
        x = self.res_block3(x)
        x = self.pyramidpolling(x)
        x = self.flatten(self.L3_dropout(x))

        x = self.L4(x)
        x = self.L4_dropout(x)

        c_x = self.c_L4(x)
        c_x = self.c_L4_dropout(c_x)

        p1_x = self.p1_L4(x)
        p1_x = self.p1_L4_dropout(p1_x)

        p2_x = self.p2_L4(x)
        p2_x = self.p2_L4_dropout(p2_x)
        # import pdb; pdb.set_trace()


        c_l5_1_dropout = self.L5_1_dropout(self.L5_1(c_x))
        c_l5_2_dropout = self.L5_2_dropout(self.L5_2(c_x))
        c_l5_3_dropout = self.L5_3_dropout(self.L5_3(c_x))
        c_l5_4_dropout = self.L5_4_dropout(self.L5_4(c_x))

        c_y_gt21_logits = self.softmax(self.Y_gt21_logits(c_l5_1_dropout))
        c_y_genotype_logits = self.softmax(self.Y_genotype_logits(c_l5_2_dropout))
        c_y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(c_l5_3_dropout))
        c_y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(c_l5_4_dropout))


        p1_l5_1_dropout = self.L5_1_dropout(self.L5_1(p1_x))
        p1_l5_2_dropout = self.L5_2_dropout(self.L5_2(p1_x))
        p1_l5_3_dropout = self.L5_3_dropout(self.L5_3(p1_x))
        p1_l5_4_dropout = self.L5_4_dropout(self.L5_4(p1_x))

        p1_y_gt21_logits = self.softmax(self.Y_gt21_logits(p1_l5_1_dropout))
        p1_y_genotype_logits = self.softmax(self.Y_genotype_logits(p1_l5_2_dropout))
        p1_y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(p1_l5_3_dropout))
        p1_y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(p1_l5_4_dropout))


        p2_l5_1_dropout = self.L5_1_dropout(self.L5_1(p2_x))
        p2_l5_2_dropout = self.L5_2_dropout(self.L5_2(p2_x))
        p2_l5_3_dropout = self.L5_3_dropout(self.L5_3(p2_x))
        p2_l5_4_dropout = self.L5_4_dropout(self.L5_4(p2_x))

        p2_y_gt21_logits = self.softmax(self.Y_gt21_logits(p2_l5_1_dropout))
        p2_y_genotype_logits = self.softmax(self.Y_genotype_logits(p2_l5_2_dropout))
        p2_y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(p2_l5_3_dropout))
        p2_y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(p2_l5_4_dropout))

        if self.predict:
            return tf.concat([c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
                              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
                              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2], axis=1)

        return [c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2]


class Clair3_Trio_V_o2(tf.keras.Model):
    # Residual CNN model for clair3 full alignment input
    def __init__(self, add_indel_length=False, predict=False):
        super(Clair3_Trio_V_o2, self).__init__()
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        # legacy
        self.add_indel_length = add_indel_length

        self.predict = predict

        self.conv1 = BasicConv2D(filters=64,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same",)

        self.res_block1 = make_basic_block_layer(filter_num=64,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv3 = BasicConv2D(filters=128,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block2 = make_basic_block_layer(filter_num=128,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv5 = BasicConv2D(filters=256,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block3 = make_basic_block_layer(filter_num=256,
                                            blocks=1, stride=1)

        self.pyramidpolling = PyramidPolling()

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.flatten = tf.keras.layers.Flatten()


        # self.L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        # self.L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.c_L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.c_L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.p1_L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.p1_L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.p2_L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.p2_L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_2_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

        self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()

    def call(self, inputs):

        x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM
    
        # import pdb; pdb.set_trace()

        x = self.conv1(x)
        x = self.res_block1(x)
        x = self.conv3(x)
        x = self.res_block2(x)
        x = self.conv5(x)
        x = self.res_block3(x)
        x = self.pyramidpolling(x)
        x = self.flatten(self.L3_dropout(x))

        # x = self.L4(x)
        # x = self.L4_dropout(x)

        c_x = self.c_L4(x)
        c_x = self.c_L4_dropout(c_x)

        p1_x = self.p1_L4(x)
        p1_x = self.p1_L4_dropout(p1_x)

        p2_x = self.p2_L4(x)
        p2_x = self.p2_L4_dropout(p2_x)
        # import pdb; pdb.set_trace()


        c_l5_1_dropout = self.L5_1_dropout(self.L5_1(c_x))
        c_l5_2_dropout = self.L5_2_dropout(self.L5_2(c_x))
        c_l5_3_dropout = self.L5_3_dropout(self.L5_3(c_x))
        c_l5_4_dropout = self.L5_4_dropout(self.L5_4(c_x))

        c_y_gt21_logits = self.softmax(self.Y_gt21_logits(c_l5_1_dropout))
        c_y_genotype_logits = self.softmax(self.Y_genotype_logits(c_l5_2_dropout))
        c_y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(c_l5_3_dropout))
        c_y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(c_l5_4_dropout))


        p1_l5_1_dropout = self.L5_1_dropout(self.L5_1(p1_x))
        p1_l5_2_dropout = self.L5_2_dropout(self.L5_2(p1_x))
        p1_l5_3_dropout = self.L5_3_dropout(self.L5_3(p1_x))
        p1_l5_4_dropout = self.L5_4_dropout(self.L5_4(p1_x))

        p1_y_gt21_logits = self.softmax(self.Y_gt21_logits(p1_l5_1_dropout))
        p1_y_genotype_logits = self.softmax(self.Y_genotype_logits(p1_l5_2_dropout))
        p1_y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(p1_l5_3_dropout))
        p1_y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(p1_l5_4_dropout))


        p2_l5_1_dropout = self.L5_1_dropout(self.L5_1(p2_x))
        p2_l5_2_dropout = self.L5_2_dropout(self.L5_2(p2_x))
        p2_l5_3_dropout = self.L5_3_dropout(self.L5_3(p2_x))
        p2_l5_4_dropout = self.L5_4_dropout(self.L5_4(p2_x))

        p2_y_gt21_logits = self.softmax(self.Y_gt21_logits(p2_l5_1_dropout))
        p2_y_genotype_logits = self.softmax(self.Y_genotype_logits(p2_l5_2_dropout))
        p2_y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(p2_l5_3_dropout))
        p2_y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(p2_l5_4_dropout))

        if self.predict:
            return tf.concat([c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
                              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
                              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2], axis=1)

        # return [c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
        #       p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
        #       p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2]
        return [c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
               p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
               p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2, \
                tf.concat([c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
                              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
                              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2], axis=1)]




class Clair3_Trio_V_o1_rres(tf.keras.Model):
    # Residual CNN model for clair3 full alignment input
    def __init__(self, add_indel_length=False, predict=False):
        super(Clair3_Trio_V_o1_rres, self).__init__()
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        # legacy
        self.add_indel_length = add_indel_length

        self.predict = predict

        self.conv1 = BasicConv2D(filters=64,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same",)

        self.res_block1 = make_basic_block_layer(filter_num=64,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv3 = BasicConv2D(filters=128,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block2 = make_basic_block_layer(filter_num=128,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv5 = BasicConv2D(filters=256,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block3 = make_basic_block_layer(filter_num=256,
                                            blocks=1, stride=1)

        self.conv6 = BasicConv2D(filters=512,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block4 = make_basic_block_layer(filter_num=512,
                                            blocks=1, stride=1)


        self.pyramidpolling = PyramidPolling()

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.flatten = tf.keras.layers.Flatten()


        # self.L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        # self.L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.c_L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.c_L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.p1_L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.p1_L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.p2_L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.p2_L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_2_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

        self.Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()

    def call(self, inputs):
        x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM
    
        # import pdb; pdb.set_trace()

        _input_sz = x.shape[1] // 3
        x1 = x[:,_input_sz * 0 :_input_sz * 1,:,:]
        x2 = x[:,_input_sz * 1 :_input_sz * 2,:,:]
        x3 = x[:,_input_sz * 2 :_input_sz * 3,:,:]

        x1 = self.conv1(x1)
        x1 = self.res_block1(x1)
        x1 = self.conv3(x1)
        x1 = self.res_block2(x1)
        x1 = self.conv5(x1)
        x1 = self.res_block3(x1)
        # x1 = self.pyramidpolling(x1)
        # x1 = self.flatten(self.L3_dropout(x1))

        x2 = self.conv1(x2)
        x2 = self.res_block1(x2)
        x2 = self.conv3(x2)
        x2 = self.res_block2(x2)
        x2 = self.conv5(x2)
        x2 = self.res_block3(x2)
        # x2 = self.pyramidpolling(x2)
        # x2 = self.flatten(self.L3_dropout(x2))

        x3 = self.conv1(x3)
        x3 = self.res_block1(x3)
        x3 = self.conv3(x3)
        x3 = self.res_block2(x3)
        x3 = self.conv5(x3)
        x3 = self.res_block3(x3)
        # x3 = self.pyramidpolling(x3)
        # x3 = self.flatten(self.L3_dropout(x3))

        x = tf.concat([x1, x2, x3], 1)

        # import pdb; pdb.set_trace()

        x = self.conv6(x)
        x = self.res_block4(x)
        x = self.pyramidpolling(x)
        x = self.flatten(self.L3_dropout(x))

        # x = self.L4(x)
        # x = self.L4_dropout(x)

        c_x = self.c_L4(x)
        c_x = self.c_L4_dropout(c_x)

        p1_x = self.p1_L4(x)
        p1_x = self.p1_L4_dropout(p1_x)

        p2_x = self.p2_L4(x)
        p2_x = self.p2_L4_dropout(p2_x)
        # import pdb; pdb.set_trace()


        c_l5_1_dropout = self.L5_1_dropout(self.L5_1(c_x))
        c_l5_2_dropout = self.L5_2_dropout(self.L5_2(c_x))
        c_l5_3_dropout = self.L5_3_dropout(self.L5_3(c_x))
        c_l5_4_dropout = self.L5_4_dropout(self.L5_4(c_x))

        c_y_gt21_logits = self.softmax(self.Y_gt21_logits(c_l5_1_dropout))
        c_y_genotype_logits = self.softmax(self.Y_genotype_logits(c_l5_2_dropout))
        c_y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(c_l5_3_dropout))
        c_y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(c_l5_4_dropout))


        p1_l5_1_dropout = self.L5_1_dropout(self.L5_1(p1_x))
        p1_l5_2_dropout = self.L5_2_dropout(self.L5_2(p1_x))
        p1_l5_3_dropout = self.L5_3_dropout(self.L5_3(p1_x))
        p1_l5_4_dropout = self.L5_4_dropout(self.L5_4(p1_x))

        p1_y_gt21_logits = self.softmax(self.Y_gt21_logits(p1_l5_1_dropout))
        p1_y_genotype_logits = self.softmax(self.Y_genotype_logits(p1_l5_2_dropout))
        p1_y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(p1_l5_3_dropout))
        p1_y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(p1_l5_4_dropout))


        p2_l5_1_dropout = self.L5_1_dropout(self.L5_1(p2_x))
        p2_l5_2_dropout = self.L5_2_dropout(self.L5_2(p2_x))
        p2_l5_3_dropout = self.L5_3_dropout(self.L5_3(p2_x))
        p2_l5_4_dropout = self.L5_4_dropout(self.L5_4(p2_x))

        p2_y_gt21_logits = self.softmax(self.Y_gt21_logits(p2_l5_1_dropout))
        p2_y_genotype_logits = self.softmax(self.Y_genotype_logits(p2_l5_2_dropout))
        p2_y_indel_length_logits_1 = self.softmax(self.Y_indel_length_logits_1(p2_l5_3_dropout))
        p2_y_indel_length_logits_2 = self.softmax(self.Y_indel_length_logits_2(p2_l5_4_dropout))

        if self.predict:
            return tf.concat([c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
                              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
                              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2], axis=1)

        return [c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2]



# Clair3 Trio, using shared res_block
class Clair3_Trio_V_res(tf.keras.Model):
    def __init__(self, add_indel_length=False, predict=False):
        super(Clair3_Trio_V_res, self).__init__()
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        # legacy
        self.add_indel_length = add_indel_length

        self.predict = predict

        self.conv1 = BasicConv2D(filters=64,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same",)

        self.res_block1 = make_basic_block_layer(filter_num=64,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv3 = BasicConv2D(filters=128,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block2 = make_basic_block_layer(filter_num=128,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv5 = BasicConv2D(filters=256,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block3 = make_basic_block_layer(filter_num=256,
                                            blocks=1, stride=1)

        self.pyramidpolling = PyramidPolling()

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.flatten = tf.keras.layers.Flatten()

        self.L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_2_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

        self.c_Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.c_Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.c_Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.c_Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.p1_Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p1_Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p1_Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.p1_Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.p2_Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p2_Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p2_Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.p2_Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()

    def call(self, inputs):

        x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM
    
        # import pdb; pdb.set_trace()

        _input_sz = x.shape[1] // 3
        x1 = x[:,_input_sz * 0 :_input_sz * 1,:,:]
        x2 = x[:,_input_sz * 1 :_input_sz * 2,:,:]
        x3 = x[:,_input_sz * 2 :_input_sz * 3,:,:]

        x1 = self.conv1(x1)
        x1 = self.res_block1(x1)
        x1 = self.conv3(x1)
        x1 = self.res_block2(x1)
        x1 = self.conv5(x1)
        x1 = self.res_block3(x1)
        x1 = self.pyramidpolling(x1)
        x1 = self.flatten(self.L3_dropout(x1))

        x2 = self.conv1(x2)
        x2 = self.res_block1(x2)
        x2 = self.conv3(x2)
        x2 = self.res_block2(x2)
        x2 = self.conv5(x2)
        x2 = self.res_block3(x2)
        x2 = self.pyramidpolling(x2)
        x2 = self.flatten(self.L3_dropout(x2))

        x3 = self.conv1(x3)
        x3 = self.res_block1(x3)
        x3 = self.conv3(x3)
        x3 = self.res_block2(x3)
        x3 = self.conv5(x3)
        x3 = self.res_block3(x3)
        x3 = self.pyramidpolling(x3)
        x3 = self.flatten(self.L3_dropout(x3))

        # import pdb; pdb.set_trace()
        x = tf.concat([x1, x2, x3], 1)

        x = self.L4(x)
        x = self.L4_dropout(x)

        l5_1_dropout = self.L5_1_dropout(self.L5_1(x))

        l5_2_dropout = self.L5_2_dropout(self.L5_2(x))

        l5_3_dropout = self.L5_3_dropout(self.L5_3(x))

        l5_4_dropout = self.L5_4_dropout(self.L5_4(x))


        c_y_gt21_logits = self.softmax(self.c_Y_gt21_logits(l5_1_dropout))
        c_y_genotype_logits = self.softmax(self.c_Y_genotype_logits(l5_2_dropout))
        c_y_indel_length_logits_1 = self.softmax(self.c_Y_indel_length_logits_1(l5_3_dropout))
        c_y_indel_length_logits_2 = self.softmax(self.c_Y_indel_length_logits_2(l5_4_dropout))


        p1_y_gt21_logits = self.softmax(self.p1_Y_gt21_logits(l5_1_dropout))
        p1_y_genotype_logits = self.softmax(self.p1_Y_genotype_logits(l5_2_dropout))
        p1_y_indel_length_logits_1 = self.softmax(self.p1_Y_indel_length_logits_1(l5_3_dropout))
        p1_y_indel_length_logits_2 = self.softmax(self.p1_Y_indel_length_logits_2(l5_4_dropout))

        p2_y_gt21_logits = self.softmax(self.p2_Y_gt21_logits(l5_1_dropout))
        p2_y_genotype_logits = self.softmax(self.p2_Y_genotype_logits(l5_2_dropout))
        p2_y_indel_length_logits_1 = self.softmax(self.p2_Y_indel_length_logits_1(l5_3_dropout))
        p2_y_indel_length_logits_2 = self.softmax(self.p2_Y_indel_length_logits_2(l5_4_dropout))

        if self.predict:
            return tf.concat([c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
                              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
                              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2], axis=1)

        return [c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2]


# Clair3 Trio, using shared res_block
class Clair3_Trio_V_rres(tf.keras.Model):
    def __init__(self, add_indel_length=False, predict=False, is_padding=False):
        super(Clair3_Trio_V_rres, self).__init__()
        self.output_gt21_shape = params['output_gt21_shape']
        self.output_genotype_shape = params['output_genotype_shape']
        self.output_indel_length_shape_1 = params['output_indel_length_shape_1']
        self.output_indel_length_shape_2 = params['output_indel_length_shape_2']

        self.L3_dropout_rate = params['L3_dropout_rate']
        self.L4_num_units = params['L4_num_units']
        self.L4_dropout_rate = params['L4_dropout_rate']
        self.L5_1_num_units = params['L5_1_num_units']
        self.L5_1_dropout_rate = params['L5_1_dropout_rate']
        self.L5_2_num_units = params['L5_2_num_units']
        self.L5_2_dropout_rate = params['L5_2_dropout_rate']
        self.L5_3_num_units = params['L5_3_num_units']
        self.L5_3_dropout_rate = params['L5_3_dropout_rate']
        self.L5_4_num_units = params['L5_4_num_units']
        self.L5_4_dropout_rate = params['L5_4_dropout_rate']

        self.output_label_split = [
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2,
            self.output_gt21_shape,
            self.output_genotype_shape,
            self.output_indel_length_shape_1,
            self.output_indel_length_shape_2
        ]

        # legacy
        self.add_indel_length = add_indel_length
        self.is_padding = is_padding

        self.predict = predict

        self.conv1 = BasicConv2D(filters=64,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same",)

        self.res_block1 = make_basic_block_layer(filter_num=64,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv3 = BasicConv2D(filters=128,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block2 = make_basic_block_layer(filter_num=128,
                                            blocks=1, stride=1, SeparableConv=False)

        self.conv5 = BasicConv2D(filters=256,
                                 kernel_size=(3, 3),
                                 strides=2,
                                 padding="same")

        self.res_block3 = make_basic_block_layer(filter_num=256,
                                            blocks=1, stride=1)


        self.pyramidpolling = PyramidPolling()

        self.L3_dropout = tf.keras.layers.Dropout(rate=self.L3_dropout_rate)

        self.flatten = tf.keras.layers.Flatten()

        self.L4 = tf.keras.layers.Dense(units=self.L4_num_units, activation='selu',kernel_regularizer=L2_regularizers)
        self.L4_dropout = tf.keras.layers.Dropout(rate=self.L4_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_1 = tf.keras.layers.Dense(units=self.L5_1_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_1_dropout = tf.keras.layers.Dropout(rate=self.L5_1_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_2 = tf.keras.layers.Dense(units=self.L5_2_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_2_dropout = tf.keras.layers.Dropout(rate=self.L5_2_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_3 = tf.keras.layers.Dense(units=self.L5_3_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_3_dropout = tf.keras.layers.Dropout(rate=self.L5_3_dropout_rate, seed=param.OPERATION_SEED)

        self.L5_4 = tf.keras.layers.Dense(units=self.L5_4_num_units, activation='selu', kernel_regularizer=L2_regularizers)
        self.L5_4_dropout = tf.keras.layers.Dropout(rate=self.L5_4_dropout_rate, seed=param.OPERATION_SEED)

        self.c_Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.c_Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.c_Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.c_Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.p1_Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p1_Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p1_Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.p1_Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.p2_Y_gt21_logits = tf.keras.layers.Dense(units=self.output_gt21_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p2_Y_genotype_logits = tf.keras.layers.Dense(units=self.output_genotype_shape, activation='selu', kernel_regularizer=L2_regularizers)
        self.p2_Y_indel_length_logits_1 = tf.keras.layers.Dense(units=self.output_indel_length_shape_1, activation='selu',kernel_regularizer=L2_regularizers)
        self.p2_Y_indel_length_logits_2 = tf.keras.layers.Dense(units=self.output_indel_length_shape_2, activation='selu',kernel_regularizer=L2_regularizers)

        self.softmax = tf.keras.layers.Softmax()

    def call(self, inputs):

        if self.is_padding:
            i1 = inputs[:,:,:,:-1]
            i2 = inputs[:,:,:,-1:]
            i1 = tf.cast(i1, tf.float32) / param.NORMALIZE_NUM
            i2 = tf.cast(i2, tf.float32)
            x = tf.concat([i1, i2], -1)
        else:
            x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM

        # x = tf.cast(inputs, tf.float32) / param.NORMALIZE_NUM
    

        _input_sz = x.shape[1] // 3
        x1 = x[:,_input_sz * 0 :_input_sz * 1,:,:]
        x2 = x[:,_input_sz * 1 :_input_sz * 2,:,:]
        x3 = x[:,_input_sz * 2 :_input_sz * 3,:,:]
        # import pdb; pdb.set_trace()

        x1 = self.conv1(x1)
        x1 = self.res_block1(x1)
        x1 = self.conv3(x1)
        x1 = self.res_block2(x1)
        # x1 = self.conv5(x1)
        # x1 = self.res_block3(x1)
        # x1 = self.pyramidpolling(x1)
        # x1 = self.flatten(self.L3_dropout(x1))

        x2 = self.conv1(x2)
        x2 = self.res_block1(x2)
        x2 = self.conv3(x2)
        x2 = self.res_block2(x2)
        # x2 = self.conv5(x2)
        # x2 = self.res_block3(x2)
        # x2 = self.pyramidpolling(x2)
        # x2 = self.flatten(self.L3_dropout(x2))

        x3 = self.conv1(x3)
        x3 = self.res_block1(x3)
        x3 = self.conv3(x3)
        x3 = self.res_block2(x3)
        # x3 = self.conv5(x3)
        # x3 = self.res_block3(x3)
        # x3 = self.pyramidpolling(x3)
        # x3 = self.flatten(self.L3_dropout(x3))

        x = tf.concat([x1, x2, x3], 1)

        # import pdb; pdb.set_trace()

        x = self.conv5(x)
        x = self.res_block3(x)
        x = self.pyramidpolling(x)
        x = self.flatten(self.L3_dropout(x))

        # import pdb; pdb.set_trace()

        x = self.L4(x)
        x = self.L4_dropout(x)

        l5_1_dropout = self.L5_1_dropout(self.L5_1(x))

        l5_2_dropout = self.L5_2_dropout(self.L5_2(x))

        l5_3_dropout = self.L5_3_dropout(self.L5_3(x))

        l5_4_dropout = self.L5_4_dropout(self.L5_4(x))


        c_y_gt21_logits = self.softmax(self.c_Y_gt21_logits(l5_1_dropout))
        c_y_genotype_logits = self.softmax(self.c_Y_genotype_logits(l5_2_dropout))
        c_y_indel_length_logits_1 = self.softmax(self.c_Y_indel_length_logits_1(l5_3_dropout))
        c_y_indel_length_logits_2 = self.softmax(self.c_Y_indel_length_logits_2(l5_4_dropout))


        p1_y_gt21_logits = self.softmax(self.p1_Y_gt21_logits(l5_1_dropout))
        p1_y_genotype_logits = self.softmax(self.p1_Y_genotype_logits(l5_2_dropout))
        p1_y_indel_length_logits_1 = self.softmax(self.p1_Y_indel_length_logits_1(l5_3_dropout))
        p1_y_indel_length_logits_2 = self.softmax(self.p1_Y_indel_length_logits_2(l5_4_dropout))

        p2_y_gt21_logits = self.softmax(self.p2_Y_gt21_logits(l5_1_dropout))
        p2_y_genotype_logits = self.softmax(self.p2_Y_genotype_logits(l5_2_dropout))
        p2_y_indel_length_logits_1 = self.softmax(self.p2_Y_indel_length_logits_1(l5_3_dropout))
        p2_y_indel_length_logits_2 = self.softmax(self.p2_Y_indel_length_logits_2(l5_4_dropout))

        if self.predict:
            return tf.concat([c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
                              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
                              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2], axis=1)

        return [c_y_gt21_logits, c_y_genotype_logits, c_y_indel_length_logits_1, c_y_indel_length_logits_2, \
              p1_y_gt21_logits, p1_y_genotype_logits, p1_y_indel_length_logits_1, p1_y_indel_length_logits_2, \
              p2_y_gt21_logits, p2_y_genotype_logits, p2_y_indel_length_logits_1, p2_y_indel_length_logits_2]



class WarmUp(tf.keras.optimizers.schedules.LearningRateSchedule):
  """Applies a warmup schedule on a given learning rate decay schedule."""

  def __init__(self,
               initial_learning_rate,
               decay_schedule_fn,
               warmup_steps,
               power=1.0,
               name=None):
    super(WarmUp, self).__init__()
    self.initial_learning_rate = initial_learning_rate
    self.warmup_steps = warmup_steps
    self.power = power
    self.decay_schedule_fn = decay_schedule_fn
    self.name = name

  def __call__(self, step):
    with tf.name_scope(self.name or 'WarmUp') as name:
      # Implements polynomial warmup. i.e., if global_step < warmup_steps, the
      # learning rate will be `global_step/num_warmup_steps * init_lr`.
      global_step_float = tf.cast(step, tf.float32)
      warmup_steps_float = tf.cast(self.warmup_steps, tf.float32)
      warmup_percent_done = global_step_float / warmup_steps_float
      warmup_learning_rate = (
          self.initial_learning_rate *
          tf.math.pow(warmup_percent_done, self.power))
      return tf.cond(
          global_step_float < warmup_steps_float,
          lambda: warmup_learning_rate,
          lambda: self.decay_schedule_fn(step),
          name=name)


class AdamWeightDecay(tf.keras.optimizers.Adam):
  """Adam enables L2 weight decay and clip_by_global_norm on gradients.
  Just adding the square of the weights to the loss function is *not* the
  correct way of using L2 regularization/weight decay with Adam, since that will
  interact with the m and v parameters in strange ways.
  Instead we want to decay the weights in a manner that doesn't interact with
  the m/v parameters. This is equivalent to adding the square of the weights to
  the loss with plain (non-momentum) SGD.
  """

  def __init__(self,
               learning_rate=0.001,
               beta_1=0.9,
               beta_2=0.999,
               epsilon=1e-7,
               amsgrad=False,
               weight_decay_rate=0.0,
               include_in_weight_decay=None,
               exclude_from_weight_decay=None,
               gradient_clip_norm=1.0,
               name='AdamWeightDecay',
               **kwargs):
    super(AdamWeightDecay, self).__init__(learning_rate, beta_1, beta_2,
                                          epsilon, amsgrad, name, **kwargs)
    self.weight_decay_rate = weight_decay_rate
    self.gradient_clip_norm = gradient_clip_norm
    self._include_in_weight_decay = include_in_weight_decay
    self._exclude_from_weight_decay = exclude_from_weight_decay
    logging.info('gradient_clip_norm=%f', gradient_clip_norm)

  @classmethod
  def from_config(cls, config):
    """Creates an optimizer from its config with WarmUp custom object."""
    custom_objects = {'WarmUp': WarmUp}
    return super(AdamWeightDecay, cls).from_config(
        config, custom_objects=custom_objects)

  def _prepare_local(self, var_device, var_dtype, apply_state):
    super(AdamWeightDecay, self)._prepare_local(var_device, var_dtype,
                                                apply_state)
    apply_state[(var_device, var_dtype)]['weight_decay_rate'] = tf.constant(
        self.weight_decay_rate, name='adam_weight_decay_rate')

  def _decay_weights_op(self, var, learning_rate, apply_state):
    do_decay = self._do_use_weight_decay(var.name)
    if do_decay:
      return var.assign_sub(
          learning_rate * var *
          apply_state[(var.device, var.dtype.base_dtype)]['weight_decay_rate'],
          use_locking=self._use_locking)
    return tf.no_op()

  def apply_gradients(self,
                      grads_and_vars,
                      name=None,
                      experimental_aggregate_gradients=True):
    grads, tvars = list(zip(*grads_and_vars))
    if experimental_aggregate_gradients and self.gradient_clip_norm > 0.0:
      # when experimental_aggregate_gradients = False, apply_gradients() no
      # longer implicitly allreduce gradients, users manually allreduce gradient
      # and passed the allreduced grads_and_vars. For now, the
      # clip_by_global_norm will be moved to before the explicit allreduce to
      # keep the math the same as TF 1 and pre TF 2.2 implementation.
      (grads, _) = tf.clip_by_global_norm(
          grads, clip_norm=self.gradient_clip_norm)
    return super(AdamWeightDecay, self).apply_gradients(
        zip(grads, tvars),
        name=name,
        experimental_aggregate_gradients=experimental_aggregate_gradients)

  def _get_lr(self, var_device, var_dtype, apply_state):
    """Retrieves the learning rate with the given state."""
    if apply_state is None:
      return self._decayed_lr_t[var_dtype], {}

    apply_state = apply_state or {}
    coefficients = apply_state.get((var_device, var_dtype))
    if coefficients is None:
      coefficients = self._fallback_apply_state(var_device, var_dtype)
      apply_state[(var_device, var_dtype)] = coefficients

    return coefficients['lr_t'], dict(apply_state=apply_state)

  def _resource_apply_dense(self, grad, var, apply_state=None):
    lr_t, kwargs = self._get_lr(var.device, var.dtype.base_dtype, apply_state)
    decay = self._decay_weights_op(var, lr_t, apply_state)
    with tf.control_dependencies([decay]):
      return super(AdamWeightDecay,
                   self)._resource_apply_dense(grad, var, **kwargs)

  def _resource_apply_sparse(self, grad, var, indices, apply_state=None):
    lr_t, kwargs = self._get_lr(var.device, var.dtype.base_dtype, apply_state)
    decay = self._decay_weights_op(var, lr_t, apply_state)
    with tf.control_dependencies([decay]):
      return super(AdamWeightDecay,
                   self)._resource_apply_sparse(grad, var, indices, **kwargs)

  def get_config(self):
    config = super(AdamWeightDecay, self).get_config()
    config.update({
        'weight_decay_rate': self.weight_decay_rate,
    })
    return config

  def _do_use_weight_decay(self, param_name):
    """Whether to use L2 weight decay for `param_name`."""
    if self.weight_decay_rate == 0:
      return False

    if self._include_in_weight_decay:
      for r in self._include_in_weight_decay:
        if re.search(r, param_name) is not None:
          return True

    if self._exclude_from_weight_decay:
      for r in self._exclude_from_weight_decay:
        if re.search(r, param_name) is not None:
          return False
    return True

def create_BERT_optimizer(init_lr,
                     num_train_steps,
                     num_warmup_steps,
                     end_lr=0.0,
                     beta_1=0.9):
    lr_schedule = tf.keras.optimizers.schedules.PolynomialDecay(
      initial_learning_rate=init_lr,
      decay_steps=num_train_steps,
      end_learning_rate=end_lr)
    if num_warmup_steps:
        lr_schedule = WarmUp(
            initial_learning_rate=init_lr,
            decay_schedule_fn=lr_schedule,
            warmup_steps=num_warmup_steps)

    logging.info('using Adamw optimizer')
    optimizer = AdamWeightDecay(
        learning_rate=lr_schedule,
        weight_decay_rate=0.01,
        beta_1=beta_1,
        beta_2=0.999,
        epsilon=1e-6,
        exclude_from_weight_decay=['LayerNorm', 'layer_norm', 'bias'])

    return optimizer

