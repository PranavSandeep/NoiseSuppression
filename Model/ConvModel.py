from tensorflow.keras.layers import Input, Conv2D, Conv2DTranspose, LeakyReLU, BatchNormalization, Concatenate
from tensorflow.keras.models import Model
import tensorflow.keras.backend as K
import tensorflow as tf

def red_conv_block(x, filters):
    x = Conv2D(filters, 3, padding="same")(x)
    x = BatchNormalization()(x)
    x = LeakyReLU()(x)
    return x

def red_deconv_block(x, filters):
    x = Conv2DTranspose(filters, 3, strides=1, padding="same")(x)
    x = BatchNormalization()(x)
    x = LeakyReLU()(x)
    return x

def build_rced_model(input_shape=(257, 251, 1)):
    inputs = Input(shape=input_shape)
    x = inputs

    encoding_filters = [32, 64, 128]
    decoding_filters = [128, 64, 32]  # reverse for decoder

    skips = []

    # Encoder
    for f in encoding_filters:
        x = red_conv_block(x, f)
        skips.append(x)

    # Bottleneck
    x = red_conv_block(x, 256)

    # Decoder
    for f, skip in zip(decoding_filters, reversed(skips)):
        x = red_deconv_block(x, f)
        x = Concatenate()([x, skip])

    # Output
    output = Conv2D(1, kernel_size=1, activation="linear")(x)
    return Model(inputs, output)

def sdr_loss(y_true, y_pred, epsilon=1e-8):
    signal_power = K.sum(K.square(y_true), axis=[1, 2, 3]) + epsilon
    noise_power = K.sum(K.square(y_pred - y_true), axis=[1, 2, 3]) + epsilon
    sdr = 10.0 * tf.math.log(signal_power / noise_power) / tf.math.log(tf.constant(10.0))
    return -sdr

model = build_rced_model()
model.compile(optimizer="adam", loss=sdr_loss)
model.summary()

# Learning rate scheduler
from tensorflow.keras.callbacks import ReduceLROnPlateau
reduce_lr = ReduceLROnPlateau(monitor='loss', factor=0.5, patience=2, min_lr=1e-6, verbose=1)

