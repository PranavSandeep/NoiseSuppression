from tensorflow.keras.layers import (Input, Conv2D, SeparableConv2D,
                                     LeakyReLU, BatchNormalization,
                                     Concatenate, GRU, Dense, Reshape)
from tensorflow.keras.models import Model
import tensorflow as tf


class ReduceMean(tf.keras.layers.Layer):
    def __init__(self, axis, **kwargs):
        super().__init__(**kwargs)
        self.axis = axis

    def call(self, x):
        return tf.reduce_mean(x, axis=self.axis)

    def compute_output_shape(self, input_shape):
        shape = list(input_shape)
        del shape[self.axis]
        return tuple(shape)

    def get_config(self):
        config = super().get_config()
        config.update({"axis": self.axis})
        return config


#Went from Conv1D to Conv2D, didn't commit it tho.

def red_conv_block(x, filters):
    x = SeparableConv2D(filters, (3, 1), padding="same")(x)
    x = BatchNormalization()(x)
    x = LeakyReLU()(x)
    return x

def red_deconv_block(x, filters):
    x = Conv2D(filters, (3, 1), padding="same")(x)
    x = BatchNormalization()(x)
    x = LeakyReLU()(x)
    return x




def build_rced_model(input_shape=(257, None, 1), gru_units=256, batch_size=1, stateful=True):
    F, T = input_shape[0], input_shape[1]

    if stateful:
        inputs = Input(batch_shape=(batch_size, F, T))
    else:
        inputs = Input(shape=(F, T))


    #   (B, F=257, T=5) → (B, F=257, T=5, 1)
    x = Reshape((F, T, 1))(inputs)

    encoding_filters = [32, 64, 128]
    decoding_filters = [128, 64, 32]
    skips = []

    # ── Encoder: Conv2D(3,1) slides over F=257, T=5 is untouched ─────────────
    for f in encoding_filters:
        x = red_conv_block(x, f)
        skips.append(x)

    # Bottleneck
    x = red_conv_block(x, 256)


    x_gru = ReduceMean(axis=1)(x)        # (B, T=5, 256)  ← F collapsed

    x_gru = GRU(gru_units,
                return_sequences=True,
                stateful=stateful,
                name="gru")(x_gru)       # (B, T=5, gru_units)

    x_gru = Dense(256)(x_gru)            # (B, T=5, 256)
    x_gru = Reshape((1, T, 256))(x_gru)  # (B, 1, T=5, 256) ← broadcast over F

    # Inject temporal context back into every freq bin
    x = x + x_gru                        # (B, F=257, T=5, 256)
    x = Conv2D(256, (1, 1))(x)           # pointwise mix after addition

    # ── Decoder ───────────────────────────────────────────────────────────────
    for f, skip in zip(decoding_filters, reversed(skips)):
        x = red_deconv_block(x, f)        # (B, F=257, T=5, f)
        x = Concatenate()([x, skip])      # (B, F=257, T=5, f + skip_f)

    # ── Mask output ───────────────────────────────────────────────────────────
    mask = Conv2D(1, (1, 1), activation="sigmoid")(x)  # (B, F=257, T=5, 1)
    mask = Reshape((F, T))(mask)                        # (B, F=257, T=5)

    inputs_4d = Reshape((F, T, 1))(inputs)
    mask_4d   = Reshape((F, T, 1))(mask)
    output = mask_4d * inputs_4d                        # (B, F=257, T=5, 1)

    model = Model(inputs, output)
    return model


def asymmetric_magnitude_loss(alpha=3.0, beta=1.0, eps=1e-8):
    def loss_fn(y_true, y_pred):
        residual = y_pred - y_true
        under_suppressed = tf.maximum(residual, 0.0)
        over_suppressed  = tf.maximum(-residual, 0.0)
        loss = alpha * tf.square(under_suppressed) + beta * tf.square(over_suppressed)
        return tf.reduce_mean(loss)
    return loss_fn


model = build_rced_model(stateful=False)
#These values work well
model.compile(optimizer="adam", loss=asymmetric_magnitude_loss(alpha=4.0, beta=1.0))
model.summary()