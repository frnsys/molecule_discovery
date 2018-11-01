import random
import numpy as np
import tensorflow as tf
from tqdm import trange

class NNet:
    def __init__(self, n_features, n_classes, layers=[8,8], learning_rate=0.01):
        self.keep_prob = tf.placeholder(tf.float32)
        self.x = self._input(n_features)
        net = self._net(n_features, layers)


        # output layer
        logits = tf.layers.dense(net, n_classes, activation=None)
        probs = tf.nn.sigmoid(logits)
        self.threshold = tf.placeholder(tf.float32, shape=())
        self.predict = tf.cast(tf.greater(probs, self.threshold), tf.int64)

        self.y = tf.placeholder(tf.float32, shape=(None, n_classes))

        self.loss_op = tf.losses.sigmoid_cross_entropy(self.y, logits)
        self.train_op = tf.train.AdamOptimizer(learning_rate).minimize(self.loss_op)

        init = tf.group(tf.global_variables_initializer(), tf.local_variables_initializer())
        self.sess = tf.Session()
        self.sess.run(init)

    def _input(self, n_features):
        """Define input placeholder for the network"""
        return tf.placeholder(tf.float32, shape=(None, n_features))

    def _net(self, n_features, layers):
        """Define the network"""
        net = tf.layers.dense(self.x, units=n_features, activation=None)
        for i, units in enumerate(layers):
            net = tf.nn.dropout(net, self.keep_prob)
            net = tf.layers.dense(net, units=units, activation=tf.nn.relu)
        return net

    def train(self, training, epochs=400, keep_prob=0.5):
        losses = []
        it = trange(epochs)
        for e in it:
            # Assemble batch
            random.shuffle(training)
            x = np.array([s[0] for s in training])
            y = np.array([s[1] for s in training])
            _, err = self.sess.run(
                [self.train_op, self.loss_op],
                feed_dict={
                    self.x: x,
                    self.y: y,
                    self.keep_prob: keep_prob
                }
            )
            losses.append(err)
            it.set_postfix(
                loss=err,
                mean_loss=np.mean(losses[-10:]) if losses else None)
        return losses

    def run(self, X, threshold=0.3):
        return self.sess.run(self.predict, feed_dict={
            self.x: X,
            self.threshold: threshold,
            self.keep_prob: 1.0
        })
