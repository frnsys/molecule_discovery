"""
Model to predict ATC code for a compound
"""

import os
import json
import random
import numpy as np
import tensorflow as tf
from tqdm import trange
from rdkit import Chem
from rdkit.Chem import AllChem


def process_smile(smi):
    """Generate fingerprint from SMILES
    string for the ATC model"""
    mol = Chem.MolFromSmiles(smi)
    fpr = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return np.array([int(i) for i in fpr.ToBitString()])


class ATCModel:
    def __init__(self, n_features, n_classes, layers=[8,8], learning_rate=0.01):
        self._conf = {k: v for k, v in locals().items() if k not in ['self']}

        self.keep_prob = tf.placeholder(tf.float32)
        self.x = tf.placeholder(tf.float32, shape=(None, n_features))
        net = self._net(n_features, layers)

        # output layer
        logits = tf.layers.dense(net, n_classes, activation=None)
        probs = tf.nn.softmax(logits)
        self.predict = tf.argmax(probs, axis=1)

        self.y = tf.placeholder(tf.int64, shape=(None,))
        self.loss_op = tf.losses.sparse_softmax_cross_entropy(self.y, logits)
        self.train_op = tf.train.AdamOptimizer(learning_rate).minimize(self.loss_op)

        correct_pred = tf.equal(tf.argmax(probs, 1), self.y)
        self.acc_op = tf.reduce_mean(tf.cast(correct_pred, tf.float32))

        self.sess = tf.Session()
        self.saver = tf.train.Saver()

        tf.summary.scalar('loss', self.loss_op)
        tf.summary.scalar('accuracy', self.acc_op)
        self.summary = tf.summary.merge_all()
        self.writer = tf.summary.FileWriter('/tmp/testing', self.sess.graph)

        init = tf.global_variables_initializer()
        self.sess.run(init)

    def _net(self, n_features, layers):
        """Define the network"""
        net = tf.layers.dense(self.x, units=n_features, activation=None)
        for i, units in enumerate(layers):
            net = tf.nn.dropout(net, self.keep_prob)
            net = tf.layers.dense(net, units=units, activation=tf.nn.relu)
        return net

    def train(self, training, epochs=400, keep_prob=0.5):
        losses = []
        accuracies = []
        it = trange(epochs)
        for e in it:
            # Assemble batch
            random.shuffle(training)
            x = np.array([s[0] for s in training])
            y = np.array([s[1] for s in training])
            summary, _, err, acc = self.sess.run(
                [self.summary, self.train_op, self.loss_op, self.acc_op],
                feed_dict={
                    self.x: x,
                    self.y: y,
                    self.keep_prob: keep_prob
                }
            )
            losses.append(err)
            accuracies.append(acc)
            it.set_postfix(
                loss=np.mean(losses[-10:]) if losses else None,
                acc=np.mean(accuracies[-10:]) if losses else None)
            self.writer.add_summary(summary, e)
        return losses

    def predict(self, smis):
        X = [process_smile(smi) for smi in smis]
        return self.sess.run(self.predict, feed_dict={
            self.x: X,
            self.keep_prob: 1.0
        })

    def save(self, path):
        if not os.path.exists(path):
            os.makedirs(path)

        # Save config
        with open(os.path.join(path, 'conf.json'), 'w') as f:
            json.dump(self._conf, f)

        # Save trained params
        ckpt_path = os.path.join(path, 'model.ckpt')
        return self.saver.save(self.sess, ckpt_path)

    @classmethod
    def load(cls, path):
        # Load model with config
        with open(os.path.join(path, 'conf.json'), 'r') as f:
            conf = json.load(f)

        model = cls(**conf)

        # Load trained params
        ckpt_path = os.path.join(path, 'model.ckpt')
        model.saver.restore(model.sess, ckpt_path)
        return model