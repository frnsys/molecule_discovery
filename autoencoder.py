"""
LSTM autoencoder for learning word and document embeddings

Adapted from: <https://github.com/ilblackdragon/tf_examples/tree/master/seq2seq>
"""

import tensorflow as tf
from tensorflow.contrib import layers

import json
# max_n = 100

def stream():
    with open('data/pubmed.dat', 'r') as f:
        for i, line in enumerate(f):
            # if i >= max_n:
            #     break
            yield json.loads(line)

# Build vocab
from collections import defaultdict
min_vocab_count = 1
docs = [d for d in stream()]
vocab = defaultdict(int)
for doc in docs:
    toks = doc['toks']
    toks = [t.lower() for t in toks['title'] + toks['abstract']]
    for tok in toks:
        vocab[tok] += 1
vocab = ['UNK', '</S>'] + [t for t, c in vocab.items() if c >= min_vocab_count]
vocab2id = {v: i for i, v in enumerate(vocab)}
with open('vocab.dat', 'w') as f:
    f.write('\n'.join(vocab))

UNK = 0
END = 1

def gen():
    for doc in docs:
        toks = doc['toks']
        toks = [t.lower() for t in toks['title'] + toks['abstract']]
        doc = [vocab2id.get(tok, UNK) for tok in toks] + [END]
        yield doc

model_dir = 'testing_save'
epochs = 10
params = {
    'vocab_size': len(vocab),
    'batch_size': 64,
    'embed_dim': 128,
    'num_units': 256,
    'optimizer': 'Adam',
    'learning_rate': 0.001
}

ds = tf.data.Dataset.from_generator(gen, output_types=tf.int32)
ds = ds.shuffle(buffer_size=params['batch_size']*10)
ds = ds.padded_batch(params['batch_size'],
                     padded_shapes=[None],
                     padding_values=END,
                     drop_remainder=True)
ds = ds.repeat(epochs)
it = ds.make_one_shot_iterator()
inputs = it.get_next()
outputs = inputs

def decoder(helper, encoder_outputs, params, scope='decoder', reuse=None):
    with tf.variable_scope(scope, reuse=reuse):
        # Prepare attention mechanism
        attention_mechanism = tf.contrib.seq2seq.BahdanauAttention(
            num_units=params['num_units'], memory=encoder_outputs)
        cell = tf.contrib.rnn.LSTMCell(num_units=params['num_units'])
        attn_cell = tf.contrib.seq2seq.AttentionWrapper(
            cell, attention_mechanism, attention_layer_size=params['num_units']/2)
        out_cell = tf.contrib.rnn.OutputProjectionWrapper(
            attn_cell, params['vocab_size'])

        # Prepare the decoder with the attention cell
        decoder = tf.contrib.seq2seq.BasicDecoder(
            cell=out_cell, helper=helper,
            initial_state=out_cell.zero_state(
                dtype=tf.float32, batch_size=params['batch_size']))
        final_outputs, final_state, final_sequence_lengths = tf.contrib.seq2seq.dynamic_decode(
            decoder=decoder, impute_finished=True)
        return final_outputs


# Embeddings layer for the input
inputs = layers.embed_sequence(
    inputs,
    vocab_size=params['vocab_size'],
    embed_dim=params['embed_dim'],
    scope='embed')

# Encoder
cell = tf.contrib.rnn.LSTMCell(num_units=params['num_units'])

# Regular LSTM
# encoder_outputs, encoder_final_state = tf.nn.dynamic_rnn(cell, inputs, dtype=tf.float32)

# Bidirectional LSTM
encoder_outputs, encoder_final_state = tf.nn.bidirectional_dynamic_rnn(cell, cell, inputs, dtype=tf.float32)
encoder_outputs = tf.concat(encoder_outputs, 2)

# Prepare output
start_tokens = tf.zeros([params['batch_size']], dtype=tf.int32)
train_output = tf.concat([tf.expand_dims(start_tokens, 1), outputs], 1)
output_lengths = tf.reduce_sum(tf.to_int32(tf.not_equal(train_output, 1)), 1)

# Use same embeddings variable as the input
output_embed = layers.embed_sequence(
    train_output,
    vocab_size=params['vocab_size'],
    embed_dim=params['embed_dim'],
    scope='embed', reuse=True)

# Setup decoder for use during training
train_helper = tf.contrib.seq2seq.TrainingHelper(output_embed, output_lengths)
train_outputs = decoder(train_helper, encoder_outputs, params)

# Prioritize examples that the model was wrong on,
# by setting weight=1 to any example where the prediction was not 1,
# i.e. incorrect
weights = tf.to_float(tf.not_equal(train_output[:, :-1], 1))

# Training loss and op
loss_op = tf.contrib.seq2seq.sequence_loss(
    train_outputs.rnn_output, outputs, weights=weights)
train_op = layers.optimize_loss(
    loss_op, global_step=tf.train.get_or_create_global_step(),
    optimizer=params['optimizer'],
    learning_rate=params['learning_rate'],
    summaries=['loss', 'learning_rate'])

# Access embeddings directly
with tf.variable_scope('embed', reuse=True):
    embeddings = tf.get_variable('embeddings')

# Setup inference/prediction decoder, reusing
# variables from the training decoder
pred_helper = tf.contrib.seq2seq.GreedyEmbeddingHelper(
    embeddings, start_tokens=start_tokens, end_token=END)
pred_outputs = decoder(pred_helper, encoder_outputs, params, reuse=True)

step_counter = tf.train.StepCounterHook(every_n_steps=100)
logger = tf.train.LoggingTensorHook({'loss': loss_op}, every_n_iter=10)
tf.logging.set_verbosity(tf.logging.INFO)
with tf.train.MonitoredTrainingSession(
    is_chief=True,
    checkpoint_dir=model_dir,
    save_checkpoint_secs=600,
    hooks=[step_counter, logger]) as sess:
    while not sess.should_stop():
        sess.run(train_op)
