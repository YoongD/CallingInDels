from tensorflow.contrib import layers
from tensorflow.contrib import rnn
from sklearn.utils import shuffle
from ast import literal_eval

import tensorflow as tf
import numpy as np
import math
import sys
import csv

#tf.set_random_seed(777)  # for reproducibility
learning_rate = 0.0001

path = sys.argv[1]
name = sys.argv[2]

t_name = path + "/" + name + "_train_dataset.csv"
v_name = path + "/" + name + "_val_dataset.csv"

### X,Y data pre-processing ###

X_train = []
Y_train = []
with open(t_name, newline='') as csvfile:
    csvReader = csv.reader(csvfile, delimiter=',')
    for row in csvReader:
        new_row = np.array([literal_eval(feature) for feature in row])
        X_train.append(new_row[:,:-1].tolist())
        Y_train.append(new_row[-1,[-1]].tolist())

X_val = []
Y_val = []
with open(v_name, newline='') as csvfile:
    csvReader = csv.reader(csvfile, delimiter=',')
    for row in csvReader:
        new_row = np.array([literal_eval(feature) for feature in row])
        X_val.append(new_row[:,:-1].tolist())
        Y_val.append(new_row[-1,[-1]].tolist())

### define values for dimension

seq_length = 11  # [-5,+5,-4,...,-1,+1,0]
data_dim = 9     # [ins%,dels%,none%,cov,len,ins_mq,dels_mq,none_mq,bq]
output_dim = 1   # 0~2
nb_classes = 3   # ins,dels,none

# define placeholder X, Y
X = tf.placeholder(tf.float32, [None, seq_length, data_dim], name='X')
Y = tf.placeholder(tf.int32, [None, output_dim], name='Y')  # 0 ~ 2

# one-hot encoding
Y_one_hot = tf.one_hot(Y, nb_classes)  # one hot
Y_one_hot = tf.reshape(Y_one_hot, [-1, nb_classes])

hidden_dim = 128

cell = tf.contrib.rnn.BasicLSTMCell(num_units=hidden_dim, state_is_tuple=True)
outputs, _states = tf.nn.dynamic_rnn(cell, X, dtype=tf.float32)

fc_out = tf.contrib.layers.fully_connected(outputs[:,-1], hidden_dim, activation_fn=None)

W = tf.get_variable(shape=[hidden_dim, nb_classes], initializer=tf.contrib.layers.xavier_initializer(), name='weight')
b = tf.Variable(tf.constant(0., shape=[nb_classes]), name='bias')

logit = tf.matmul(fc_out, W) + b

# Cross entropy cost/loss
cost = tf.nn.softmax_cross_entropy_with_logits(logits=logit, labels=Y_one_hot)
loss = tf.reduce_mean(cost, name='loss')
train = tf.train.AdamOptimizer(learning_rate=learning_rate, name='train').minimize(loss)

hypothesis = tf.nn.softmax(logit, name = 'hypo')
prediction = tf.argmax(hypothesis, 1, name='pred')
correct_prediction = tf.equal(prediction, tf.argmax(Y_one_hot, 1), name='correct_pred')
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32), name='acc')

# using saver
saver = tf.train.Saver(max_to_keep=None)

train_dir = path + "/plot_train_0.0001_2"
val_dir = path + "/plot_val_0.0001_2"
model_dir = path + "/train_model_0.0001_2"

tensor_loss = tf.Variable(0.0)
tf.summary.scalar("loss_", tensor_loss)

# Launch graph
with tf.Session() as sess:

    sess.run(tf.global_variables_initializer())

    merged = tf.summary.merge_all()
    train_writer = tf.summary.FileWriter(train_dir, sess.graph)
    val_writer = tf.summary.FileWriter(val_dir, sess.graph)

    train_epochs = 20000
    batch_size = 128

    for epoch in range(train_epochs):

        loss_list = []

        # shuffle train data befor start 1-epoch 
        X_, Y_ = shuffle(X_train, Y_train)

        for batch in range(int(math.ceil(len(X_train)/batch_size))):

            # make batch size train data        
            x_batch, y_batch = X_[batch*batch_size : batch*batch_size+batch_size], Y_[batch*batch_size : batch*batch_size+batch_size]

            # train
            _, t_loss = sess.run([train, loss], feed_dict={X: x_batch, Y: y_batch})
            loss_list.append(t_loss)

        if epoch % 100 == 0:

            v_loss, v_acc, v_hypo = sess.run([loss, accuracy, hypothesis], feed_dict={X: X_val, Y: Y_val})
            print("Epoch: {:5}\tTrain_Loss: {:.3f}\tVal_Loss: {:.3f}\tVal_Acc: {:.2%}".format(epoch, t_loss, v_loss, v_acc))
            #print(v_loss)
            #print(v_acc)
            #print(v_hypo)

            saver.save(sess, model_dir + "/train", epoch)

            tc = np.mean(loss_list)
            t_summ = sess.run(merged, {tensor_loss: tc})
            v_summ = sess.run(merged, {tensor_loss: v_loss})

            train_writer.add_summary(t_summ, global_step=epoch)
            train_writer.flush()

            val_writer.add_summary(v_summ, global_step=epoch)
            val_writer.flush()



