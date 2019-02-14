from sklearn.utils import shuffle

import tensorflow as tf
import numpy as np
import math

model_path = './resplit_0916/model/'
writer_path = './resplit_0916/logs/'


tf.set_random_seed(777)  # for reproducibility
learning_rate = 0.0001

#xy = np.loadtxt('./oldData/nonoverlap3_trainingData_20000_re.csv', delimiter=',', dtype=np.float32)
#xy = np.loadtxt('./newtrain/notrandom_trainingData_20000.csv', delimiter=',', dtype=np.float32)
#val = np.loadtxt('./test/newtrain_trainingData_test.csv', delimiter=',', dtype=np.float32)


#xy = np.loadtxt('/home/yoong/chrom/20000/newtrain/newtrain_trainingData_20000_3_2.csv', delimiter=',', dtype=np.float32)
in_data = np.loadtxt('./input_split_file/y_0.csv', delimiter=',', dtype=np.float32)
del_data = np.loadtxt('./input_split_file/y_1.csv', delimiter=',', dtype=np.float32)
none_data = np.loadtxt('./input_split_file/y_2.csv', delimiter=',', dtype=np.float32)

# in_data = xy[:20000]
# del_data = xy_1
# none_data = xy[20000:]
in_data = shuffle(in_data)
del_data = shuffle(del_data)
none_data = shuffle(none_data)

none_split_size = int((len(in_data)+len(del_data))/2)
none_data = none_data[:none_split_size]

in_split = int(len(in_data)/4)
del_split = int(len(del_data)/4)
none_split = int(none_split_size/4)

print(in_split,'/',del_split,'/',none_split)

train_data = in_data[:(in_split*3)]
train_data = np.append(train_data, del_data[:(del_split*3)],0)
train_data = np.append(train_data, none_data[:(none_split*3)],0)


val_data = in_data[(in_split*3):]
val_data = np.append(val_data, del_data[(del_split*3):], 0)
val_data = np.append(val_data, none_data[(none_split*3):], 0)

train_data = shuffle(train_data)
val_data = shuffle(val_data)

X_train = train_data[:, 0:-1]
Y_train = train_data[:, [-1]]

X_val = val_data[:, 0:-1]
Y_val = val_data[:, [-1]]



#x_data = xy[:, 0:-1]
#y_data = xy[:, [-1]]
#X_train, Y_train = shuffle(x_data,y_data)


# indel_data = xy[:20000]
# none_data = xy[20000:]
# indel_data = shuffle(indel_data)
# none_data = shuffle(none_data)
#
# train_data = indel_data[:15000]
# train_data = np.append(train_data, none_data[:15000],0)
#
# val_data = indel_data[15000:]
# val_data = np.append(val_data, none_data[15000:], 0)
#
# train_data = shuffle(train_data)
# val_data = shuffle(val_data)
#
# X_train = train_data[:, 0:-1]
# Y_train = train_data[:, [-1]]
#
# X_val = val_data[:, 0:-1]
# Y_val = val_data[:, [-1]]


nb_classes = 3  # 0 ~ 2

accuracy_list = []

input_size = 8
X = tf.placeholder(tf.float32, [None, input_size])
Y = tf.placeholder(tf.int32, [None, 1])  # 0 ~ 2

Y_one_hot = tf.one_hot(Y, nb_classes)  # one hot
Y_one_hot = tf.reshape(Y_one_hot, [-1, nb_classes])

hidden1_size = 128
hidden2_size = 128

#dropout_rate = 0.6
with tf.name_scope("layer1") as scope:
    W1 = tf.get_variable("W1", shape=[input_size, hidden1_size], initializer=tf.contrib.layers.xavier_initializer())
    b1 = tf.Variable(tf.constant(0.0, shape=[hidden1_size]), name='bias1')

    layer1 = tf.nn.relu(tf.matmul(X, W1) + b1)

with tf.name_scope("layer2") as scope:
    W2 = tf.get_variable("W2", shape=[hidden1_size, hidden2_size], initializer=tf.contrib.layers.xavier_initializer())
    b2 = tf.Variable(tf.constant(0.0, shape=[hidden2_size]), name='bias2')

    layer2 = tf.nn.relu(tf.matmul(layer1, W2) + b2)

with tf.name_scope("layer3") as scope:
    W3 = tf.get_variable("W3", shape=[hidden2_size, nb_classes], initializer=tf.contrib.layers.xavier_initializer())
    b3 = tf.Variable(tf.constant(0.0, shape=[nb_classes]), name='bias3')

    hypothesis = tf.matmul(layer2, W3) + b3

tf.add_to_collection('hypo', hypothesis)
tf.add_to_collection('input', X)
tf.add_to_collection('input', Y)

tf.add_to_collection('vars', W1)
tf.add_to_collection('vars', W2)
tf.add_to_collection('vars', W3)
tf.add_to_collection('vars', b1)
tf.add_to_collection('vars', b2)
tf.add_to_collection('vars', b3)

# using saver
saver = tf.train.Saver(max_to_keep=None)

# Cross entropy cost/loss
cost = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=hypothesis, labels=Y_one_hot))
optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)

prediction = tf.argmax(hypothesis, 1)
correct_prediction = tf.equal(prediction, tf.argmax(Y_one_hot, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

tensor_cost = tf.Variable(0.0)
tf.summary.scalar("cost", tensor_cost)

# Launch graph
with tf.Session() as sess:

    sess.run(tf.global_variables_initializer())

    write_cost = tf.summary.merge_all()
    writer_train = tf.summary.FileWriter(writer_path+'plot_train', sess.graph)
    writer_val = tf.summary.FileWriter(writer_path+'plot_val', sess.graph)

    train_epochs = 200000
    batch_size = 128

    for epoch in range(train_epochs):
        
        # shuffle train data befor start 1-epoch 
        X_, Y_ = shuffle(X_train, Y_train)

        # initialize train_costs list
        train_costs = []

        for batch in range(int(math.ceil(len(X_)/batch_size))):

            # make batch size train data        
            x_batch, y_batch = X_[batch*batch_size : batch*batch_size+batch_size], Y_[batch*batch_size : batch*batch_size+batch_size]

            # train
            _, train_loss = sess.run([optimizer, cost], feed_dict={X: x_batch, Y: y_batch})
            train_costs.append(train_loss)

        if epoch % 100 == 0:

            val_loss, acc = sess.run([cost, accuracy], feed_dict={X: X_val, Y: Y_val})
            print("Epoch: {:5}\tVal_Loss: {:.3f}\tVal_Acc: {}".format(epoch, val_loss, acc))
#            print("Epoch: {:5}\tVal_Loss: {:.3f}\tVal_Acc: {:.2%}".format(epoch, val_loss, acc))
            saver.save(sess, model_path+'train', epoch)
         
            tc = np.sum(train_costs)/len(train_costs)
            train_cost = sess.run(write_cost, {tensor_cost: tc})
            val_cost = sess.run(write_cost, {tensor_cost: val_loss})
            writer_train.add_summary(train_cost, global_step=epoch)
            writer_train.flush()
            print(tc)
            writer_val.add_summary(val_cost, global_step=epoch)
            writer_val.flush()
            print(val_loss)

