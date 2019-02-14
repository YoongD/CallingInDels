import numpy as np
import pysam
import pandas as pd
import tensorflow as tf
import os

bamfile = '/home/yoong/chrom/20000/test/realigned_test.bam'
#bamfile = './newtrain/realigned_NewTrain.bam'
#bamfile = "./forTest/TestRealigned.bam"

nb_classes = 3

input_size = 8
X = tf.placeholder(tf.float32, [None, input_size])
Y = tf.placeholder(tf.int32, [None, 1])  # 0 ~ 2

Y_one_hot = tf.one_hot(Y, nb_classes)  # one hot
Y_one_hot = tf.reshape(Y_one_hot, [-1, nb_classes])

hidden1_size = 128
hidden2_size = 128

"""
dropout_rate = 0.6
with tf.name_scope("layer1") as scope:
    W1 = tf.get_variable("W1", shape=[input_size, hidden1_size], initializer=tf.contrib.layers.xavier_initializer())
    b1 = tf.Variable(tf.constant(0., shape=[hidden1_size]), name='bias1')

    layer1 = tf.nn.relu(tf.matmul(X, W1) + b1)

with tf.name_scope("layer2") as scope:
    W2 = tf.get_variable("W2", shape=[hidden1_size, hidden2_size], initializer=tf.contrib.layers.xavier_initializer())
    b2 = tf.Variable(tf.constant(0., shape=[hidden2_size]), name='bias2')

    layer2 = tf.nn.relu(tf.matmul(layer1, W2) + b2)

with tf.name_scope("layer3") as scope:
    W3 = tf.get_variable("W3", shape=[hidden2_size, nb_classes], initializer=tf.contrib.layers.xavier_initializer())
    b3 = tf.Variable(tf.constant(0., shape=[nb_classes]), name='bias3')


    logits = tf.matmul(layer2, W3) + b3
    hypothesis = tf.nn.softmax(tf.nn.relu(logits))

# Cross entropy cost/loss
cost
"""




"""
prediction = tf.argmax(hypothesis, 1)
correct_prediction = tf.equal(prediction, tf.argmax(Y_one_hot, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))
"""



""" restore the ckpt file to saver """
sess = tf.Session()
sess.run(tf.global_variables_initializer())
SAVER_DIR = './resplit_0916/model/'#'./model/'
#SAVER_DIR = "./change0830/model/"
#saver = tf.train.Saver()
checkpoint_path = os.path.join(SAVER_DIR, "train")
ckpt = tf.train.get_checkpoint_state(SAVER_DIR)
#saver = tf.Saver()
saver = tf.train.import_meta_graph(SAVER_DIR+'train-1500.meta')
#if ckpt and ckpt.model_checkpoint_path:
saver.restore(sess, SAVER_DIR+'train-1500')




#print('x',tf.get_collection('input'))
X = tf.get_collection('input')[0]
Y = tf.get_collection('input')[1]

Y_one_hot = tf.one_hot(Y, nb_classes)  # one hot
Y_one_hot = tf.reshape(Y_one_hot, [-1, nb_classes])

all_vars = tf.get_collection('vars')
W1 = all_vars[0]
W2=all_vars[1]
W3=all_vars[2]

b1=all_vars[3]
b2=all_vars[4]
b3=all_vars[5]

#layer1 = tf.sigmoid(tf.matmul(X, W1)+b1)
#layer2 = tf.sigmoid(tf.matmul(layer1, W2)+b2)
#logits = tf.matmul(layer2, W3)+b3
layer1 = tf.nn.relu(tf.matmul(X, W1) + b1)
layer2 = tf.nn.relu(tf.matmul(layer1, W2) + b2)
hypothesis = tf.matmul(layer2, W3) + b3
#logits = tf.matmul(layer2, W3) + b3
#hypothesis = tf.nn.softmax(tf.nn.relu(logits))
#hypothesis = tf.get_collection('hypo')[0]

# Cross entropy cost/loss
cost_i = tf.nn.softmax_cross_entropy_with_logits(logits=hypothesis, labels=Y_one_hot)
cost = tf.reduce_mean(cost_i)
#optimizer = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(cost)

prediction = tf.argmax(hypothesis, 1)

correct_prediction = tf.equal(prediction, tf.argmax(Y_one_hot, 1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))


#write_cost = tf.summary.merge_all()
#writer_test = tf.summary.FileWriter('./logs_test/plot_test', sess.graph)

#print("w1",W1,"w2",W2,"w3",W3, accuracy, cost)

""" making input feature here """
infile = pysam.AlignmentFile(bamfile, "rb")
max_pos = infile.lengths

#log = pd.read_table('yLog20000_nonoverlap_test.txt')
#log = pd.read_table('./newtrain/yLog20000_nnewtrain.txt')
log = pd.read_table('/home/yoong/chrom/20000/test/yLog20000_test.txt')
log = log.sort_values(by=['position'], axis=0)

## added
forward_baseQ = [0, 0, 0, 0, 0, 0]  # for baseQuality
backward_baseQ = [0, 0, 0, 0, 0, 0]
f_bq = 0.
b_bq = 0.
idList = []

delflag = 0

mapQual = np.array([0, 0, 0])  # mappingQuality
indelLen = 0

noCnt = 0
total = 0
ins = 0
dels = 0

y_data = 2


totalCnt = 0
accu = 0


resultDict = {"true_pos":0, "true_neg":0, "false_pos":0, "false_neg":0}


indelResult = [['pos','y_pred','data']]

### csv writer
#csvfile = open('./nonoverlap2/nonlaped2_inputdata20000.csv', 'w')
#wrcsv = csv.writer(csvfile)
logIndex = 0

#for logIndex in range(len(log.index)):



# for normalize
max_bq = 42.0
max_mq = 70.0
max_len = 30.0
max_cov = 100.0



frontLogFlag = -1
backLogFlag = -1

forward_baseQ = [0, 0, 0, 0, 0, 0]  # for baseQuality
backward_baseQ = [0, 0, 0, 0, 0, 0]
f_bq = 0.
b_bq = 0.
idList = []

delflag = 0

mapQual = np.array([0, 0, 0])  # mappingQuality
indelLen = 0

noCnt = 0
total = 0
ins = 0
dels = 0


### if you want to set specific region, set start='startpos', stop='stoppos'
for pileupcolumn in infile.pileup("chr21", truncate=True):

    frontLogFlag = -1
    backLogFlag = -1

    #if logIndex == len(log.index)-1:
    #    pass

    if log.iloc[logIndex][1] - 45 < pileupcolumn.pos < log.iloc[logIndex][1] + 45:
        for i in range(logIndex - 1, -1, -1):
            if (int(log.iloc[logIndex][1]) - int(log.iloc[i][1])) >= 50:
                break
            else:
                frontLogFlag = i

        for i in range(logIndex + 1, len(log.index)):
            if (int(log.iloc[i][1]) - int(log.iloc[logIndex][1])) >= 50:
                break
            else:
                backLogFlag = i

    elif pileupcolumn.pos >= log.iloc[logIndex][1] + 45 and logIndex<len(log.index)-1\
        and log.iloc[logIndex+1][1] - 50 <= pileupcolumn.pos <= log.iloc[logIndex+1][1] + 50:
        logIndex += 1
#        frontLogFlag = -1
#        backLogFlag = -1
        for i in range(logIndex - 1, -1, -1):
            if (int(log.iloc[logIndex][1]) - int(log.iloc[i][1])) >= 50:
                break
            else:
                frontLogFlag = i

        for i in range(logIndex + 1, len(log.index)):
            if (int(log.iloc[i][1]) - int(log.iloc[logIndex][1])) >= 50:
                break
            else:
                backLogFlag = i


    lengthDict = dict()
    arrayt = []
#    if pileupcolumn.pos < log.iloc[logIndex][1] - 50:
#        continue
#    elif pileupcolumn.pos > log.iloc[logIndex][1] + 50:
#        break

    if delflag > 0:
        delflag -= 1
#        print(pileupcolumn.pos)
        continue

    for pileupread in pileupcolumn.pileups:

        ### insert
        if pileupread.indel > 0:
            mapQual[0] += pileupread.alignment.mapping_quality
            ins += 1

            ### add maximum indel_length in the column.
            #                if pileupread.indel > indelLen:
            #                    indelLen = pileupread.indel
            if pileupread.indel in lengthDict.keys():
                lengthDict[pileupread.indel] += 1
            else:
                lengthDict[pileupread.indel] = 1

        ### delete
        elif pileupread.indel < 0:  # modfied: elif (X)
            mapQual[1] += pileupread.alignment.mapping_quality
            dels += 1

            ### add maximum indel_length in the column.
            #                if abs(pileupread.indel) > indelLen:
            #                    indelLen = -pileupread.indel
            delLen = -(pileupread.indel)
            if delLen in lengthDict.keys():
                lengthDict[delLen] += 1
            else:
                lengthDict[delLen] = 1
            

        ### not indel
        elif pileupread.indel == 0:
            mapQual[2] += pileupread.alignment.mapping_quality
            indelLen = 0

        ### add baseQuality
        if pileupread.is_del:
            noCnt += 1
            continue

        """ when deletion occured, back_base quality is needed to be the futher more back."""
        if total < 6:
            forward_baseQ[total] += pileupread.alignment.query_qualities[pileupread.query_position]
        else:
            forward_baseQ[-1] += pileupread.alignment.query_qualities[pileupread.query_position]

    ### for base quality
    if pileupcolumn.n - noCnt != 0:
        if total < 6:
            forward_baseQ[total] /= float(pileupcolumn.n - noCnt)
            if total != 0:
                f_bq = np.mean(forward_baseQ[:total]) /max_bq
        else:
            forward_baseQ[-1] /= float(pileupcolumn.n - noCnt)
            f_bq = np.mean(forward_baseQ[:-1]) /max_bq
    else:
        forward_baseQ[-1] = max_bq
        f_bq = np.mean(forward_baseQ[:-1])/max_bq  # must be normalized!!!

    ### if deletion, .indel is negative value of length
    #        if indelLen < 0:
    #            indelLen = -indelLen

    if len(lengthDict.keys()) == 0:
        indelLen = 0
    else:
        indelLen = max(zip(lengthDict.values(), lengthDict.keys()))[1]


    #        if indelLen > 10:
    #            print('malLength_pos:', pileupcolumn.pos)

    ### the first column has the highest value(42) of forward_baseQualities
    ### and the last column has the highest value(42) of backward_baseQualities

    arrayt = [(float(ins) / pileupcolumn.n), (float(dels) / pileupcolumn.n),
              (float((pileupcolumn.n - (ins + dels))) / pileupcolumn.n), pileupcolumn.n/max_cov, indelLen /max_len]

    ### add mapping quality in list
    if (ins > 0):
        arrayt.append((float(mapQual[0]) / ins)/max_mq)
    elif (dels > 0):
        #print(indelLen)
        arrayt.append((float(mapQual[1]) / dels)/max_mq)
    else:
        arrayt.append((float(mapQual[2]) / pileupcolumn.n)/max_mq)

    ### forward base quality
    if total == 0:
        arrayt.append(0.)   #not 1. # not 42.
    else:
        arrayt.append(f_bq)

    ### if deletion, skip (indel_length-1) length backward.
    if arrayt[1] > 0 and indelLen>1:
        delflag = indelLen

    # append y_data here
    if log.iloc[logIndex][2] == 0 and arrayt[0] >= 0.2 and log.iloc[logIndex][3] == int(arrayt[4]*max_len):
        y_data = 0
    elif log.iloc[logIndex][2] == 1 and arrayt[1] >= 0.2 and log.iloc[logIndex][3] == int(arrayt[4]*max_len):
        y_data = 1
    else:
        y_data = 2




    ### Look at the side positions of current and fix the y_data
    ### (if the position is overlapped, these two flags gonna be positive)
    if frontLogFlag > -1 and y_data == 2:
        for k in range(logIndex - frontLogFlag):
            if log.iloc[logIndex - k - 1][2] == 0 and arrayt[0] >= 0.2 and float(
                    log.iloc[logIndex - k - 1][3]) /max_len == arrayt[4]:
                y_data = 0
                break
            elif log.iloc[logIndex - k - 1][2] == 1 and arrayt[1] >= 0.2 and float(
                    log.iloc[logIndex - k - 1][3]) /max_len == arrayt[4]:
                y_data = 1
                break

    if backLogFlag > -1 and y_data == 2:
        for k in range(backLogFlag - logIndex):
            if log.iloc[logIndex + k + 1][2] == 0 and arrayt[0] >= 0.2 and float(
                    log.iloc[logIndex + k + 1][3]) /max_len == arrayt[4]:
                y_data = 0
                break
            elif log.iloc[logIndex + k + 1][2] == 1 and arrayt[1] >= 0.2 and float(
                    log.iloc[logIndex + k + 1][3]) /max_len == arrayt[4]:
                y_data = 1
                break

    
#    if y_data == 0:
#        print('y:0')
#    elif y_data == 1:
#        print('y:1')
    
    arrayt.append(y_data)
    idList.append(arrayt)



    if total >= 10:
# change baseQ max to 60.0
#        b_bq = np.mean(forward_baseQ[1:])/42.0
        b_bq = np.mean(forward_baseQ[1:])/max_bq
        idList[0].insert(-1, b_bq)
        #            wrcsv.writerow(idList[0])
        accu += sess.run(accuracy, feed_dict={X:[idList[0][:-1]] ,Y: [[idList[0][-1]]]})
#        y_pred = sess.run(prediction, feed_dict={X:[idList[0][:-1]]})
#        print(y_pred)
        totalCnt += 1
        if total % 500 == 0:
            print("test_Acc: ", float(accu) / totalCnt)
#            print("test_Acc: {:.2%}".format(accu/totalCnt))
        #                print('testAccuracy: %.2f'%(accu/totalCnt))

        """ sess.run with only x and compare its prediction with y_data in here """
        y_pred = sess.run(prediction, feed_dict={X: [idList[0][:-1]]})

#        if y_pred[0] == idList[0][-1]:
#            print(pileupcolumn.pos,'  data: ',idList[0],' y_pred: ',y_pred)
#            indelResult.append([pileupcolumn.pos,y_pred[0],idList[0]])
        
#        if y_pred[0] != idList[0][-1]:
#            print(pileupcolumn.pos,'  data: ',idList[0],' y_pred: ',y_pred)
#            indelResult.append([pileupcolumn.pos,y_pred[0],idList[0]])
#            print('here again ', logIndex)
        if idList[0][2] <= 0.8 and idList[0][-1] == 2:
            print(pileupcolumn.pos-5, '  ', idList[0])
#            indelResult.append([pileupcolumn.pos,idList[0][-1],y_pred[0]])

#        if y_data != 2 :
#            indelResult.append([y_data, y_pred])
#            indelResult.append(str(y_data)+'\t'+str(y_pred))
#        print('y_data: ', y_data, ' y_pred: ',y_pred)

        #            print(y_pred)
        if y_pred[0] == idList[0][-1]:
            if idList[0][-1] == 2:
                # true-negative
                resultDict['true_neg'] += 1
            else:
                # true-positive
                resultDict['true_pos'] += 1

        else:
            print(pileupcolumn.pos,'  data: ',idList[0],' y_pred: ',y_pred)
            indelResult.append([pileupcolumn.pos,y_pred[0],idList[0]])

            if idList[0][-1] == 2:
                # false_positive
                resultDict['false_pos'] += 1
            else:
                # false_negative
                resultDict['false_neg'] += 1


    if total>=5:
        ### initialize
        del forward_baseQ[0]
        forward_baseQ.append(0)
        del idList[0]


    # Initialize
    total += 1  # for count all the columns
    ins = 0  # insert count in a column
    dels = 0  # delete count in a column
    mapQual = [0, 0, 0]  # save sum of mapping quality separately ins/dels/none in each column.
    indelLen = 0  # if indel is found, save indel length
    y_data = 2  # defalut y_data = 2 (0:ins, 1:dels, 2:none)
    f_bq = 1.  # not 42.
    b_bq = 1.  # not 42.
    noCnt = 0


"""
### last 5 element's backward base quality adding ###
for i in range(4):
# change baseQ max to 60.0
#    b_bq = np.mean(forward_baseQ[(i + 1):-1])/42.0  # the last element is 0, so have to be excluded
    b_bq = np.mean(forward_baseQ[(i + 1):-1])/max_bq
    idList[i].insert(-1, b_bq)
    #        wrcsv.writerow(idList[i])
    accu += sess.run(accuracy, feed_dict={X: [idList[i][:-1]], Y: [[idList[i][-1]]]})
    totalCnt += 1

    # sess.run with only x and compare its prediction with y_data in here 
    y_pred = sess.run(prediction, feed_dict={X: [idList[i][:-1]]})

#    y_pred = sess.run(prediction, feed_dict={X: [idList[0][:-1]]})
#    if idList[i][-1] != 2 :
#        indelResult.append([idList[i][-1], y_pred])
#        indelResult.write(str(idList[i][-1])+ '\t'+str(y_pred))

    if y_pred[0] == idList[i][-1]:
        if y_data == 2:
            # true-negative
            resultDict['true_neg'] += 1
        else:
            # true-positive
            resultDict['true_pos'] += 1

    else:
        if y_data == 2:
            # false_positive
            resultDict['false_pos'] += 1
        else:
            # false_negative
            resultDict['false_neg'] += 1


#        print('testAccuracy: %.2f' % (accu / totalCnt))
### the last column's backward base quality adding
#idList[-1].insert(-1, 42.)
idList[-1].insert(-1, 1.)  #max value have to be 1!!!
#    wrcsv.writerow(idList[-1])
accu += sess.run(accuracy, feed_dict={X:[idList[-1][:-1]] ,Y: [[idList[-1][-1]]]})
totalCnt += 1

# sess.run with only x and compare its prediction with y_data in here 
y_pred = sess.run(prediction, feed_dict={X:[idList[-1][:-1]]})

#if idList[-1][-1] != 2:
#    indelResult.append([idList[-1][-1],y_pred])
#    indelResult.write(str(idList[-1][-1])+'\t'+str(y_pred))

#    print(y_pred)
if y_pred[0] == y_data:
    if y_data == 2:
        # true-negative
        resultDict['true_neg'] += 1
    else:
        # true-positive
        resultDict['true_pos'] += 1

else:
    if y_data == 2:
        # false_positive
        resultDict['false_pos'] += 1
    else:
        # false_negative
        resultDict['false_neg'] += 1
#    print('testAccuracy: %.2f' % (accu / totalCnt))
#    malCnt = 0

# print(mal_list)
# print('mal_size: ',len(mal_list))s
#print('total accuracy: ', accu/float(totalCnt))
#print("total_Accuracy: {:.2%}".format(accu/totalCnt))
"""

print("total_Accuracy: ", float(accu) / totalCnt)
infile.close()
#csvfile.close()

#indel_f = open("/home/yoong/chrom/20000/test/result/indels.txt",'w')
#indel_f = pd.DataFrame(indelResult)
#indel_f.to_csv('./result/indels.csv')

#indel_f.close()

tname = "./result/testResult.txt"
txtfile = open(tname, 'w')

txtfile.write('total_accuracy:{:.2%}\n'.format(accu/totalCnt))
txtfile.write("TP\tTN\tFP\tFN\n")
txtfile.write(str(resultDict["true_pos"])+'\t'+str(resultDict['true_neg'])+'\t'+
              str(resultDict['false_pos'])+'\t'+str(resultDict['false_neg'])+
              "\n")
txtfile.close()

resultfile = open('./result/mis_indels.txt', 'w')
for row in indelResult:
    resultfile.write(str(row)+'\n')
resultfile.close()
