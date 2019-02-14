import numpy as np
import pysam
import pandas as pd
import tensorflow as tf
import os
import sys

path_model = sys.argv[1]
path_test = sys.argv[2]
name = sys.argv[3]
num = int(sys.argv[4]) # 20000 or 40000 now
SAVER_DIR = sys.argv[5] #"./train_model/"
SAVER_DIR = path_model + SAVER_DIR

f_name = path_test + "/hg19_chr21_" + name + ".fa"
bamfile = path_test + "/realigned_" + name + ".bam"
logfile = path_test + "/ylLog_" + name + "_nucleo.txt"
tl_name = path_test + "/testLog_" + name + ".txt"

### max value for normalize
indel_max_length = 30.0
coverage_max = 100.0
mappingQ_max = 70.0
baseQ_max = 42.0

indel_range = 60

tlfile = open(tl_name, 'w')

### read log and bam file
log = pd.read_table(logfile)
infile = pysam.AlignmentFile(bamfile, "rb")

""" restore the ckpt file to saver """
sess = tf.Session()

saver = tf.train.import_meta_graph(SAVER_DIR+'/train-100.meta')
saver.restore(sess, SAVER_DIR+'/train-100')

X = tf.get_default_graph().get_tensor_by_name("X:0")
Y = tf.get_default_graph().get_tensor_by_name("Y:0")

xList = []
yList = []

accu = 0
totalCnt = 0

frontLogFlag = 0
backLogFlag = 0

show_progress = 0
def printProgress (iteration, total, prev_percent, accu_mid, prefix = '', suffix = '', decimals = 1, barLength = 100):
    global show_progress

    formatStr = "{0:." + str(decimals) + "f}"
    percent = formatStr.format(100 * (iteration / float(total)))
    accu_mid = "{0:.7f}".format(accu_mid)

    filledLength = int(round(barLength * iteration / float(total)))
    bar = '#' * filledLength + '-' * (barLength - filledLength)
    if percent != prev_percent:
        show_progress = 1
    else:
        show_progress = 0

    if show_progress:
        sys.stdout.write('\r%s |%s| %s%s %s %s\n' % (prefix, bar, percent, '%', suffix, accu_mid))
        print("true_neg: ", true_neg, " / false_neg: ", false_neg, " / false_pos: ", false_pos, " / true_pos: ", true_pos, "\n")

    if iteration == total:
        sys.stdout.write('\n')
        sys.stdout.flush()

    return percent

def set_nList():

    global start_n, end_n

    with open(f_name) as f:
        f.readline()
        nList = list(f.read().upper().replace('\n', ''))
        
    total = len(nList)
    
    # start after N's sequence
    start_n = 0
    while True:
        if nList[start_n] != 'N': break
        else:
            start_n += 1
    
    end_n = total-1
    while True:
        if nList[end_n] != 'N': break
        else:
            end_n -= 1

start_n, end_n = (0,0)
set_nList()
process_total = end_n - start_n

out_of_log = 0
error_cnt = 0

true_pos = 0
true_neg = 0
false_pos = 0
false_neg = 0

formatStr = "{0:." + str(1) + "f}"
prev_percent = formatStr.format(100 * (0 / float(process_total)))
for pileupcolumn in infile.pileup("chr21",truncate=True):

    ### Initialize
    ins = 0  # insert count in a column
    dels = 0  # delete count in a column
    noCnt = 0

    mapQual = [0, 0, 0]  # save sum of mapping quality separately ins/dels/none in each column.
    baseQ = 0.0

    lengthDict = dict()
    currPos = pileupcolumn.pos
    arrayt = []

    if totalCnt == 0:
        accu_mid = 0
    else:
        accu_mid = float(accu) / float(totalCnt)

    prev_percent = printProgress(currPos-start_n, process_total, prev_percent, accu_mid, 'Process', 'Complete', 1, 50)

    coverage = len(pileupcolumn.pileups)
    if coverage == 0:
        idList = []
        continue

    ### range(frontLogFlag, backLogFlag)
    while (out_of_log==0) and (currPos-indel_range > log.iloc[frontLogFlag][1]):
        frontLogFlag += 1
        if frontLogFlag == num:
            #frontLogFlag -= 1
            out_of_log = 1

    backLogFlag = frontLogFlag
    while (out_of_log==0) and (currPos+indel_range > log.iloc[backLogFlag][1]):
        backLogFlag += 1
        if backLogFlag == num:
            #backLogFlag -= 1
            break

    for pileupread in pileupcolumn.pileups:

        ### insert
        if pileupread.indel > 0:
            mapQual[0] += pileupread.alignment.mapping_quality
            ins += 1

            if pileupread.indel in lengthDict.keys():
                lengthDict[pileupread.indel] += 1
            else:
                lengthDict[pileupread.indel] = 1

        ### delete
        elif pileupread.indel < 0:
            mapQual[1] += pileupread.alignment.mapping_quality
            dels += 1

            if -pileupread.indel in lengthDict.keys():
                lengthDict[pileupread.indel] += 1
            else:
                lengthDict[pileupread.indel] = 1

        ### not indel
        elif pileupread.indel == 0:
            mapQual[2] += pileupread.alignment.mapping_quality
            noCnt += 1

        ### add baseQuality
        if pileupread.is_del or pileupread.is_refskip:
            baseQ += baseQ_max
        else:
            baseQ += pileupread.alignment.query_qualities[pileupread.query_position]

    ### for base quality
    baseQ /= float(coverage)

    ### setting indel length
    if len(lengthDict.keys()) == 0:
        indelLen = 0
    else:
       indelLen = max(zip(lengthDict.values(), lengthDict.keys()))[1]

    ### add mapping quality in list
    if (ins > 0) and (dels > 0):
        print("ins and dels : ", pileupcolumn.pos)
        
    if ins>0:
         ins_mq = float(mapQual[0] / ins)
    else:
         ins_mq = 0.

    if dels>0:
         del_mq = float(mapQual[1] / dels)
    else:
         del_mq = 0.
    
    if noCnt>0:
         none_mq = float(mapQual[2] / noCnt)
    else:
         none_mq = 0.

    arrayt = [ ins/coverage, dels/coverage, (coverage - (ins + dels))/coverage, 
                  coverage/coverage_max, abs(indelLen)/indel_max_length, ins_mq/mappingQ_max, 
                     del_mq/mappingQ_max, none_mq/mappingQ_max, baseQ/baseQ_max ]

    # defalut y_data = 2 (0:ins, 1:dels, 2:none)
    y_data = 2
    if (out_of_log==0) and ((arrayt[0]>=0.2) or (arrayt[1]>=0.2)):

        for currLogFlag in range(frontLogFlag, backLogFlag):

            if log.iloc[currLogFlag][1]-indel_range <= currPos <= log.iloc[currLogFlag][1]+indel_range:

                if indelLen == (log.iloc[currLogFlag][3]) * pow(-1, log.iloc[currLogFlag][2]):
                    if log.iloc[currLogFlag][2]==0:
                        y_data = 0
                        break
                    elif log.iloc[currLogFlag][2]==1:
                        y_data = 1
                        break


    xList.append(arrayt)
    yList.append([y_data])

    ## 앞,뒤로 5개의 arrayt가 나올 수 없는 구간은 test하지 못함
    if len(xList) == 11:
        
        tmpX = [xList[0], xList[10], xList[1], xList[9], xList[2], 
                xList[8], xList[3], xList[7], xList[4], xList[6], xList[5]]
        
        totalCnt += 1
        #accu_tmp = sess.run('acc:0', feed_dict={X:[tmpX], Y:[yList[5]]})
        pred = sess.run('pred:0', feed_dict={X:[tmpX]})
        #accu += accu_tmp

        if pred != yList[5]:
            print("false : ", pred, " : ", yList, "/", currPos, "/", accu/totalCnt, "/", accu, "/", totalCnt, "\n", arrayt)
            tl_content = str(currPos) + ", " + str(yList[5]) + ", " + str(pred) + "\n" + str(arrayt)
            tlfile.write(tl_content)
        else:
            accu += 1

        if yList[5] == [2]:
            if pred == [2]: true_neg  += 1
            else:           false_neg += 1

        else:
            if pred == [2]: false_pos += 1
            else:           true_pos  += 1

        del xList[0]
        del yList[0]

print("total_Accuracy:" + str(accu/totalCnt))
infile.close()
tlfile.close()
