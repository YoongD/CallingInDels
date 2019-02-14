import random

import numpy as np
import csv
import pysam
import pandas as pd

#### set 'path' , 'infile', 'lname'
# path = "../../../chrom/0511For10000/cover50/"

# infile = pysam.AlignmentFile(path + "realigned10000.bam", "rb")
infile = pysam.AlignmentFile("/home/yoong/chrom/20000/newtrain/realigned_NewTrain.bam", "rb")
#infile = pysam.AlignmentFile("/home/yoong/chrom/20000/test/realigned_test.bam", "rb")

max_pos = infile.lengths
#print(max_pos)

#log = pd.read_table('/home/yoong/chrom/20000/test/yLog20000_test.txt')
log = pd.read_table('/home/yoong/chrom/20000/newtrain/yLog20000_nnewtrain.txt')
log = log.sort_values(by=['position'], axis=0)


max_pos = infile.lengths
logList = list(log['position'])
""" random position pick  => #20000 """
nonNum = 0
nonList = []

while nonNum <60000 :
    s = random.randrange(9411204, max_pos[0]+1)
    logList.append(s)
    logList.sort()
    sIndex = logList.index(s)

    if sIndex == 0:
        if logList[sIndex+1] - s <= 100:
            logList.pop(sIndex)
            continue

    elif sIndex == len(logList)-1:
        if s - logList[sIndex-1] <= 100:
            logList.pop(sIndex)
            continue

    else:
        if logList[sIndex+1] - s <= 100 \
        or s - logList[sIndex-1] <= 100:
            logList.pop(sIndex)
            continue

    if s not in nonList:
        nonList.append(s)
#        nonIndexDict[s] = sIndex
        nonNum += 1




# for normalize
max_bq = 42.0
max_mq = 70.0
max_len = 30.0
max_cov = 100.0



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
# idList = []  # insert%, delete%, none%, #coverage, indel_length, mapQual, y_data
y_data = 2

### csv writer
csvfile = open('./111_inputdata_rand.csv', 'w')
wrcsv = csv.writer(csvfile)


for logIndex in range(len(log.index)):

    frontLogFlag = -1
    backLogFlag = -1

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

    """
    if (logIndex > 0) and (int(log.iloc[logIndex][1])-int(log.iloc[logIndex-1][1])) <=50:
        frontLogFlag = 1 
    if (logIndex < len(log.index)-1) and (int(log.iloc[logIndex+1][1])-int(log.iloc[logIndex][1])) <=50:
        backLogFlag = 1
    """

    forward_baseQ = [0, 0, 0, 0, 0, 0]  # for baseQuality
    backward_baseQ = [0, 0, 0, 0, 0, 0]
    f_bq = 0.
    b_bq = 0.
    idList = []

    delflag = 0
    #    malCnt = 0

    mapQual = np.array([0, 0, 0])  # mappingQuality
    indelLen = 0

    noCnt = 0
    total = 0
    ins = 0
    dels = 0

    # writeCount = -1

    nonIndelPos = 0
    nonIndelFlag = False


    ### if you want to set specific region, set start='startpos', stop='stoppos'
    for pileupcolumn in infile.pileup("chr21", int(log.iloc[logIndex][1]) - 60, int(log.iloc[logIndex][1]) + 60, truncate=True):

#        last_pileupcolumn = pileupcolumn.pos
        #        lengthDict = {0:0}
        lengthDict = dict()
        arrayt = []

        # if pileupcolumn.pos < log.iloc[logIndex][1] - 50:
        #     continue
        # elif pileupcolumn.pos > log.iloc[logIndex][1] + 50:
        #     break

        if delflag > 0:
            delflag -= 1

        for pileupread in pileupcolumn.pileups:

            ### insert
            if pileupread.indel > 0:
                mapQual[0] += pileupread.alignment.mapping_quality
                ins += 1

                ### add maximum indel_length in the column.
                if pileupread.indel in lengthDict.keys():
                    lengthDict[pileupread.indel] += 1
                else:
                    lengthDict[pileupread.indel] = 1

            ### delete
            elif pileupread.indel < 0:  # modfied: elif (X)
                mapQual[1] += pileupread.alignment.mapping_quality
                dels += 1

                ### add maximum indel_length in the column.
                if -pileupread.indel in lengthDict.keys():
                    lengthDict[-pileupread.indel] += 1
                else:
                    lengthDict[-pileupread.indel] = 1

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
                    f_bq = np.mean(forward_baseQ[:total])
            else:
                forward_baseQ[-1] /= float(pileupcolumn.n - noCnt)
                f_bq = np.mean(forward_baseQ[:-1])
        else:
            forward_baseQ[-1] = 42.
            f_bq = np.mean(forward_baseQ[:-1])

        ### if deletion, .indel is negative value of length
        #        if dels > 0:
        #        if indelLen < 0:
        #            indelLen = -indelLen
        #        indelLen = max(lengthDict.values())
        #        print(lengthDict)
        if len(lengthDict.keys()) == 0:
            indelLen = 0
        else:
            indelLen = max(zip(lengthDict.values(), lengthDict.keys()))[1]

        if indelLen > 10:
            print('malLength_pos:', pileupcolumn.pos)

        ### the first column has the highest value(42) of forward_baseQualities
        ### and the last column has the highest value(42) of backward_baseQualities

        arrayt = [(float(ins) / pileupcolumn.n), (float(dels) / pileupcolumn.n),
                  (float((pileupcolumn.n - (ins + dels))) / pileupcolumn.n), pileupcolumn.n/max_cov, indelLen/max_len]

        ### add mapping quality in list
        if (ins > 0):
            arrayt.append(float(mapQual[0] / ins) / max_mq)
        elif (dels > 0):
            arrayt.append(float(mapQual[1] / dels) / max_mq)
        else:
            arrayt.append(float(mapQual[2] / pileupcolumn.n) / max_mq)

        ### forward base quality
        if total == 0:
            arrayt.append(0.)  #42.???
        else:
            arrayt.append(f_bq/max_bq)

        ### if deletion, skip (indel_length-1) length backward.
#        if arrayt[1] >= 0.2:
        if arrayt[1] > 0 and indelLen>1:
            delflag = indelLen

        # append y_data here
        if log.iloc[logIndex][2] == 0 and arrayt[0] >= 0.2 and int(log.iloc[logIndex][3])/max_len == arrayt[4]:
            y_data = 0
        elif log.iloc[logIndex][2] == 1 and arrayt[1] >= 0.2 and int(log.iloc[logIndex][3])/max_len == arrayt[4]:
            y_data = 1
        else:
            y_data = 2


        ### Look at the side positions of current and fix the y_data
        ### (if the position is overlapped, these two flags gonna be positive)
        if frontLogFlag > -1 and y_data == 2:
            for k in range(logIndex - frontLogFlag):
                if log.iloc[logIndex - k - 1][2] == 0 and arrayt[0] >= 0.2 and int(log.iloc[logIndex - k - 1][3])/max_len == \
                        arrayt[4]:
                    y_data = 0
                    break
                elif log.iloc[logIndex - k - 1][2] == 1 and arrayt[1] >= 0.2 and int(log.iloc[logIndex - k - 1][3])/max_len == \
                        arrayt[4]:
                    y_data = 1
                    break

        if backLogFlag > -1 and y_data == 2:
            for k in range(backLogFlag - logIndex):
                if log.iloc[logIndex + k + 1][2] == 0 and arrayt[0] >= 0.2 and int(log.iloc[logIndex + k + 1][3])/max_len == \
                        arrayt[4]:
                    y_data = 0
                    break
                elif log.iloc[logIndex + k + 1][2] == 1 and arrayt[1] >= 0.2 and int(log.iloc[logIndex + k + 1][3])/max_len == \
                        arrayt[4]:
                    y_data = 1
                    break


        arrayt.append(y_data)
        idList.append(arrayt)

        if total >= 10:
            b_bq = np.mean(forward_baseQ[1:]) / max_bq
            idList[0].insert(-1, b_bq)
            wrcsv.writerow(idList[0])


            ### initialize
            del forward_baseQ[0]
            forward_baseQ.append(0)
            del idList[0]


        # if writeCount > 0:
        #     writeCount-=1

        # Initialize
        total += 1  # for count all the columns
        ins = 0  # insert count in a column
        dels = 0  # delete count in a column
        mapQual = [0, 0, 0]  # save sum of mapping quality separately ins/dels/none in each column.
        indelLen = 0  # if indel is found, save indel length
        y_data = 2  # defalut y_data = 2 (0:ins, 1:dels, 2:none)
        f_bq = 0. #42
        b_bq = 0. #42
        noCnt = 0

    """ no need to add here"""
    ### last 5 element's backward base quality adding ###
    # for i in range(4):
    #     b_bq = np.mean(forward_baseQ[(i + 1):-1])  # the last element is 0, so have to be excluded
    #     idList[i].insert(-1, b_bq)
    #
    #     # if writeCount==0:
    #     wrcsv.writerow(idList[i])
    #
    # ### the last column's backward base quality adding
    # idList[-1].insert(-1, 42.)
    # wrcsv.writerow(idList[-1])



for non_pos in nonList:


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

    nonIndelPos = 0
    nonIndelFlag = False

    last_pileupcolumn = 0

    ### if you want to set specific region, set start='startpos', stop='stoppos'
    for pileupcolumn in infile.pileup("chr21", non_pos -25, non_pos +25, truncate=True):

        lengthDict = dict()
        arrayt = []


        # if pileupcolumn.pos < non_pos - 20:
        #     continue
        # elif pileupcolumn.pos > non_pos + 20:
        #     break


        if delflag > 0:
            delflag -= 1
            continue

        last_pileupcolumn = pileupcolumn.pos


        for pileupread in pileupcolumn.pileups:

            ### insert
            if pileupread.indel > 0:
                mapQual[0] += pileupread.alignment.mapping_quality
                ins += 1

                ### add maximum indel_length in the column.
                if pileupread.indel in lengthDict.keys():
                    lengthDict[pileupread.indel] += 1
                else:
                    lengthDict[pileupread.indel] = 1

            ### delete
            elif pileupread.indel < 0:  # modfied: elif (X)
                mapQual[1] += pileupread.alignment.mapping_quality
                dels += 1

                ### add maximum indel_length in the column.
                if -pileupread.indel in lengthDict.keys():
                    lengthDict[-pileupread.indel] += 1
                else:
                    lengthDict[-pileupread.indel] = 1

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
                    f_bq = np.mean(forward_baseQ[:total])
            else:
                forward_baseQ[-1] /= float(pileupcolumn.n - noCnt)
                f_bq = np.mean(forward_baseQ[:-1])
        else:
            forward_baseQ[-1] = 42.
            f_bq = np.mean(forward_baseQ[:-1])

        ### if deletion, .indel is negative value of length
        if len(lengthDict.keys()) == 0:
            indelLen = 0
        else:
            indelLen = max(zip(lengthDict.values(), lengthDict.keys()))[1]

        if indelLen > 10:
            print('malLength_pos:', pileupcolumn.pos)

        # add detected indels but not manipulated by us
        if nonIndelPos == 0 and indelLen > 0:
            nonIndelPos = pileupcolumn.pos
            nonIndelFlag = True


        ### the first column has the highest value(42) of forward_baseQualities
        ### and the last column has the highest value(42) of backward_baseQualities

        arrayt = [(float(ins) / pileupcolumn.n), (float(dels) / pileupcolumn.n),
                  (float((pileupcolumn.n - (ins + dels))) / pileupcolumn.n), pileupcolumn.n/max_cov, indelLen/max_len]

        ### add mapping quality in list
        if (ins > 0):
            arrayt.append(float(mapQual[0] / ins)/max_mq)
        elif (dels > 0):
            arrayt.append(float(mapQual[1] / dels)/max_mq)
        else:
            arrayt.append(float(mapQual[2] / pileupcolumn.n)/max_mq)

        ### forward base quality
        if total == 0:
            arrayt.append(0.)  #42.??  -> not added anyway in the later
        else:
            arrayt.append(f_bq/max_bq)

        ### if deletion, skip (indel_length-1) length backward.
        if arrayt[1] >= 0.2:
            delflag = indelLen

        # append y_data here
        # append y_data here
        if log.iloc[logIndex][2]==0 and arrayt[0]>=0.2 and float(log.iloc[logIndex][3])/max_len == arrayt[4]:
            y_data = 0
        elif log.iloc[logIndex][2]==1 and arrayt[1]>=0.2 and float(log.iloc[logIndex][3])/max_len == arrayt[4]:
            y_data = 1
        else:
            y_data = 2


        arrayt.append(y_data)
        idList.append(arrayt)

        if total >= 10:
            b_bq = np.mean(forward_baseQ[1:])/max_bq
            idList[0].insert(-1, b_bq)
            if pileupcolumn.pos == non_pos+5 and y_data ==2:
                wrcsv.writerow(idList[0])

            elif nonIndelFlag and pileupcolumn.pos == nonIndelPos+5: #and y_data ==2:
                wrcsv.writerow(idList[0])
                nonIndelPos = 0
                nonIndelFlag = False


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
        f_bq = 0.  #42
        b_bq = 0.  #42
        noCnt = 0



infile.close()
csvfile.close()
