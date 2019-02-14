import random

import numpy as np
import csv
import pysam
import pandas as pd

#### set 'path' , 'infile', 'lname'
# path = "../../../chrom/0511For10000/cover50/"

# infile = pysam.AlignmentFile(path + "realigned10000.bam", "rb")
#infile = pysam.AlignmentFile("/home/yoong/chrom/20000/newtrain/realigned_NewTrain.bam", "rb")
infile = pysam.AlignmentFile("/home/yoong/chrom/20000/test/realigned_test.bam", "rb")

max_pos = infile.lengths
#print(max_pos)

log = pd.read_table('/home/yoong/chrom/20000/test/yLog20000_test.txt')
log = log.sort_values(by=['position'], axis=0)


max_pos = infile.lengths
logList = list(log['position'])
""" random position pick  => #20000 """
nonNum = 0
nonList = []
#nonIndexDict = dict()
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


ydataList = []

# Loading position of indel in log file
# lList = np.array([])
# lCount = 0

"""
#lname = path + "yLog10000.txt"  ###logfile name
lname = "../../../chrom/0511For10000/yLog10000.txt"  ###logfile name
with open(lname) as log:
    log.readline()
    while True:
        tmp = list(log.readline().split('\t'))
        if tmp == ['']:    break
        lList.__add__(np.array([int(tmp[1]), int(tmp[2])]))  # positon([1]), ins/dels([2])
        lCount += 1  # how many rows==how many indels

"""

"""
### open text for log
posFile = "input_position10000.txt"
pos_file = open(posFile, 'w')
posArr = []
"""
# print(lList[:,0])
### csv writer
csvfile = open('./test/newTrain_inputdata_test.csv', 'w')
wrcsv = csv.writer(csvfile)

# mal_list = []

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

    writeCount = -1

    nonIndelPos = 0
    nonIndelFlag = False


    ### if you want to set specific region, set start='startpos', stop='stoppos'
    for pileupcolumn in infile.pileup("chr21", int(log.iloc[logIndex][1]) - 55, int(log.iloc[logIndex][1]) + 55):

#        last_pileupcolumn = pileupcolumn.pos
        #        lengthDict = {0:0}
        lengthDict = dict()
        arrayt = []

        if pileupcolumn.pos < log.iloc[logIndex][1] - 50:
            continue
        elif pileupcolumn.pos > log.iloc[logIndex][1] + 50:
            break

        if delflag > 0:
            delflag -= 1
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
                  (float((pileupcolumn.n - (ins + dels))) / pileupcolumn.n), pileupcolumn.n, indelLen]

        ### add mapping quality in list
        if (ins > 0):
            arrayt.append(float(mapQual[0] / ins))
        elif (dels > 0):
            arrayt.append(float(mapQual[1] / dels))
        else:
            arrayt.append(float(mapQual[2] / pileupcolumn.n))

        ### forward base quality
        if total == 0:
#            arrayt.append(60.)
            arrayt.append(42.)
        else:
            arrayt.append(f_bq)

        ### if deletion, skip (indel_length-1) length backward.
#        if arrayt[1] >= 0.2:
        if arrayt[1] > 0 and indelLen>1:
            delflag = indelLen

        # append y_data here
        if log.iloc[logIndex][2] == 0 and arrayt[0] >= 0.2 and int(log.iloc[logIndex][3]) == arrayt[4]:
            y_data = 0
        elif log.iloc[logIndex][2] == 1 and arrayt[1] >= 0.2 and int(log.iloc[logIndex][3]) == arrayt[4]:
            y_data = 1
        else:
            y_data = 2

        """
        ## modify y_data from the logIndex
        if arrayt[2]>=0.8 and y_data==2:
            for k in range(-4,5):
                if k == 0:
                    continue
                if log.iloc[logIndex+k][2]==0 and arrayt[0]>=0.2 and float(log.iloc[logIndex+k][3])/10.0 == arrayt[4]:
                    y_data = 0
                    break
                elif log.iloc[logIndex+k][2]==1 and arrayt[1]>=0.2 and float(log.iloc[logIndex+k][3])/10.0 == arrayt[4]:
                    y_data = 1
                    break
        """

        ### Look at the side positions of current and fix the y_data
        ### (if the position is overlapped, these two flags gonna be positive)
        if frontLogFlag > -1 and y_data == 2:
            for k in range(logIndex - frontLogFlag):
                if log.iloc[logIndex - k - 1][2] == 0 and arrayt[0] >= 0.2 and int(log.iloc[logIndex - k - 1][3]) == \
                        arrayt[4]:
                    y_data = 0
                    break
                elif log.iloc[logIndex - k - 1][2] == 1 and arrayt[1] >= 0.2 and int(log.iloc[logIndex - k - 1][3]) == \
                        arrayt[4]:
                    y_data = 1
                    break

        if backLogFlag > -1 and y_data == 2:
            for k in range(backLogFlag - logIndex):
                if log.iloc[logIndex + k + 1][2] == 0 and arrayt[0] >= 0.2 and int(log.iloc[logIndex + k + 1][3]) == \
                        arrayt[4]:
                    y_data = 0
                    break
                elif log.iloc[logIndex + k + 1][2] == 1 and arrayt[1] >= 0.2 and int(log.iloc[logIndex + k + 1][3]) == \
                        arrayt[4]:
                    y_data = 1
                    break

        """            
        ### add y_data in list.
        for i in range(lCount):
            if pileupcolumn.pos + 1 == lList[i][0]:  ### if next position is in log data
                if lList[i][1] == 0:  # insert
                    y_data = 0
                    break
                elif lList[i][1] == 1:  # delete
                    y_data = 1
                    break
        ### not indels is a default set to 2.
        """
        """
        if arrayt[2] <=0.8 and y_data ==2:
            if pileupcolumn.pos in mal_list:
                continue
            else:
                mal_list.append(pileupcolumn.pos)
                malCnt += 1
        """
        #            print('pos',pileupcolumn.pos,' y_index',log.iloc[logIndex][0] ,' in: ', arrayt[0], ' del: ', arrayt[1], ' ', arrayt[4], 'vs', float(log.iloc[logIndex][3])/10.0, ' T/F ', arrayt[4]==float(log.iloc[logIndex][3])/10.0)
#        nonFlag = 0
        if y_data !=2 and writeCount == -1 and (pileupcolumn.pos not in ydataList):
#            if pileupcolumn.pos in ydataList:
#                nonFlag = 1
#            else:
            ydataList.append(pileupcolumn.pos)
#                nonFlag = 2
#                if writeCount == -1:
            writeCount = 5


#        if nonIndelPos > 0 and indelLen > 0:
#            nonIndelPos = pileupcolumn.pos
#            nonIndelFlag = True



        arrayt.append(y_data)
        idList.append(arrayt)

        if total >= 5:
            b_bq = np.mean(forward_baseQ[1:])
            idList[0].insert(-1, b_bq)
#            if nonFlag == 2 or (nonFlag == 0 and indelLen > 0):
            if writeCount == 0:
                wrcsv.writerow(idList[0])
                writeCount = -1


            ### initialize
            del forward_baseQ[0]
            forward_baseQ.append(0)
            del idList[0]


        if writeCount > 0:
            writeCount-=1

        # Initialize
        total += 1  # for count all the columns
        ins = 0  # insert count in a column
        dels = 0  # delete count in a column
        mapQual = [0, 0, 0]  # save sum of mapping quality separately ins/dels/none in each column.
        indelLen = 0  # if indel is found, save indel length
        y_data = 2  # defalut y_data = 2 (0:ins, 1:dels, 2:none)
        f_bq = 42
        b_bq = 42
        noCnt = 0

    ### last 5 element's backward base quality adding ###
    for i in range(4):
        b_bq = np.mean(forward_baseQ[(i + 1):-1])  # the last element is 0, so have to be excluded
        idList[i].insert(-1, b_bq)

        if writeCount==0:
            wrcsv.writerow(idList[i])
            writeCount = -1
        if writeCount > 0:
            writeCount -= 1
#        if (y_data !=2 and ((last_pileupcolumn+i+1) not in ydataList)) \
#                or (y_data==2 and indelLen > 0):
#        if (idList[i][-1] !=2 and ((last_pileupcolumn+i+1) not in ydataList)) \
#                or (idList[i][-1]==2 and idList[i][4] > 0):
#            wrcsv.writerow(idList[i])

    ### the last column's backward base quality adding
    idList[-1].insert(-1, 42.)
    if writeCount == 0:
        wrcsv.writerow(idList[-1])
        writeCount = -1

#    if (idList[-1][-1] != 2 and ((last_pileupcolumn + i + 1) not in ydataList)) \
#            or (idList[-1][-1] == 2 and idList[-1][4] > 0):
#        wrcsv.writerow(idList[-1])
#    malCnt = 0


# print(mal_list)
# print('mal_size: ',len(mal_list))




for non_pos in nonList:

    """
    frontLogFlag = -1
    backLogFlag = -1

    for i in range(nonIndexDict[non_pos] - 1, -1, -1):
        if (non_pos - int(log.iloc[i][1])) >= 50:
            break
        else:
            frontLogFlag = i

    for i in range(nonIndexDict[non_pos], len(log.index)):
        if (int(log.iloc[i][1]) - non_pos) >= 50:
            break
        else:
            backLogFlag = i

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

    nonIndelPos = 0
    nonIndelFlag = False

    last_pileupcolumn = 0

    ### if you want to set specific region, set start='startpos', stop='stoppos'
    for pileupcolumn in infile.pileup("chr21", non_pos -20, non_pos +20):

        lengthDict = dict()
        arrayt = []


        if pileupcolumn.pos < non_pos - 20:
            continue
        elif pileupcolumn.pos > non_pos + 20:
            break


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


        # add detected indels but not manipulated by us
        if nonIndelPos == 0 and indelLen > 0:
            nonIndelPos = pileupcolumn.pos
            nonIndelFlag = True


        ### the first column has the highest value(42) of forward_baseQualities
        ### and the last column has the highest value(42) of backward_baseQualities

        arrayt = [(float(ins) / pileupcolumn.n), (float(dels) / pileupcolumn.n),
                  (float((pileupcolumn.n - (ins + dels))) / pileupcolumn.n), pileupcolumn.n, indelLen]

        ### add mapping quality in list
        if (ins > 0):
            arrayt.append(float(mapQual[0] / ins))
        elif (dels > 0):
            arrayt.append(float(mapQual[1] / dels))
        else:
            arrayt.append(float(mapQual[2] / pileupcolumn.n))

        ### forward base quality
        if total == 0:
            arrayt.append(42.)
        else:
            arrayt.append(f_bq)

        ### if deletion, skip (indel_length-1) length backward.
        if arrayt[1] >= 0.2:
            delflag = indelLen

        # append y_data here
        # append y_data here
        if log.iloc[logIndex][2]==0 and arrayt[0]>=0.2 and float(log.iloc[logIndex][3]) == arrayt[4]:
            y_data = 0
        elif log.iloc[logIndex][2]==1 and arrayt[1]>=0.2 and float(log.iloc[logIndex][3]) == arrayt[4]:
            y_data = 1
        else:
            y_data = 2


        arrayt.append(y_data)
        idList.append(arrayt)

        if total >= 5:
            b_bq = np.mean(forward_baseQ[1:])
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
        f_bq = 42
        b_bq = 42
        noCnt = 0

    """
    for i in range(4):
        b_bq = np.mean(forward_baseQ[(i + 1):-1])  # the last element is 0, so have to be excluded
        idList[i].insert(-1, b_bq)

        if nonIndelFlag and nonIndelPos+5 == last_pileupcolumn+i+1:
            wrcsv.writerow(idList[i])
            nonIndelFlag = False
#        if (y_data !=2 and ((last_pileupcolumn+i+1) not in ydataList)) \
#                or (y_data==2 and indelLen > 0):
#        if (idList[i][-1] !=2 and ((last_pileupcolumn+i+1) not in ydataList)) \
#                or (idList[i][-1]==2 and idList[i][4] > 0):
#            wrcsv.writerow(idList[i])

    ### the last column's backward base quality adding
    idList[-1].insert(-1, 42.)
    if nonIndelFlag and nonIndelPos + 5 == last_pileupcolumn + 5:
#    if writeCount == 0:
        wrcsv.writerow(idList[-1])
        nonIndelFlag = False
    """


infile.close()
csvfile.close()
