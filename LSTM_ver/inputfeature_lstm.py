import numpy as np
import csv
import pysam
import pandas as pd
import random
import sys

from time import strftime

#### set 'path', 'name'
path = sys.argv[1]
name = sys.argv[2]
indel_num = int(sys.argv[3])

#### set file name
f_name = path + "/hg19_chr21_" + name + ".fa"
i_name = path + "/realigned_" + name + ".bam"
l_name = path + "/ylLog_" + name + "_nucleo.txt"
csv_name0 = path + "/" + name + "_features_lstm_0.csv"
csv_name1 = path + "/" + name + "_features_lstm_1.csv"
csv_name2 = path + "/" + name + "_features_lstm_2.csv"

### max value for normalize
indel_max_length = 30.0
coverage_max = 100.0
mappingQ_max = 70.0
baseQ_max = 42.0

### gap for finding mis-aligned indel
indel_range = 60

### csv writer
csvfile0 = open(csv_name0, 'w')
csvfile1 = open(csv_name1, 'w')
csvfile2 = open(csv_name2, 'w')

wrcsv0 = csv.writer(csvfile0)
wrcsv1 = csv.writer(csvfile1)
wrcsv2 = csv.writer(csvfile2)

### read log and bam file
log = pd.read_table(l_name)
infile = pysam.AlignmentFile(i_name, "rb")

print("start : " + strftime("%y%m%d-%H%M%S"))

## logIndex : 봐야할 첫 로그 인덱스
## backLogFlag : logIndex를 기점으로 indel_range를 넘지않는 범위까지 확장
## [frontLogFlag, backLogFlag] == range(frontLogFlag, backLogFlag+1)

### initializing for while-loop
logIndex = 0
exclude_list = []

while True:

    if logIndex == indel_num: break

    frontLogFlag = logIndex
    backLogFlag = logIndex
    for i in range(logIndex+1, indel_num):
        if (int(log.iloc[i][1]) - int(log.iloc[backLogFlag][1])) > indel_range+5:
            break
        else:
            backLogFlag = i

    if backLogFlag+1 != indel_num:
        logIndex = backLogFlag+1
    else:
        logIndex = indel_num

    exclude_list += list(range(log.iloc[frontLogFlag][1]-indel_range, log.iloc[backLogFlag][1]+indel_range))

    idList = []
    ### if you want to set specific region, set start='startpos', stop='stoppos'
    for pileupcolumn in infile.pileup("chr21", int(log.iloc[frontLogFlag][1])-(indel_range+5), int(log.iloc[backLogFlag][1])+(indel_range+5)+1, truncate=True):

        ### Initialize
        ins = 0  # insert count in a column
        dels = 0  # delete count in a column
        noCnt = 0

        mapQual = [0, 0, 0]  # save sum of mapping quality separately ins/dels/none in each column.
        baseQ = 0.0

        lengthDict = dict()
        currPos = pileupcolumn.pos
        arrayt = []
        
        coverage = len(pileupcolumn.pileups)
        if coverage == 0:
            print("coverage is 0!! : ", currPos)

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
 
                if pileupread.indel in lengthDict.keys():
                    lengthDict[pileupread.indel] += 1
                else:
                    lengthDict[pileupread.indel] = 1

            ### not indel
            elif pileupread.indel == 0:
                mapQual[2] += pileupread.alignment.mapping_quality
                noCnt += 1

            ### query_position is None if is_del or is_refskip is set
            if pileupread.is_del or pileupread.is_refskip:
                baseQ += 0
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
        if (arrayt[0]>=0.2) or (arrayt[1]>=0.2):

            for currLogFlag in range(frontLogFlag, backLogFlag+1):

                if log.iloc[currLogFlag][1]-indel_range <= currPos <= log.iloc[currLogFlag][1]+indel_range:

                    if indelLen == (log.iloc[currLogFlag][3]) * pow(-1, log.iloc[currLogFlag][2]):
                        if log.iloc[currLogFlag][2]==0:
                            y_data = 0
                            break
                        elif log.iloc[currLogFlag][2]==1:
                            y_data = 1
                            break

        arrayt.append(y_data)
        idList.append(arrayt)

        if len(idList) == 11:

            tmpList = [idList[0], idList[10], idList[1], idList[9], idList[2], 
                       idList[8], idList[3], idList[7], idList[4], idList[6], idList[5]]

            if idList[5][-1] == 0:   wrcsv0.writerow(tmpList)
            elif idList[5][-1] == 1: wrcsv1.writerow(tmpList)
            else:                    wrcsv2.writerow(tmpList)

            del idList[0]

print("finish1 : " + strftime("%y%m%d-%H%M%S"))

print("NONE FEATURE START")
#######################################################################################
#                                                                                     #
#                                 none feature                                        #
#                                                                                     #
#######################################################################################

### calculate start_n and end_n
start_n, end_n = (0,0)

with open(f_name) as f:
    f.readline()
    nList = list(f.read().upper().replace('\n', ''))
        
nList_len = len(nList)
    
# start after N's sequence
start_n = 0
while True:
    if nList[start_n] != 'N': break
    else:
        start_n += 1
    
end_n = nList_len-1
while True:
    if nList[end_n] != 'N': break
    else:
        end_n -= 1

### pick position randomly
random_list = []
while len(random_list)<indel_num/2+20 :
    tmp = random.randrange(start_n+100, end_n-100)
    if (not tmp in random_list) and (not tmp in exclude_list) and (not 'N' in nList[tmp-10:tmp+10]):
        random_list.append(tmp)

print(len(random_list))
print("finish picking random position : " + strftime("%y%m%d-%H%M%S"))

for randPos in random_list:

    idList = []

    ### if you want to set specific region, set start='startpos', stop='stoppos'
    for pileupcolumn in infile.pileup("chr21", randPos-5, randPos+6, truncate=True):

        #cnt += 1

        ### Initialize
        ins = 0  # insert count in a column
        dels = 0  # delete count in a column
        noCnt = 0

        mapQual = [0, 0, 0]  # save sum of mapping quality separately ins/dels/none in each column.
        baseQ = 0.0

        lengthDict = dict()
        currPos = pileupcolumn.pos
        
        coverage = len(pileupcolumn.pileups)
        if coverage == 0:
            print("coverage is 0!! : ", currPos)
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
 
                if pileupread.indel in lengthDict.keys():
                    lengthDict[pileupread.indel] += 1
                else:
                    lengthDict[pileupread.indel] = 1

            ### not indel
            elif pileupread.indel == 0:
                mapQual[2] += pileupread.alignment.mapping_quality
                noCnt += 1

            ### query_position is None if is_del or is_refskip is set
            if pileupread.is_del or pileupread.is_refskip:
                baseQ += 0
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
                         del_mq/mappingQ_max, none_mq/mappingQ_max, baseQ/baseQ_max, 2 ]

        idList.append(arrayt)
        if len(idList) == 11:
            tmpList = [idList[0], idList[10], idList[1], idList[9], idList[2], 
                       idList[8], idList[3], idList[7], idList[4], idList[6], idList[5]]
            wrcsv2.writerow(tmpList)

infile.close()
csvfile0.close()
csvfile1.close()
csvfile2.close()

print("finish2 : " + strftime("%y%m%d-%H%M%S"))
