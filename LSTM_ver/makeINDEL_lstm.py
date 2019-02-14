
# coding: utf-8

# In[348]:

import sys
import random
import numpy as np
from time import strftime

path = sys.argv[1]
name = sys.argv[2]
aname = path + "/hg19_chr21_" + name + ".fa"
yname = path + "/ylLog_" + name + "_nucleo.txt"
lname = path + "/ylLog_" + name + ".txt"

## define variables

# set indel's num and length
numofINDEL = 40000
lenofINDEL = 10

# set gap between indels
gap = 10

## fa file's all sequences
nList = []

## nList's length(48,129,895)
total = 0
modi = 0

total_len = 0

## first/last postion that start/end with "not N"
start_n, end_n = (0,0)

## define info(nuc_list, pos_list, len_list, type_list(0/1))

nuc_list = []
pos_list = []
pos_arr = []
len_list = []
type_list = []

# In[349]:


def INSERT(iLen, pos):
    
    global nList, nuc_list
    
    nucleo = ['A', 'G', 'C', 'T']
    
    tempChar = ''
    tempString = ""

    for i in range(iLen):
        tempChar = random.choice(nucleo)
        tempString += tempChar
        nList.insert(pos+i, tempChar)
        
    nuc_list.append(tempString)


# In[350]:


def DELETE(iLen, pos):
    
    global nList, nuc_list
    
    tempString = ""
    for i in range(iLen):
        if nList[pos] == 'N':
            iLen = i
            break
        tempString += nList.pop(pos)
    
    nuc_list.append(tempString)
    len_list[-1] = iLen


# In[352]:


## set nList from fa file
def set_nList():
    
    global nList, start_n, end_n, total
    
    pname = "./hg19_chr21.fa"
    with open(pname) as f:
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


# In[ ]:


def init_position():
   
    global pos_list, pos_arr 
    
    # 50 if for nList[iPos]=='N'
    candi_pos = [0]*len(nList)
    while len(pos_list) != numofINDEL+50: 
        #iPos = random.randrange(start_n+20, end_n-20)
        iPos = random.randrange(start_n, end_n-2)
        if candi_pos[iPos] or nList[iPos] == 'N': continue
        pos_list.append(iPos)
        candi_pos[iPos-gap-lenofINDEL:iPos+gap+lenofINDEL+1] = [1]*(2*gap+2*lenofINDEL+1)

    pos_list.sort()
    pos_arr = np.array(pos_list)

    for i in range(len(pos_list)-1):
        if pos_list[i+1]-pos_list[i] < gap:
            print("=== gap initial error ===")
            print(pos_list[i:i+2])

    pos_list = []

# In[355]:


if __name__ == "__main__":
    
    ## insert : delete = 1 : 1
    numofINS = numofINDEL/2
    numofDEL = numofINS
    if numofINDEL%2 != 0:
        numofINS = numofINDEL//2 + 1
        numofDEL = numofINDEL//2
    
    set_nList()
    print("set_nList : " + strftime("%y%m%d-%H%M%S"))
    
    init_position()
    print("set pos_list : " + strftime("%y%m%d-%H%M%S"))
        
    ## nIns is num of insertion / nDel is num of deletion
    nIns = 0
    nDel = 0

    while len(pos_list) != numofINDEL:
        
        # if modified position is 'N'
        if nList[pos_arr[0]] == 'N':

            print("N!!! continue.." + str(pos_arr[0]))
            pos_arr = pos_arr[1:]
            continue

        pos_list.append(pos_arr[0])
        pos_arr = pos_arr[1:]        

        # distance before insert/delete
        dis = pos_arr[0]-pos_list[-1]
        

       # choose type randomly
        while True:
            iType = random.randrange(0, 2)

            if iType==0 and nIns != numofINS: break
            elif iType==1 and nDel != numofDEL: break

        if dis <= gap:
            print("dis < gap !!error!!")
            print(pos_list[-1])
            print(new_value)
            print(pos_arr[0])
        
        iLen = random.randrange(0, lenofINDEL)+1
        
        type_list.append(iType)
        len_list.append(iLen)
        
        # if insert --> tmp is -1*iLen
        # if delete --> tmp is 1*iLen
        tmp = iLen * pow(-1, iType+1)
        total_len -= tmp
    
        if iType:
            nDel += 1
            DELETE(len_list[-1], pos_list[-1])
        else:
            nIns += 1
            INSERT(len_list[-1], pos_list[-1]) 

        pos_list[-1] += modi
        modi += tmp

    #aname = "./data/hg19_chr21_" + name + ".fa"
    with open(aname, 'w') as mseq:
        mseq.write(">chr21\n")
        nList.reverse()
        while nList:
            for i in range(50):
                if nList:
                    mseq.write(nList.pop())
                else:
                    break
                    
            mseq.write("\n")

    #yname = "./data/ylLog_" + name + "_nucleo.txt"
    with open(yname, 'w') as ylog:
        ylog.write("#	position	type	len	nucleo\n")
        for i in range(numofINDEL):
            ylog.write(str(i+1)+"\t"+str(pos_list[i])+"\t"+str(type_list[i])+"\t"+str(len_list[i])+"\t"+nuc_list[i]+"\n")


    #lname = "./data/ylLog_" + name + ".txt"
    with open(lname, 'w') as log: 
        log.write("#	num	position	type	len\n")
        for i in range(numofINDEL):
            log.write("%d	%d	%d	%d\n" %(i, pos_list[i]+1, type_list[i], len_list[i]))


