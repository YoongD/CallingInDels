import random

## position, indel(0/1), length, nucleo
posArr = []
indelArr = []
lenArr = []
#nucArr = []



def indexModifier(start, cnt):
	tmp = start;
	for i in range(cnt):
		if tmp > posArr[i]:
			if indelArr[i] == 0:	#indel
				tmp -= lenArr[i]
			elif indelArr[i] == 1:
				tmp += lenArr[i]
	return tmp


##########################################################################################
######################################### INSERT #########################################
##########################################################################################

def INSERT(nList, cnt, start, iLen) :
	nucleo = ['A', 'G', 'C', 'T']

	if nList[start-1] == 'N' and nList[start] == 'N' :
		return 0

	##### yoong modified #####
	tmp = indexModifier(start, cnt-1)
	
	## insert in array ##
	posArr.append(tmp)
	lenArr.append(iLen)
	indelArr.append(0)

#	log.write("%d	%d	" %(cnt, tmp))

	# insert random nucleotide ( length : iLen )
#	tempString = ""
#	tempChar = ''
	for j in range(iLen) :
#		tempChar = random.choice(nucleo)
#		tempString += tempChar
#		nList.insert(start+j, tempChar)
		nList.insert(start+j, random.choice(nucleo))
#	nucArr.append(tempString)
	
#	log.write("0	%d\n" %iLen)
	return iLen

##########################################################################################
######################################### DELETE #########################################
##########################################################################################

def DELETE(nList, cnt, start, iLen) :
	if nList[start] == 'N' :
		return 0


        ##### yoong modified #####
	tmp=indexModifier(start, cnt-1)

        ## insert in array ##
	posArr.append(tmp)
	indelArr.append(1)	#delete

#	log.write("%d   %d      " %(cnt, tmp))
	
#	tempString = ""
	for j in range(iLen) :
		if nList[start] == 'N' : 
			iLen = j			
			break	        	
#		tempString += nList.pop(start)
 		nList.pop(start)   
#	nucArr.append(tempString)
	lenArr.append(iLen)
	
#	log.write("1	%d\n" %iLen)
	return iLen


##########################################################################################
########################################## MAIN ##########################################
##########################################################################################

# 962597*50 + 45 = 48129895
pname = "../../hg19_chr21.fa"
aname = "hg19_chr21_20000.fa"
#lname = "logm20000_.txt"

# fa file to list
with open(pname) as f:
	f.readline()
	nList = list(f.read().upper().replace('\n', ''))

# open modified .fa file
mseq = open(aname, 'w')

# open text for log
#log = open(lname, 'w') 

# initialize variables
total = 48129895

cnt = 0
nIns = 0
nDel = 0
mod = 0

numofINDEL = 20000
lenofINDEL = 10

sList = []
while len(sList) < numofINDEL:
	s = random.randrange(0, total)
	
	if nList[s] == 'N':
        	if s == 0 or nList[s-1] == 'N':
        	    continue
        	if s+1 >= total or nList[s+1] == 'N':
        	    continue
	
	if s-lenofINDEL in nList or s+lenofINDEL in nList:
		continue
    	
	sList.append(s)
#sList.sort()

numofINS = numofINDEL/2
numofDEL = numofINDEL/2
if(numofINDEL%2 != 0):
	numofINS = numofINDEL//2 + 1
	numofDEL = numofINDEL//2
print("hi/..?\n")
# total num of INDEL : cnt
# in this case we make 100 indels
# nIns is num of insertion / nDel is num of deletion
# in this case we make 50 insertions and 50 deleteion
while cnt < numofINDEL :
	#start = sList[cnt]+mod 
	iLen = random.randrange(0, lenofINDEL)+1       
 	
	Type = random.randrange(0, 2)
	if Type == 0 and nIns<numofINS :	
		tmp=INSERT(nList, cnt+1, sList[cnt], iLen)
		if tmp > 0:	
			nIns+=1
			total+=tmp
			#mod+=tmp
	
	elif nDel < numofDEL :
		tmp=DELETE(nList, cnt+1, sList[cnt], iLen)
		if tmp > 0:
			nDel+=1
			total-=tmp
			#mod-=tmp
	cnt = nIns+nDel
    
# print result
# compare var total with length of list
print("TOTAL : ")
print(total)
print(len(nList))

# [hg19_chr21.fa 's modified version]
mseq.write(">chr21\n")
nList.reverse()
while nList :
	for i in range(50) :
		if nList :
			mseq.write(nList.pop())
		else :
			break
	mseq.write("\n")



#### write yoong log ####
### yoong log ###
# open text for log
yname = "yLog20000.txt"
ylog = open(yname, 'w')

ylog.write("#	position	indel(0/1)	len\n")
for i in range(numofINDEL):
	ylog.write(str(i+1)+"\t"+str(posArr[i])+"\t"+str(indelArr[i])+"\t"+str(lenArr[i])+"\n")

ylog.close()
mseq.close()
#log.close()


