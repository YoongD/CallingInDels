import csv
import random
import pandas as pd


"""
 read all the position of indels from the log file and compare neighbor
 The maximum neigborhood size is [-45, 45]
"""

#log = pd.read_table('yLog20000.txt')

df = pd.read_csv('./newtrain/notrandom_newtrain_inputdata20000.csv', names=['ins_ratio', 'del_ratio', 'normal_ratio',
                                              'coverage', 'length', 'mQual', 'f_bQual', 'b_bQual', 'y_data'])


# => delete in the next files
df['length'] *= 10


### normalize ###
# addition : length#
indelLen = df['length']
max_len = max(indelLen)
indelLen /= max_len
df['length'] = indelLen



cover_column = df['coverage']
mq_column = df['mQual']

# get max coverage and mapping quaility
max_cover = max(cover_column)
max_mq = max(mq_column)

print('max_len', max_len)
print('max_coverage', max_cover)
print('max_mappingQ',max_mq)

# divide by max value
cover_column /= max_cover
mq_column /= max_mq

# edit df by normalized column
df['coverage'] = cover_column
df['mQual'] = mq_column

# base Quality normalized
df['f_bQual'] /= 42
df['b_bQual'] /= 42


"""
 Remove the elements which have any duplicated feature
 and also remove the columns that don't need for the features
 Then, save the features of in/del/none in the same portion 
"""
#df = df.drop(['total', 'position'], 1)
df = df.drop_duplicates()  # remove all the duplicants in df
print(df[df.normal_ratio < 0.2])

### add 50, 50, 50
f = open('./newtrain/notrandom_trainingData_20000.csv', 'w')
wr = csv.writer(f)

indelCount = 0
countNone = 0
for i in range(len(df.index)):
    if df.iloc[i, 2] == 1.0:
        continue

    elif df.iloc[i, 0] > 0.2 or df.iloc[i, 1] > 0.2:
        wr.writerow(df.iloc[i, :])
        indelCount += 1

    elif countNone < 20000:
        wr.writerow(df.iloc[i, :])
        countNone += 1

print("INDELcnt : " + str(indelCount))
print("NONEcnt : " + str(countNone))

""" indel 20000, none 20000"""
#total = indelCount/2 - countNone
while countNone < 20000:  #mod#
    randPos = random.randrange(0, len(df.index))
    if df.iloc[randPos, 2] == 1.0:
        wr.writerow(df.iloc[randPos, :])
        countNone += 1  #mod# increase iteration number

f.close()
