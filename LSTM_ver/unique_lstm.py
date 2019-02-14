import pandas as pd
import numpy as np
import math
import sys

import csv
from ast import literal_eval

path = sys.argv[1]
name = sys.argv[2]
indel_num = int(sys.argv[3])

ins_name = path + '/' + name + "_features_lstm_0.csv"
dels_name = path + '/' + name + "_features_lstm_1.csv"
none_name = path + '/' + name + "_features_lstm_2.csv"

ins_df = pd.read_csv(ins_name, names=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'])
dels_df = pd.read_csv(dels_name, names=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'])
none_df = pd.read_csv(none_name, names=['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11'])

# remove all the duplicants in df
ins_df = ins_df.drop_duplicates()
dels_df = dels_df.drop_duplicates()
none_df = none_df.drop_duplicates()

indel_len = 0
if len(ins_df) >= len(dels_df):
    indel_len = len(dels_df)
elif len(ins_df) < len(dels_df):
    indel_len = len(ins_df)

if indel_len%2 != 0:
    indel_len = indel_len-1

ins_df = ins_df.sample(frac=1).reset_index(drop=True)[:indel_len]
dels_df = dels_df.sample(frac=1).reset_index(drop=True)[:indel_len]
none_df = none_df.sample(frac=1).reset_index(drop=True)[:indel_len]

train_size = int(indel_len/4*3)

train_df = pd.concat([ins_df[:train_size], dels_df[:train_size], none_df[:train_size]]).reset_index(drop=True)
val_df = pd.concat([ins_df[train_size:], dels_df[train_size:], none_df[train_size:]]).reset_index(drop=True)

train_df = train_df.sample(frac=1).reset_index(drop=True)
val_df = val_df.sample(frac=1).reset_index(drop=True)

train_df.to_csv(path + "/" + name + "_train_dataset.csv", header=False, index=False)
val_df.to_csv(path + "/" + name + "_val_dataset.csv", header=False, index=False)

