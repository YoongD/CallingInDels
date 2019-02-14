import numpy as np

import pandas as pd

y_path = './input_split_file/'
df = pd.read_csv('./111_inputdata_rand.csv', names=['ins_ratio', 'del_ratio', 'normal_ratio',
                                              'coverage', 'length', 'mQual', 'f_bQual', 'b_bQual', 'y_data'])

df = df.sort_values(by=['y_data'], axis=0)

y_0 = df[df.y_data == 0]
y_1 = df[df.y_data == 1]
y_2 = df[df.y_data == 2]

y_0.to_csv(y_path+'y_0.csv', header=False, index=False)
y_1.to_csv(y_path+'y_1.csv', header=False, index=False)
y_2.to_csv(y_path+'y_2.csv', header=False, index=False)
