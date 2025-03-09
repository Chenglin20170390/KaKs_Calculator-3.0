#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@File    :   22.1.kaks.matrix.py
@Time    :   2025/03/09 23:17:49
@Author  :   Lin Cheng 
@Version :   1.0
@Contact :   chenglin_solab@163.com
@License :   (C)Copyright 2022-2023, CAAS ShenZhen
@Desc    :   None
'''

import datetime, sys
import pandas as pd
import numpy as np

# here put the import lib
print('***********************************************************')
print('Start Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')

infile=sys.argv[1]
outfile=sys.argv[2]
f_in=open(infile,'r')
f_out=open(outfile,'w')

# 定义包含文本数据的多行字符串（你也可以从文件中读取）
# 将每行数据解析为字典列表
# 用于记录行名和列名的出现顺序

# 记录行名、列名，并保持输入顺序
rows, cols = [], []
data_dict = {}

for line in f_in:
    parts = line.split()
    key = parts[0]
    value = float(parts[4])
    row_name, col_name = key.split('&', 1)  # 按 '&' 拆分

    # 维护行、列顺序
    if row_name not in rows:
        rows.append(row_name)
    if col_name not in cols:
        cols.append(col_name)

    # 存储数据
    data_dict[(row_name, col_name)] = value

# 创建 DataFrame，初始化所有值为 0
df = pd.DataFrame(0, index=rows, columns=cols, dtype=float)

# 填充数据
for (r, c), v in data_dict.items():
    df.loc[r, c] = v  # 确保用 .loc 进行赋值

# 将矩阵写入 TXT 文件（使用制表符分隔，并保留行索引）
df.to_csv(f_out, sep="\t")
print(f"矩阵已写入文件：{outfile}")

f_in.close()
f_out.close()
print('***********************************************************')
print('Stop Time:	'+ str(datetime.datetime.now()))
print('***********************************************************')



