#!/usr/bin/env python
# -*-coding:utf-8 -*-
'''
@File    :   main.py
@Time    :   2023/04/09 10:49:57
@Author  :   Bio-Maw 
@Version :   1.0
@Contact :   martint.maw@gmail.com
@Desc    :   None
'''

import pandas as pd
from pathlib import Path

INPUT_DIR = Path('/mnt/c/Projects/Data/ev_sizes')
file_list= list(INPUT_DIR.glob('*.txt'))
