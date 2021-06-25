#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 20:32:42 2021

@author: givanna
"""

# Split data into test and train based on LOOCV
import pandas as pd
from sklearn.model_selection import LeaveOneOut

data = pd.read_csv("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/patient_groups.csv")
patients = data['PatientName']
loocv = LeaveOneOut()
data_splits = loocv.split(patients)

train_dats = []
test_dats = []

for train_idx, test_idx in data_splits:
    train_dat, test_dat = patients[train_idx], patients[test_idx]
    train_dats.append(train_dat.to_numpy())
    test_dats.append(test_dat.to_numpy()[0])
    
fold_idx = 1
for train_dat, test_dat in zip(train_dats, test_dats):
    data_fold = pd.DataFrame(train_dat, columns=['PatientName'])
    data_fold['set'] = 'train'

    data_test = pd.DataFrame({'PatientName': [test_dat], 'set': ['test']})
    data_fold = data_fold.append(data_test)
    data_fold.to_csv("~/Documents/phd/tracksom_differential/cmv_complete_tp_only/fold_" + str(fold_idx) + ".csv", index=False)
    fold_idx = fold_idx + 1
