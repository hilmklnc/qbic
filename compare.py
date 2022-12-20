#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 16:14:27 2022

@author: hilmi
"""



# ZNF24 --> Mus_musculus|NA|Unpublished|Zfp24.txt 

# Homo_sapiens|NA|Shen2018|E2F1.txt --> bak 
# python3 olskmer.py /home/hilmi/qbic/PBM_inputs/Homo_sapiens_NA_Shen2018_E2F1.txt /home/hilmi/qbic/ -k 6 -d 128 -g 0 -p 0

# Homo_sapiens|NA|Shen2018|MYC-MAX.txt --> myc
# Homo_sapiens|NA|Unpublished|GATA1.txt --> gata1

import pandas as pd

def compare_results(control_data, predict_data):
    control_data.columns = ["diff","t"]
    df_bool = pd.DataFrame(predicted_data.dropna() == control_data.dropna())
    print("Identical Shapes:",control_data.shape == predicted_data.shape) 
    print("\nIdentical Descriptive Statistics:\n",control_data.describe() == predicted_data.describe())
    print("\nAll rows and columns are identical:\n",df_bool.all())

# H.sapiens, E2F1 TF: 
control_data = pd.read_csv("/home/hilmi/qbic/control_predicteds/prediction6mer.Homo_sapiens_NA_Shen2018_E2F1.txt", sep=" ")
predicted_data = pd.read_csv("/home/hilmi/qbic/my_predictions/prediction6mer.Homo_sapiens_NA_Shen2018_E2F1.txt", sep= " ")

compare_results(control_data, predicted_data)

p_control = pd.read_csv("/home/hilmi/qbic/control_predicteds/pval6mer.Homo_sapiens_NA_Shen2018_E2F1.csv",header=None)    
p_predict = pd.read_csv("/home/hilmi/qbic/control_predicteds/pval6mer.Homo_sapiens_NA_Shen2018_E2F1.csv",header=None)

df_p = pd.DataFrame(p_control.dropna()==p_predict.dropna())
print(df_p.all())


# H.sapiens, MYC TF:
control_data = pd.read_csv("/home/hilmi/qbic/control_predicteds/prediction6mer.Homo_sapiens_NA_Shen2018_MYC-MAX.txt", sep=" ")
predicted_data = pd.read_csv("/home/hilmi/qbic/my_predictions/prediction6mer.Homo_sapiens_NA_Shen2018_MYC-MAX.txt", sep= " ")

compare_results(control_data, predicted_data)

# D.melanogaster datasets
control_data = pd.read_csv("/home/hilmi/qbic/control_predicteds/prediction6mer.Drosophila_melanogaster_M01487_1.94d_Zoo_01_3110.txt", sep=" ")
predicted_data = pd.read_csv("/home/hilmi/qbic/my_predictions/prediction6mer.Drosophila_melanogaster_M01487_1.94d_Zoo_01_3110.txt", sep= " ")

compare_results(control_data, predicted_data)

