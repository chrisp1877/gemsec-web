# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 12:42:53 2019

@author: knity
"""

import pandas as pd
from sklearn.decomposition import PCA
import numpy as np
from src.MD_Analysis import *

file = eg.fileopenbox(msg = "Select the pdb to be analyzed")

save_dir = eg.diropenbox(msg = "Select the directory to save to")

phi_psi = md.Angle_Calc.get_phi_psi(flow_variables['file'])
sin_cos_df = md.Angle_Calc.get_sin_cos(phi_psi)

sin_cos_df.reset_index(drop=True, inplace=True)      

##Choose an analysis
msg ="Choose an Analysis to run"
title = "MD Analysis"
choices = ["Cos/Sin PCA", "Hilbert PCA", "Double Hilbert PCA"]
choice = choicebox(msg, title, choices)

df = sin_cos_df

if choice == "Double Hilbert PCA":
	df = Hilbert_Transform.hilb_collapse(df) ##is a dataframe
elif choice == "Hilbert PCA":
	Hilbert_20d = Hilbert_Transform.hilb_collapse(df) ##is a dataframe
	df = Hilbert_Transform.hilb_collapse(Hilbert_20d) ##is a dataframe
else:
	df

pca = Dim_Reduction.pca(df, file, save_dir, "2d") ##is a tuple

prep_var = np.round(pca[1].explained_variance_ratio_ * 100, decimals = 1)

##Choose a plot
msg ="Choose a Plot"
title = "MD Analysis"
choices = ["Probability Map", "Density Plot", "PCA Gif"]
choice = choicebox(msg, title, choices)

if choice == "Probability Map":
   md.PCA_Components.PC_prob_map(input_table_1, prep_var, save_dir)
elif choice == "PCA Gif":
   md.PCA_Components.gen_2d_PCA_gif(input_table_1, prep_var, save_dir)
else: #plot == "Density Plot":
   md.PCA_Components.PC_den_plt(input_table_1, save_dir)

