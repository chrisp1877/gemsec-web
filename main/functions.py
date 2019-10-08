import main.src.MD_Analysis as md
import pandas as pd
from sklearn.decomposition import PCA
import numpy as np


''' PARAMETERS '''
''' analysis = 'double_hilbert', 'hilbert', 'cos_sin_pca' '''
''' visualization = 'prob_map', 'pca_gif', 'density plot' '''
def md_vis(pdb_filepath, analysis, visualization):
   filepath = pdb_filepath
   with open(filepath) as file:
      phi_psi = md.Angle_Calc.get_phi_psi(filepath)
      sin_cos_df = md.Angle_Calc.get_sin_cos(phi_psi)
      df = sin_cos_df
      phi_psi = md.Angle_Calc.get_phi_psi(filepath)
      sin_cos_df = md.Angle_Calc.get_sin_cos(phi_psi)
      sin_cos_df.reset_index(drop=True, inplace=True)    
      if analysis == "double_hilbert":
         df = md.Hilbert_Transform.hilb_collapse(df) ##is a dataframe
      elif analysis == "hilbert":
         Hilbert_20d = md.Hilbert_Transform.hilb_collapse(df) ##is a dataframe
         df = md.Hilbert_Transform.hilb_collapse(Hilbert_20d) ##is a dataframe

      pca = md.Dim_Reduction.pca(df, filepath, '.', "2d") ##is a tuple
      prep_var = np.round(pca[1].explained_variance_ratio_ * 100, decimals = 1)

      if visualization == "prob_map":
         plt_path = md.PCA_Components.PC_prob_map(pca[0], per_var=prep_var, wd= 'media/')
      elif visualization == "pca_gif":
         plt_path = md.PCA_Components.gen_2d_PCA_gif(pca[0], per_var=prep_var, wd= 'media/')
      else:
         plt_path = md.PCA_Components.PC_den_plt(pca[0], 'media/')

      return "/" + plt_path



