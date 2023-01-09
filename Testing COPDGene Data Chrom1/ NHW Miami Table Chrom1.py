import numpy as np
import random
import scipy as sp
import pandas as pd
import seaborn as sns
import sklearn.linear_model
import matplotlib.pyplot as plt
from statsmodels.api import OLS
import statsmodels.api as sm
from sklearn.linear_model import LinearRegression
from scipy import stats
from scipy.stats import gaussian_kde as GKDE
from scipy.stats import uniform, norm, beta,skewnorm,gamma
from mpl_toolkits import mplot3d
import scipy.io as sio
from statsmodels.distributions.empirical_distribution import ECDF
from joblib import dump, load
from astropy.table import Table, Column
import scikit_posthocs as sp
from os.path import join
from pandas_plink import read_plink
from pandas_plink import get_data_folder




(bim, fam, bed)=read_plink("/home/math/jananin/UQ/NHWResults/CG10k_NHW_hg19_Oct2017",verbose=True)

NHWdf=pd.read_csv('/home/math/jananin/UQ/NHWResults/NHW1/Manipulate/NHWdf.csv')

NHWdf=NHWdf.iloc[:,1:]

KruskalPvalueA=np.loadtxt('/home/math/jananin/UQ/NHWResults/NHW1/Manipulate/KruskalPvalueA.txt')
KruskalPvalueE=np.loadtxt('/home/math/jananin/UQ/NHWResults/NHW1/Manipulate/KruskalPvalueE.txt')

OLS_PvaluesA=np.loadtxt('/home/math/jananin/UQ/NHWResults/NHW1/Manipulate/OLS_PvaluesA.txt')
OLS_PvaluesE=np.loadtxt('/home/math/jananin/UQ/NHWResults/NHW1/Manipulate/OLS_PvaluesE.txt')


Ch=bim[bim['chrom'] == '1']


pos=Ch.loc[Ch['snp'].isin(NHWdf.columns), 'pos']

np.savetxt("pos.txt",pos,fmt='%+5.3d')

OLSManA = pd.DataFrame(columns=['CHR','BP','P','SNP'])
OLSManA.BP=pos
OLSManA.CHR=1
OLSManA.P=OLS_PvaluesA
OLSManA.SNP=NHWdf.columns

OLSManA.to_csv('OLSManA.csv')


KManA = pd.DataFrame(columns=['CHR','BP','P','SNP'])
KManA.BP=pos
KManA.CHR=1
KManA.P=KruskalPvalueA
KManA.SNP=NHWdf.columns

KManA.to_csv('KManA.csv')


OLSManE = pd.DataFrame(columns=['CHR','BP','P','SNP'])
OLSManE.BP=pos
OLSManE.CHR=1
OLSManE.P=OLS_PvaluesE
OLSManE.SNP=NHWdf.columns

OLSManE.to_csv('OLSManE.csv')


KManE = pd.DataFrame(columns=['CHR','BP','P','SNP'])
KManE.BP=pos
KManE.CHR=1
KManE.P=KruskalPvalueE
KManE.SNP=NHWdf.columns

KManE.to_csv('KManE.csv')


