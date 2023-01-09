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



gType=bed.compute()
Ch=bim[bim['chrom'] == '1']
ID=fam.iloc[:,0:2]
SNP=Ch.iloc[:,0:2]

df = pd.DataFrame(gType, columns = ID.set_index('fid')['iid'])
df=df.iloc[SNP.index.values]
dataNHW=pd.read_csv('/home/math/jananin/UQ/NHWResults/NHWfactor.txt', sep="\t")


IID=dataNHW['IID'].values
NHWdf=df[df.columns[df.columns.isin(IID)]]

NHWdf.insert(loc=0, column="snp", value=SNP['snp'])

dataNHW= dataNHW[dataNHW.IID.isin(NHWdf.columns.values)]

dataNHW=dataNHW.reset_index(drop=True)

Y=dataNHW["factor1_scores"]
X=dataNHW[['age_enroll','gender','SmokCigNow','pc1', 'pc2', 'pc3', 'pc4','pc5']]

X = sm.add_constant(X)
est = sm.OLS(Y, X)
est2 = est.fit()

dataNHW["Emphysema-PredominantAxis"]=est2.resid

Y=dataNHW["factor2scoresflipped"]
X=dataNHW[['age_enroll','gender','SmokCigNow','pc1', 'pc2', 'pc3', 'pc4','pc5']]

X = sm.add_constant(X)
est = sm.OLS(Y, X)
est2 = est.fit()

dataNHW["Airway-Predominant Axis"]=est2.resid


dataNHW.to_csv('dataNHW.csv')



NHWdf=NHWdf.T
NHWdf.columns = NHWdf.iloc[0]
NHWdf= NHWdf[1:]
NHWdf.index.name = None
NHWdf=NHWdf.reset_index(drop=True)
#NHWd=NHWdf.iloc[: , 1:]

N=np.zeros((3,NHWdf.shape[1]))

for j in range(NHWdf.shape[1]):
    for i in range(NHWdf.shape[0]):
        if NHWdf.iloc[i,j]==0 :
            N[0,j]+= 1
        elif NHWdf.iloc[i,j]==1:
            N[1,j] += 1
        else:
            N[2,j]+= 1

dropN0=np.where( N[0,:] < .01*(NHWdf.shape[0]))
NHWdf.drop(NHWdf.columns[dropN0],axis=1,inplace=True)

N=np.zeros((3,NHWdf.shape[1]))

for j in range(NHWdf.shape[1]):
    for i in range(NHWdf.shape[0]):
        if NHWdf.iloc[i,j]==0 :
            N[0,j]+= 1
        elif NHWdf.iloc[i,j]==1:
            N[1,j] += 1
        else:
            N[2,j]+= 1

dropN1=np.where( N[1,:] < .01*(NHWdf.shape[0]))
NHWdf.drop(NHWdf.columns[dropN1],axis=1,inplace=True)


MAF=np.zeros(NHWdf.shape[1])

for j in range(NHWdf.shape[1]):
    MAF[j]=(2*N[0,j]+N[1,j])/(2*NHWdf.shape[0])

np.savetxt('MAF1.txt',MAF, fmt='%1.2f')
dropSNP=np.where( MAF < 0.05)
NHWdf.drop(NHWdf.columns[dropSNP],axis=1,inplace=True)


NHWdf.to_csv('NHWdf.csv')

MAF=np.zeros(NHWdf.shape[1])

for j in range(NHWdf.shape[1]):
    MAF[j]=(2*N[0,j]+N[1,j])/(2*NHWdf.shape[0])

np.savetxt('MAF2.txt',MAF,fmt='%1.2f')


nIter=NHWdf.shape[1]
nSamp=NHWdf.shape[0]

y=dataNHW["Emphysema-PredominantAxis"].values
Xs=np.zeros((nSamp,nIter))

for j in range(nIter):
    Xs[:,j]=NHWdf.iloc[:,j].values

N=np.zeros((3,nIter),dtype=int)

for j in range(nIter):
    for i in range(nSamp):
        if Xs[i,j]==0 :
            N[0,j]+= 1
        elif Xs[i,j]==1:
            N[1,j] += 1
        else:
            N[2,j]+= 1

y_obs= np.zeros((nSamp,nIter))


for j in range(nIter):
    index0=np.where(Xs[:,j] == 0)
    index1=np.where(Xs[:,j] == 1)
    index2=np.where(Xs[:,j] == 2)
    indexnan=np.argwhere(np.isnan(Xs[:,j]))
    Xs0=Xs[:,j][index0]
    Xs1=Xs[:,j][index1]
    Xs2=Xs[:,j][index2]
    Xsnan=Xs[:,j][indexnan]
    Xs[:,j]=np.concatenate((Xs0,Xs1,Xs2,Xsnan), axis=None)
    y_obs0=y[index0]
    y_obs1=y[index1]
    y_obs2=y[index2]
    y_obsnan=y[indexnan]
    y_obs[:,j]=np.concatenate((y_obs0,y_obs1,y_obs2,y_obsnan), axis=None)


np.savetxt(' y_obs.txt', y_obs)
np.savetxt('Xs.txt',Xs)


Xs_bar = np.zeros(nIter)
SEb = np.zeros(nIter)
SEe=np.zeros(nIter)
for j in range(nIter):
    Xs_bar[j]=np.mean(Xs[:,j])
    Sxx = np.sum(Xs[:,j]*Xs[:,j])-nSamp*Xs_bar[j]*Xs_bar[j]
    sigma=.9
    SEb[j]= np.sqrt(sigma/Sxx)
    SEe[j]=.9-SEb[j]

OLS_Pvalues=np.zeros(nIter)
OLS_Beta=np.zeros(nIter)
for j in range(nIter):
    X=sm.add_constant(Xs[:,j])
    est=sm.OLS(y_obs[:,j],X,missing='drop').fit()
    OLS_Pvalues[j]=est.pvalues[1]
    OLS_Beta[j]=est2.params[1]

np.savetxt('OLS_PvaluesE.txt',OLS_Pvalues,fmt='%1.4e')
np.savetxt('OLS_BetaE.txt',OLS_Beta)


y_obs_bar = np.zeros((3,nIter))
SEy=np.zeros((3,nIter))
Sy=np.zeros(3)

for j in range(nIter):
        y_obs_bar[0,j]=np.mean(y_obs[0:N[0,j],j])
        y_obs_bar[1,j]=np.mean(y_obs[N[0,j]:N[1,j]+N[0,j],j])
        y_obs_bar[2,j]=np.mean(y_obs[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j])
        Sy[0]= np.sum(y_obs[0:N[0,j],j]*y_obs[0:N[0,j],j])-N[0,j]*y_obs_bar[0,j]*y_obs_bar[0,j]
        Sy[1]= np.sum(y_obs[N[0,j]:N[0,j]+N[1,j],j]*y_obs[N[0,j]:N[1,j]+N[0,j],j])-N[1,j]*y_obs_bar[1,j]*y_obs_bar[1,j]
        Sy[2]= np.sum(y_obs[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j]*y_obs[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j])-N[2,j]*y_obs_bar[2,j]*y_obs_bar[2,j]
        SEy[0,j]= np.sqrt(Sy[0]/(N[0,j]-1))
        SEy[1,j]= np.sqrt(Sy[1]/(N[1,j]-1))
        SEy[2,j]= np.sqrt(Sy[2]/(N[2,j]-1))

Beta = np.zeros((nSamp,nIter))
eps = np.zeros((nSamp,nIter))
y_initi= np.zeros((nSamp,nIter))
seedinit1 = 1

for j in range(nIter):
    #seedinit1 = seedinit1+j
    #np.random.seed(seedinit1)
    Beta[0:N[0,j],j] = norm.rvs(loc=y_obs_bar[0,j] , scale=SEy[0,j], size=N[0,j])
    Beta[N[0,j]:N[1,j]+N[0,j],j] = norm.rvs(loc=y_obs_bar[1,j], scale=SEy[1,j], size=N[1,j])
    Beta[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j] = norm.rvs(loc=y_obs_bar[2,j] , scale=SEy[2,j], size=N[2,j])

    eps[0:N[0,j],j]=np.random.normal(0,.34,size=N[0,j])
    eps[N[0,j]:N[1,j]+N[0,j],j]=np.random.normal(0,.34,size=N[1,j])
    eps[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j]=np.random.normal(0,.34,size=N[2,j])

y_initi=Beta+eps


r = np.zeros((nSamp,nIter))
for j in range(nIter):
            kde_y_obs0=GKDE(y_obs[0:N[0,j],j])
            kde_y_initi0=GKDE(y_initi[0:N[0,j],j])
            r[0:N[0,j],j]=np.divide(kde_y_obs0(y_initi[0:N[0,j],j]),kde_y_initi0(y_initi[0:int(N[0,j]),j]))
            kde_y_obs1=GKDE(y_obs[N[0,j]:N[1,j]+N[0,j],j])
            kde_y_initi1=GKDE(y_initi[N[0,j]:N[1,j]+N[0,j],j])
            r[N[0,j]:N[1,j]+N[0,j],j]=np.divide(kde_y_obs1(y_initi[N[0,j]:N[1,j]+N[0,j],j]),kde_y_initi1(y_initi[N[0,j]:N[1,j]+N[0,j],j]))
            kde_y_obs2=GKDE(y_obs[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j])
            kde_y_initi2=GKDE(y_initi[N[0,j]+N[1,j]:N[1,j]+N[0,j]+N[2,j],j])
            r[N[0,j]+N[1,j]:N[1,j]+N[0,j]+N[2,j],j]=np.divide(kde_y_obs2(y_initi[N[0,j]+N[1,j]:N[1,j]+N[0,j]+N[2,j],j]),
                                                                                             kde_y_initi2(y_initi[N[0,j]+N[1,j]:N[1,j]+N[0,j]+N[2,j],j]))


def rejection_sampling(r): # creating indexes for samples to keep from initial samples using rejection sampling 
    N = r.size
    #seedinit2 = seedinit1+1
    np.random.seed(seedinit1)
    reject_prob = np.random.uniform(low=0, high=1, size=N)
    r = r/np.max(r)
    idx = np.where(r >= reject_prob)[0]
    return idx



samples_to_keepSize=np.zeros((3,nIter),dtype=int)
accept_rate=np.zeros((3,nIter))
for j in range(nIter):
    samples_to_keepSize[0,j] = rejection_sampling(r[0:N[0,j],j]).size
    samples_to_keepSize[1,j]=rejection_sampling(r[N[0,j]:N[1,j]+N[0,j],j]).size
    samples_to_keepSize[2,j] = rejection_sampling(r[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j]).size    


    accept_rate[0,j] = samples_to_keepSize[0,j]/N[0,j]
    accept_rate[1,j] = samples_to_keepSize[1,j]/N[1,j]
    accept_rate[2,j] = samples_to_keepSize[2,j]/N[2,j]

samples_to_keep=np.zeros((nSamp,nIter),dtype=int)

for j in range(nIter):
    samples_to_keep[0:samples_to_keepSize[0,j],j]= rejection_sampling(r[0:N[0,j],j])
    samples_to_keep[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j]= rejection_sampling(r[N[0,j]:N[0,j]+N[1,j],j])+N[0,j]
    samples_to_keep[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j]=rejection_sampling(
            r[N[0,j]+N[1,j]:N[0,j]+N[1,j]+N[2,j],j])+N[0,j]+N[1,j]



updated_Beta=np.zeros((nSamp,nIter))
updated_y=np.zeros((nSamp,nIter))
updated_eps=np.zeros((nSamp,nIter))

for j in range(nIter):
    updated_Beta[0:samples_to_keepSize[0,j],j] =Beta[samples_to_keep[0:samples_to_keepSize[0,j],j],j]
    updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j] =Beta[samples_to_keep[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],j]
    updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j] =Beta[samples_to_keep[samples_to_keepSize[0,j]
              +samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j],j]

    updated_y[0:samples_to_keepSize[0,j],j] =y_initi[samples_to_keep[0:samples_to_keepSize[0,j],j],j]
    updated_y[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j] =y_initi[samples_to_keep[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],j]
    updated_y[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j] =y_initi[samples_to_keep[samples_to_keepSize[0,j]+
        samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j],j]

    updated_eps[0:samples_to_keepSize[0,j],j] =eps[samples_to_keep[0:samples_to_keepSize[0,j],j],j]
    updated_eps[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j] =eps[samples_to_keep[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],j]
    updated_eps[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j] =eps[samples_to_keep[samples_to_keepSize[0,j]+
        samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j],j]


KruskalStat=np.zeros(nIter)
KruskalPvalue=np.zeros(nIter)
for j in range(nIter):
    KruskalStat[j]=stats.kruskal(updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],
                   updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j])[0]
    KruskalPvalue[j]=stats.kruskal(updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],
                      updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j])[1]


np.savetxt('KruskalPvalueE.txt',KruskalPvalue)



mannwhitneyPvalueUpB0_B1=np.zeros(nIter)
mannwhitneyPvalueUpB1_B2=np.zeros(nIter)
mannwhitneyPvalueUpB0_B2=np.zeros(nIter)

for j in range(nIter):

        mannwhitneyPvalueUpB0_B1[j]=sp.posthoc_mannwhitney([updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+
                                     samples_to_keepSize[1,j],j],updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+
                                      samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j]], p_adjust = 'bonferroni').loc[1,2]
        mannwhitneyPvalueUpB1_B2[j]=sp.posthoc_mannwhitney([updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+
                                     samples_to_keepSize[1,j],j],updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+
                                         samples_to_keepSize[2,j],j]], p_adjust = 'bonferroni').loc[2,3]
        mannwhitneyPvalueUpB0_B2[j]=sp.posthoc_mannwhitney([updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+
                                     samples_to_keepSize[1,j],j],updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+
                                         samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j]], p_adjust = 'bonferroni').loc[1,3]


np.savetxt('mannwhitneyPvalueUpB0_B1E.txt',mannwhitneyPvalueUpB0_B1)
np.savetxt('mannwhitneyPvalueUpB1_B2E.txt',mannwhitneyPvalueUpB1_B2)
np.savetxt('mannwhitneyPvalueUpB0_B2E.txt',mannwhitneyPvalueUpB0_B2)



Xs=np.zeros((nSamp,nIter))
y=dataNHW["Airway-Predominant Axis"].values

for j in range(nIter):
    Xs[:,j]=NHWdf.iloc[:,j].values

N=np.zeros((3,nIter),dtype=int)

for j in range(nIter):
    for i in range(nSamp):
        if Xs[i,j]==0 :
            N[0,j]+= 1
        elif Xs[i,j]==1:
            N[1,j] += 1
        else:
            N[2,j]+= 1

y_obs= np.zeros((nSamp,nIter))


for j in range(nIter):
    index0=np.where(Xs[:,j] == 0)
    index1=np.where(Xs[:,j] == 1)
    index2=np.where(Xs[:,j] == 2)
    indexnan=np.argwhere(np.isnan(Xs[:,j]))
    Xs0=Xs[:,j][index0]
    Xs1=Xs[:,j][index1]
    Xs2=Xs[:,j][index2]
    Xsnan=Xs[:,j][indexnan]
    Xs[:,j]=np.concatenate((Xs0,Xs1,Xs2,Xsnan), axis=None)
    y_obs0=y[index0]
    y_obs1=y[index1]
    y_obs2=y[index2]
    y_obsnan=y[indexnan]
    y_obs[:,j]=np.concatenate((y_obs0,y_obs1,y_obs2,y_obsnan), axis=None)


np.savetxt(' y_obs.txt', y_obs)
np.savetxt('Xs.txt',Xs)


Xs_bar = np.zeros(nIter)
SEb = np.zeros(nIter)
SEe=np.zeros(nIter)
for j in range(nIter):
    Xs_bar[j]=np.mean(Xs[:,j])
    Sxx = np.sum(Xs[:,j]*Xs[:,j])-nSamp*Xs_bar[j]*Xs_bar[j]
    sigma=.9
    SEb[j]= np.sqrt(sigma/Sxx)
    SEe[j]=.9-SEb[j]

OLS_Pvalues=np.zeros(nIter)
OLS_Beta=np.zeros(nIter)
for j in range(nIter):
    X=sm.add_constant(Xs[:,j])
    est=sm.OLS(y_obs[:,j],X,missing='drop').fit()
    OLS_Pvalues[j]=est.pvalues[1]
    OLS_Beta[j]=est2.params[1]

np.savetxt('OLS_PvaluesA.txt',OLS_Pvalues,fmt='%1.4e')
np.savetxt('OLS_BetaA.txt',OLS_Beta)


y_obs_bar = np.zeros((3,nIter))
SEy=np.zeros((3,nIter))
Sy=np.zeros(3)

for j in range(nIter):
        y_obs_bar[0,j]=np.mean(y_obs[0:N[0,j],j])
        y_obs_bar[1,j]=np.mean(y_obs[N[0,j]:N[1,j]+N[0,j],j])
        y_obs_bar[2,j]=np.mean(y_obs[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j])
        Sy[0]= np.sum(y_obs[0:N[0,j],j]*y_obs[0:N[0,j],j])-N[0,j]*y_obs_bar[0,j]*y_obs_bar[0,j]
        Sy[1]= np.sum(y_obs[N[0,j]:N[0,j]+N[1,j],j]*y_obs[N[0,j]:N[1,j]+N[0,j],j])-N[1,j]*y_obs_bar[1,j]*y_obs_bar[1,j]
        Sy[2]= np.sum(y_obs[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j]*y_obs[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j])-N[2,j]*y_obs_bar[2,j]*y_obs_bar[2,j]
        SEy[0,j]= np.sqrt(Sy[0]/(N[0,j]-1))
        SEy[1,j]= np.sqrt(Sy[1]/(N[1,j]-1))
        SEy[2,j]= np.sqrt(Sy[2]/(N[2,j]-1))

Beta = np.zeros((nSamp,nIter))
eps = np.zeros((nSamp,nIter))
y_initi= np.zeros((nSamp,nIter))
seedinit1 = 1

for j in range(nIter):
    #seedinit1 = seedinit1+j
    #np.random.seed(seedinit1)
    Beta[0:N[0,j],j] = norm.rvs(loc=y_obs_bar[0,j] , scale=SEy[0,j], size=N[0,j])
    Beta[N[0,j]:N[1,j]+N[0,j],j] = norm.rvs(loc=y_obs_bar[1,j], scale=SEy[1,j], size=N[1,j])
    Beta[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j] = norm.rvs(loc=y_obs_bar[2,j] , scale=SEy[2,j], size=N[2,j])

    eps[0:N[0,j],j]=np.random.normal(0,.34,size=N[0,j])
    eps[N[0,j]:N[1,j]+N[0,j],j]=np.random.normal(0,.34,size=N[1,j])
    eps[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j]=np.random.normal(0,.34,size=N[2,j])

y_initi=Beta+eps


r = np.zeros((nSamp,nIter))
for j in range(nIter):
            kde_y_obs0=GKDE(y_obs[0:N[0,j],j])
            kde_y_initi0=GKDE(y_initi[0:N[0,j],j])
            r[0:N[0,j],j]=np.divide(kde_y_obs0(y_initi[0:N[0,j],j]),kde_y_initi0(y_initi[0:int(N[0,j]),j]))
            kde_y_obs1=GKDE(y_obs[N[0,j]:N[1,j]+N[0,j],j])
            kde_y_initi1=GKDE(y_initi[N[0,j]:N[1,j]+N[0,j],j])
            r[N[0,j]:N[1,j]+N[0,j],j]=np.divide(kde_y_obs1(y_initi[N[0,j]:N[1,j]+N[0,j],j]),kde_y_initi1(y_initi[N[0,j]:N[1,j]+N[0,j],j]))
            kde_y_obs2=GKDE(y_obs[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j])
            kde_y_initi2=GKDE(y_initi[N[0,j]+N[1,j]:N[1,j]+N[0,j]+N[2,j],j])
            r[N[0,j]+N[1,j]:N[1,j]+N[0,j]+N[2,j],j]=np.divide(kde_y_obs2(y_initi[N[0,j]+N[1,j]:N[1,j]+N[0,j]+N[2,j],j]),
                                                                                             kde_y_initi2(y_initi[N[0,j]+N[1,j]:N[1,j]+N[0,j]+N[2,j],j]))


def rejection_sampling(r): # creating indexes for samples to keep from initial samples using rejection sampling 
    N = r.size
    #seedinit2 = seedinit1+1
    np.random.seed(seedinit1)
    reject_prob = np.random.uniform(low=0, high=1, size=N)
    r = r/np.max(r)
    idx = np.where(r >= reject_prob)[0]
    return idx



samples_to_keepSize=np.zeros((3,nIter),dtype=int)
accept_rate=np.zeros((3,nIter))
for j in range(nIter):
    samples_to_keepSize[0,j] = rejection_sampling(r[0:N[0,j],j]).size
    samples_to_keepSize[1,j]=rejection_sampling(r[N[0,j]:N[1,j]+N[0,j],j]).size
    samples_to_keepSize[2,j] = rejection_sampling(r[N[1,j]+N[0,j]:N[1,j]+N[0,j]+N[2,j],j]).size    


    accept_rate[0,j] = samples_to_keepSize[0,j]/N[0,j]
    accept_rate[1,j] = samples_to_keepSize[1,j]/N[1,j]
    accept_rate[2,j] = samples_to_keepSize[2,j]/N[2,j]

samples_to_keep=np.zeros((nSamp,nIter),dtype=int)

for j in range(nIter):
    samples_to_keep[0:samples_to_keepSize[0,j],j]= rejection_sampling(r[0:N[0,j],j])
    samples_to_keep[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j]= rejection_sampling(r[N[0,j]:N[0,j]+N[1,j],j])+N[0,j]
    samples_to_keep[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j]=rejection_sampling(
            r[N[0,j]+N[1,j]:N[0,j]+N[1,j]+N[2,j],j])+N[0,j]+N[1,j]



updated_Beta=np.zeros((nSamp,nIter))
updated_y=np.zeros((nSamp,nIter))
updated_eps=np.zeros((nSamp,nIter))

for j in range(nIter):
    updated_Beta[0:samples_to_keepSize[0,j],j] =Beta[samples_to_keep[0:samples_to_keepSize[0,j],j],j]
    updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j] =Beta[samples_to_keep[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],j]
    updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j] =Beta[samples_to_keep[samples_to_keepSize[0,j]
              +samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j],j]

    updated_y[0:samples_to_keepSize[0,j],j] =y_initi[samples_to_keep[0:samples_to_keepSize[0,j],j],j]
    updated_y[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j] =y_initi[samples_to_keep[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],j]
    updated_y[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j] =y_initi[samples_to_keep[samples_to_keepSize[0,j]+
        samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j],j]

    updated_eps[0:samples_to_keepSize[0,j],j] =eps[samples_to_keep[0:samples_to_keepSize[0,j],j],j]
    updated_eps[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j] =eps[samples_to_keep[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],j]
    updated_eps[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j] =eps[samples_to_keep[samples_to_keepSize[0,j]+
        samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j],j]


KruskalStat=np.zeros(nIter)
KruskalPvalue=np.zeros(nIter)
for j in range(nIter):
    KruskalStat[j]=stats.kruskal(updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],
                   updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j])[0]
    KruskalPvalue[j]=stats.kruskal(updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j],j],
                      updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j])[1]


np.savetxt('KruskalPvalueA.txt',KruskalPvalue)



mannwhitneyPvalueUpB0_B1=np.zeros(nIter)
mannwhitneyPvalueUpB1_B2=np.zeros(nIter)
mannwhitneyPvalueUpB0_B2=np.zeros(nIter)

for j in range(nIter):

        mannwhitneyPvalueUpB0_B1[j]=sp.posthoc_mannwhitney([updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+
                                     samples_to_keepSize[1,j],j],updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+
                                      samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j]], p_adjust = 'bonferroni').loc[1,2]
        mannwhitneyPvalueUpB1_B2[j]=sp.posthoc_mannwhitney([updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+
                                     samples_to_keepSize[1,j],j],updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+samples_to_keepSize[1,j]+
                                         samples_to_keepSize[2,j],j]], p_adjust = 'bonferroni').loc[2,3]
        mannwhitneyPvalueUpB0_B2[j]=sp.posthoc_mannwhitney([updated_Beta[0:samples_to_keepSize[0,j],j], updated_Beta[samples_to_keepSize[0,j]:samples_to_keepSize[0,j]+
                                     samples_to_keepSize[1,j],j],updated_Beta[samples_to_keepSize[0,j]+samples_to_keepSize[1,j]:samples_to_keepSize[0,j]+
                                         samples_to_keepSize[1,j]+samples_to_keepSize[2,j],j]], p_adjust = 'bonferroni').loc[1,3]


np.savetxt('mannwhitneyPvalueUpB0_B1A.txt',mannwhitneyPvalueUpB0_B1)
np.savetxt('mannwhitneyPvalueUpB1_B2A.txt',mannwhitneyPvalueUpB1_B2)
np.savetxt('mannwhitneyPvalueUpB0_B2A.txt',mannwhitneyPvalueUpB0_B2)


