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
from scipy.stats import uniform, norm, beta,skewnorm,gamma, beta, f_oneway, wilcoxon,kruskal, ranksums
from mpl_toolkits import mplot3d
import scipy.io as sio
from statsmodels.distributions.empirical_distribution import ECDF
import scikit_posthocs as sp


# Generating SNP vector for 15% < MAF < 50% 

nSamp = 10000
seedinit1 = 1
as1=np.zeros(nSamp)
as2=np.zeros(nSamp)
np.random.seed(seedinit1)
SV=np.random.uniform(0.15,0.5,1)   
unifs1=np.random.uniform(0,1,nSamp)
for i in range(nSamp):
    if unifs1[i] <= SV:
        as1[i]=1

seedinit2 = seedinit1+1
np.random.seed(seedinit2)            
unifs2=np.random.uniform(0,1,nSamp)
for j in range(nSamp):
    if unifs2[j] <= SV:
        as2[j]=1
    
Xs=as1+as2


# Simulating the vector of phynotype for each SNP group for additive Normal Case

seedinit3 = seedinit2+1
np.random.seed(seedinit3)

# Generate the mean for each phynotype group
B0= skewnorm.rvs(.5, loc=-.2718 , scale=.8, size=N0)
B1= skewnorm.rvs(-2,loc=.8073 , scale=1, size=N1)
B2= skewnorm.rvs(2,loc=-.77574, scale=1, size=N2)

# Generate the error (noise) for each phynotype group
e0=np.random.normal(0, SEb,size=N0)
e1=np.random.normal(0, SEb,size=N1)
e2=np.random.normal(0, SEb,size=N2)

# Simulated phynotype groups using generated mean and error
y_obs0=B0+e0
y_obs1=B1+e1
y_obs2=B2+e2

# Finding indices for each SNP group
index0=np.where(Xs == 0)
index1=np.where(Xs == 1)
index2=np.where(Xs == 2)

# Separating SNP vector into three SNPs group 
Xs0=Xs[index0]
Xs1=Xs[index1]
Xs2=Xs[index2]

# Estimating the probability density function for each simulated phynotype group
kde_y_obs0= GKDE(y_obs0)
kde_y_obs1= GKDE(y_obs1)
kde_y_obs2= GKDE(y_obs2)


# Generating initial Betas for each group
seedinit4 = seedinit3+1
np.random.seed(seedinit4)
Beta0=norm.rvs(loc=y_obs0_bar, scale=SEy0,size=N0)
Beta1=norm.rvs(loc=y_obs1_bar, scale=SEy1,size=N1)
Beta2=norm.rvs(loc=y_obs2_bar, scale=SEy2,size=N2)


# Generating initial errors for each group
seedinit5 = seedinit4+1
np.random.seed(seedinit5)
eps0=np.random.normal(0,.35,size=N0)
eps1=np.random.normal(0,.35,size=N1)
eps2=np.random.normal(0,.35,size=N2)


# Initial phynotype for each group
y_initi0=Beta0+eps0
y_initi1=Beta1+eps1
y_initi2=Beta2+eps2


# Estimating the probability density function for each initial phynotype group
kde_y_initi0 = GKDE(y_initi0)
kde_y_initi1 = GKDE(y_initi1)
kde_y_initi2 = GKDE(y_initi2)


# Rejection sampling function to update initials
# Creating indexes for samples to keep from initial samples using rejection sampling 

def rejection_sampling(r):
    N = r.size 
    seedinit6 = seedinit5+1
    np.random.seed(seedinit6)
    reject_prob = np.random.uniform(low=0, high=1, size=N) 
    r = r/np.max(r)
    idx = np.where(r >= reject_prob)[0]
    return idx

# Computing RN weight

r0=np.divide(kde_y_obs0(y_initi0),kde_y_initi0(y_initi0))
r1=np.divide(kde_y_obs1(y_initi1),kde_y_initi1(y_initi1))
r2=np.divide(kde_y_obs2(y_initi2),kde_y_initi2(y_initi2))

# Performing rejection sampling
samples_to_keep0 = rejection_sampling(r0) 
samples_to_keep1 = rejection_sampling(r1)
samples_to_keep2 = rejection_sampling(r2)

# Computing acceptance rate from performing rejection sampling
accept_rate0 = samples_to_keep0.size/Beta0.shape[0]
accept_rate1 = samples_to_keep1.size/Beta1.shape[0]
accept_rate2 = samples_to_keep2.size/Beta2.shape[0]

# Samples of updated Beta
updated_Beta0 =Beta0[samples_to_keep0] 
updated_Beta1 =Beta1[samples_to_keep1]
updated_Beta2 =Beta2[samples_to_keep2]

# Updated Beta density
updated_y0 =y_initi0[samples_to_keep0] 
updated_y1 =y_initi1[samples_to_keep1] 
updated_y2 =y_initi2[samples_to_keep2] 

# Samples of updated phynotype
kde_updated_y0 = GKDE(y_initi0, weights=r0)
kde_updated_y1 = GKDE(y_initi1, weights=r1)
kde_updated_y2 = GKDE(y_initi2, weights=r2)

# Updated phynotype density
kde_updated_Beta0 = GKDE(Beta0, weights=r0)
kde_updated_Beta1 = GKDE(Beta1, weights=r1)
kde_updated_Beta2 = GKDE(Beta2, weights=r2)

# Samples of updated error
updated_eps0 =eps0[samples_to_keep0]
updated_eps1 =eps1[samples_to_keep1] 
updated_eps2 =eps2[samples_to_keep2]

# Updated error density
kde_updated_eps0 = GKDE(eps0, weights=r0)
kde_updated_eps1 = GKDE(eps1, weights=r1)
kde_updated_eps2 = GKDE(eps2, weights=r2)


# Performing kruskal wallis test on updated means
stats.kruskal(updated_Beta0, updated_Beta1, updated_Beta2)


# Performing F-test test on updated means
f_oneway(updated_Beta0, updated_Beta1, updated_Beta2)

