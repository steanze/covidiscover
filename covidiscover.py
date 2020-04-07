# IMPORT LIBRARIES

from scipy.stats import beta
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import root

# DEFINE FUNCTIONS

def sfunc(s):
  y = (M-K/(1-betacdfCovid(s)))*(1-betacdfNoCovid(s)) - N+K
  return y

def covidiscover(M,N,K,ac,bc,anc,bnc):
  betacdfCovid = beta(ac,bc).cdf
  betacdfNoCovid = beta(anc,bnc).cdf
  def sfunc(s):
    y = (M-K/(1-betacdfCovid(s)))*(1-betacdfNoCovid(s)) - N+K
    return y
  s0 = root(sfunc,0.5)
  pcovid= K/(M*(1-betacdfCovid(s0.x)))
  totcovid = M*pcovid
  return pcovid, totcovid

def covidsim(M,X,N,ac,bc,anc,bnc):
  symptoms_covid = np.random.beta(ac,bc,int(X))
  symptoms_nocovid = np.random.beta(anc,bnc,int(M-X))
  symptoms = np.hstack((symptoms_covid,symptoms_nocovid))
  i = np.argsort(-symptoms)
  positivity = np.hstack((np.ones(len(symptoms_covid)),np.zeros(len(symptoms_nocovid))))
  positivity = positivity[i]
  K = np.sum(positivity[0:int(N)])
  return K

# MAKE FIGURES OF SYMPTOM DISTRIBUTIONS

fig= plt.figure(figsize=(8,3))

acac = [2,3,4]
bcbc = [2,3,4]
x = np.linspace(0,1,1000)

ancanc = [1]
bncbnc = [5,6,7]

for i,ac in enumerate(acac):
  for j,bc in enumerate(bcbc):
    betaCovid = beta(ac,bc).pdf
    y1 = betaCovid(x)
    c1 = 0.2+0.2*i
    c2 = 0.2+0.2*j
    plt.plot(x,y1,label="covid",color = (c1,0,c2))
i = 1
for anc in ancanc:
  for bnc in bncbnc:
    betaNoCovid = beta(anc,bnc).pdf
    y2 = betaNoCovid(x)
    c = 0.2+0.2*i
    plt.plot(x,y2,label="no covid",color = (0,c,0))
    i=i+1
matplotlib.rc('xtick', labelsize=18) 
matplotlib.rc('ytick', labelsize=18) 

# RUN SIMULATIONS

MM = [1E5,1E6,1E7]
XX = [1E2,1E3,1E4]
NN = [1E2,1E3,1E4]
acac = [2,3,4]
bcbc = [2,3,4]
ancanc = [1]
bncbnc = [5,6,7]
niter = 10

Xe = np.zeros([len(MM),len(XX),len(NN),len(acac),len(bcbc),len(ancanc),len(bncbnc),niter])
for iM, M in enumerate(MM):
  for iX, X in enumerate(XX):
    for iN, N in enumerate(NN):
      for iac, ac in enumerate(acac):
        for ibc, bc in enumerate(bcbc):
          for ianc, anc in enumerate(ancanc):
            for ibnc, bnc in enumerate(bncbnc):
              for i in range(niter):
                K = covidsim(M,X,N,ac,bc,anc,bnc)
                pcovid,totcovid = covidiscover(M,N,K,ac,bc,anc,bnc)
                Xe[iM,iX,iN,iac,ibc,ianc,ibnc,i] = totcovid

# MAKE FIGURES OF SIMULATION RESULTS

for iM, M in enumerate(MM):
  for iN, N in enumerate(NN):
    for iX, X in enumerate(XX):
      for iac, ac in enumerate(acac):
        for ibc, bc in enumerate(bcbc):
          for ianc, anc in enumerate(ancanc):
            for ibnc, bnc in enumerate(bncbnc):
              c1 = 0.2+0.2*iac
              c2 = 0.2+0.2*ibc
              plt.scatter(np.ones(niter)*np.log10(X),np.log10(Xe[iM,iX,iN,iac,ibc,ianc,ibnc,:]),c = [[c1,0,c2]])
    plt.hlines([2,3,4],1,5, linestyles='dashed')
    plt.xticks(np.arange(3)+2)
    plt.yticks(np.arange(7)+1)
    fname = 'scatter_M'+str(iM)+'_N'+str(iN)
    plt.savefig(fname, dpi=300)
    plt.clf()
