import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sys
from ROOT import TH1F
from ROOT import TF1

plt.ion()

def fit_hist(counts,bin_start,bin_stop,nBins,min):
    
    counts=counts[counts>min]
    binWidth=(bin_stop-bin_start)/nBins
    bins=np.arange(bin_start,bin_stop,binWidth)
    Hist = TH1F("Hist", "Hist_x;x-axis;Frequency",int(nBins),int(bin_start),int(bin_stop))
    fit = TF1("fit", "gaus", -2,2)
    #fit = TF1("fit", "gaus", -200,10000)
    
    for i in np.arange(len(counts)):
        Hist.Fill(counts[i])

    Hist.Fit("fit",'0')

    p0=fit.GetParameter(0)
    p1=fit.GetParameter(1)
    p2=fit.GetParameter(2)

    x=np.arange(-100,10000,1)
    y=np.zeros([len(x)])
    for i in np.arange(len(x)):
        y[i]=fit.Eval(x[i])

    return p0, p1, p2, x, y

 

infile=open('resolution_info.dat','r')
info=cPickle.load(infile)
infile.close()
    
p_ratios=np.sqrt(info['p_ratio'].flatten())
d_ratios=info['d_ratio'].flatten()


bin_start=0.0
bin_stop=2.0
binWidth=0.05
#nBins=30.0

#binWidth=(bin_stop-bin_start)/nBins
nBins=int((bin_stop-bin_start)/binWidth)
bins=np.arange(bin_start,bin_stop,binWidth)

p0_p,p1_p,p2_p, x_p, y_p=fit_hist(p_ratios,bin_start,bin_stop,nBins,-1)
p0_d,p1_d,p2_d, x_d, y_d=fit_hist(d_ratios,bin_start,bin_stop,nBins,-1)


fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ax1.hist(p_ratios,bins=bins)
ax2.hist(d_ratios,bins=bins)


ax1.set_xlabel('$f_r$',fontsize=12)
ax2.set_xlabel('$f_p$',fontsize=12)

ax1.set_yscale('log')
ax2.set_yscale('log')




plt.show()
raw_input()
plt.close()
