import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sim_functions as sim

import sys
sys.path.insert(0, '../')
import tools as tools

plt.ion()

filename='compare_radiation_proton.dat'

file=open(filename,'r')
info=cPickle.load(file)
file.close()

print info.keys()

    
em_energy=info['em_energy']
energy=info['energy']
zenith=info['zenith']
azimuth=info['azimuth']
xmax=info['xmax']
alpha=info['alpha']
erad_ldf=info['erad_ldf']
erad_coreas=info['erad_coreas']
erad_int=info['erad_int']
res_all=info['res_all']

print len(em_energy)

good=res_all

diff1=2*(erad_coreas-erad_ldf)/(erad_coreas+erad_ldf)


binWidth=0.05
bin_start=-2.0
bin_stop=2.0
nBins=int((bin_stop-bin_start)/binWidth)
bins=np.arange(bin_start,bin_stop,binWidth)

print nBins

p0,p1,p2, x_hist, y_hist=tools.fit_hist(diff1,bin_start,bin_stop,nBins,-1.)


fig = plt.figure()
ax1 = fig.add_subplot(1,2,1,aspect=1)
ax2 = fig.add_subplot(1,2,2)

line=np.arange(1e3,1e8,1e7)

ax1.plot(line,line,color='black')
#ax2.plot(line,line,color='black')
#ax3.plot(line,line,color='black')

ax1.plot(erad_coreas,erad_ldf,'.',color='black',alpha=0.1)
#ax1.plot(erad_coreas,erad_ldf,'.')

ax1.set_yscale('log')
ax1.set_xscale('log')

ax2.hist(diff1,alpha=0.7,bins=bins)
ax2.plot(x_hist, y_hist,color='green',linewidth=2,label='mean {0:.2f}\n spread {1:.2f}'.format(p1,p2))
#ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel('coreas radiation energy')
ax1.set_ylabel('LDF radiation energy')
ax2.set_xlabel('2*(coreas-ldf)/(coreas+ldf)')
#ax2.set_ylabel('S$_{RD} (eV)$',fontsize=12)
#ax2.axis([2e16,3e18,7e3,1e8])

ax1.grid()
#ax2.grid()
plt.show()
raw_input()
plt.close()
