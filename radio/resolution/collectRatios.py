#from pycrtools import crdatabase as crdb
import numpy as np
from optparse import OptionParser
import os
import cPickle
import re




CR_events=np.genfromtxt(open('energy_events.txt','r'))
dir='/vol/astro7/lofar/sim/pipeline/production_mcvsmc_radio_only_oct2018/'
#dir='/Users/kmulrey/LOFAR/energy/LOFARenergy/radio/resolution/'


p_ratios=[]
d_ratios=[]


#for i in np.arange(len(CR_events)):
for i in np.arange(1,100):
    
    
    event=int(CR_events[i])
    print event
    check_file=0
    
    
    for file in os.listdir(dir):
        if str(event) in file:
            filename=dir+file
            check_file=1

    if check_file==1:
        print 'found event {0}: {2} / {1} \n'.format(event,len(CR_events),i+1)
        f=open(filename)
        mcvsmc_info = cPickle.load(f) # contains tuple

        (nsimprot, nsimiron, ev, xreco, xreal, xbest, cchi2, rchi2, p_ratio, p_ratio0, p_ratio1, d_ratio, dratio, sim_tot_power, rbf, d150, nch, mr, sa, xoff, yoff, xcore, ycore, en, zen, az, noisepower, antratio, nsim, nsim_prot) = mcvsmc_info
        
        print p_ratio
        print d_ratio
        for p in np.arange(len(p_ratio)):
            p_ratios.append(p_ratio[p])
    
        for d in np.arange(len(d_ratio)):
            d_ratios.append(d_ratio[d])

    else:
        print '---> didn\'t find event {0} \n'.format(event)


info={'p_ratio':np.asarray([p_ratios]),'d_ratio':np.asarray([d_ratios])}

outfile=open('resolution_info.dat','w')
cPickle.dump(info,outfile)
outfile.close()


