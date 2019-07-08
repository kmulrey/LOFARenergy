#from pycrtools import crdatabase as crdb
import numpy as np
from optparse import OptionParser
import os
import cPickle
import re




CR_events=np.genfromtxt(open('energy_events.txt','r'))
dir='/vol/astro7/lofar/sim/pipeline/production_mcvsmc_radio_only_oct2018/'
dir='/Users/kmulrey/LOFAR/energy/LOFARenergy/radio/resolution/'


p_ratios=[]
d_ratios=[]


    #for i in np.arange(len(CR_events)):
for i in np.arange(1):
    
    
    event=120768260
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

       
        selection_valid_reconstructions = (xreco > 1)
        nfailed = np.sum(xreco < 1)

        sel = selection_valid_reconstructions
        sigma_e = np.std(d_ratio[sel])
        sigma_e_radio = np.std(np.sqrt(p_ratio[sel]))
        sigma_logE_radio = np.std(0.5 * np.log10(p_ratio[sel]))
        sigma_logE = np.std(0.5 * np.log10(d_ratio[sel]))

        r_recon=np.sqrt(xoff**2+yoff**2)
        sigma_r=np.std(r_recon[sel])

        event_list.append(event)
        sigma_e_list.append(sigma_e)
        sigma_e_radio_list.append(sigma_e_radio)
        sigma_loge_list.append(sigma_logE)
        sigma_loge_radio_list.append(sigma_logE_radio)
        sigma_r_list.append(sigma_r)

    else:
        print '---> didn\'t find event {0} \n'.format(event)
'''

info={'event':np.asarray(event_list),'sigma_e':np.asarray(sigma_e_list),'sigma_e_radip':np.asarray(sigma_e_radio_list),'sigma_logE':np.asarray(sigma_loge_list),'sigma_logE_radio':np.asarray(sigma_loge_radio_list),'sigma_r':np.asarray(sigma_r_list)}

outfile=open('uncertainty_info.dat','w')
cPickle.dump(info,outfile)
outfile.close()

'''
