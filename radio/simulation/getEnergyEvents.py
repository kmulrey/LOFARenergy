import numpy as np
from optparse import OptionParser
import os
import cPickle
import os.path
import matplotlib.pyplot as plt
import sim_functions as sim
import scipy.interpolate as intp
from scipy.interpolate import interp1d
from matplotlib import cm
import re
from scipy import signal
from scipy.signal import hilbert
import fluence as flu
import radiation_energy as rad
import helper as helper

plt.ion()


#reco_dir='/vol/astro7/lofar/sim/pipeline/production_analysis_oct2018/'
event_path='/vol/astro7/lofar/sim/pipeline/events/'




parser = OptionParser()
parser.add_option('-a', '--start',type='int',help='line number of event to start',default=0)
#parser.add_option('-b', '--stop',type='int',help='line number of event to stop',default=0)
parser.add_option('-t', '--type',type='string',help='type',default='proton')
#parser.add_option('-o', '--output',type='string',help='output',default='event_information.dat')

(options, args) = parser.parse_args()
a=options.start
#b=options.stop
#outfilename=options.output
type=options.type

start_id=a*10
stop_id=start_id*10+10
outfilename='output/'+type+'_'+str(start_id)+'_'+str(stop_id)+'.dat'

# info to record

em_energy_list=[]
other_energy_list=[]
total_energy_list=[]
zenith_list=[]
az_list=[]
energy_list=[]
xmax_list=[]
alpha_list=[]
Erad_list=[]
clip_ratio_list=[]
charge_excess_ratio_list=[]
S_basic_list=[]
Srd_1_list=[]
Srd_2_list=[]
StoE_list=[]




filename=type+'_runs.txt'

sim_dir=[]
runnr=[]
count=0

with open(filename) as f:
    for line in f:
        sim_dir.append(line.split()[0])
        
        runnr.append(line.split()[1])
        count=count+1







for i in np.arange(start_id,stop_id):
    #for i in np.arange(2,3):
    #try:
    for r in np.arange(1):
        
        event= sim_dir[i].split('events/')[1].split('/')[0]
       
        ant_pos,times,efield,zenith,az,energy,xmax=sim.get_efield(sim_dir[i], runnr[i])
        em_energy_hold,other_energy_hold,total_energy_hold=sim.getEM(sim_dir[i], runnr[i])
        
        ant_pos_shower=sim.ground_to_shower(ant_pos,zenith,az,B=1.1837)
        atm=sim.get_atm(event)
        # returns antenna positions,efield,zenith,az_rot,energy (az=0=east, pi/2=north)
        # efield (3,nant,2, dlen)
        alpha=sim.GetAlpha(zenith,az,1.1837)
        
        filtered_efield=sim.filter(times,efield, 30.0, 80.0)


        fluence=flu.calculate_energy_fluence(filtered_efield, times, signal_window=100., remove_noise=True)

        Erad=rad.integral(fluence,ant_pos_shower)
        
        #print 'xmax: {0}'.format(xmax)
        hi=helper.get_vertical_height(xmax,atm)   # height in cm

        at=helper.get_atmosphere(hi,atm)
        rho=helper.get_density(hi, atm, allow_negative_heights=True)
        dmax=helper.get_distance_xmax_geometric(zenith, xmax, atm)

        clip_ratio=rad.get_clipping(dmax)
        
        a=rad.get_a(rho)
        #print '_________________________________________'
        #print 'xmax: {0}'.format(xmax)
        #print 'zenith: {0}'.format(zenith*180/np.pi)

        #print 'relative charge excess: {0}'.format(a)
        
        S=rad.get_S_basic(Erad, np.sin(alpha))
        Srd=rad.get_S(Erad, np.sin(alpha), rho)
        Srd_final=rad.get_radiation_energy(Srd, np.sin(alpha), rho)
        
        StoE=rad.StoEm(Srd_final)
        
        #print 'Er:    {0:.2f}'.format(Erad)
        #print 'S basic:    {0:.2f}'.format(S)
        #print 'Srd:        {0:.2f}'.format(Srd)
        #rint 'Srd final:  {0:.2f}'.format(Srd_final)
        #print 'Em Energy:  {0:.2f}'.format(em_energy_hold*1e9)
        #print 'S to Em Energy:  {0:.2f}'.format(StoE)


        em_energy_list.append(em_energy_hold*1e9)
        other_energy_list.append(other_energy_hold*1e9)
        total_energy_list.append(total_energy_hold*1e9)
        zenith_list.append(zenith)
        az_list.append(az)
        energy_list.append(energy*1e9)
        xmax_list.append(xmax)
        alpha_list.append(alpha)
        Erad_list.append(Erad)
        clip_ratio_list.append(clip_ratio)
        charge_excess_ratio_list.append(a)
        S_basic_list.append(S)
        Srd_1_list.append(Srd)
        Srd_2_list.append(Srd_final)
        StoE_list.append(StoE)



    #except:
#   print 'bad event'

info={'em_energy':np.asarray(em_energy_list),
    'other_energy':np.asarray(other_energy_list),
    'total_energy':np.asarray(total_energy_list),
    'zenith':np.asarray(zenith_list),
    'azimuth':np.asarray(az_list),
    'energy':np.asarray(energy_list),
    'xmax':np.asarray(xmax_list),
    'alpha':np.asarray(alpha_list),
    'Erad':np.asarray(Erad_list),
    'clip':np.asarray(clip_ratio_list),
    'charge_excess_ratio':np.asarray(charge_excess_ratio_list),
    'S_basic':np.asarray(S_basic_list),
    'Srd_1':np.asarray(Srd_1_list),
    'Srd_2':np.asarray(Srd_2_list),
    'StoEM':np.asarray(StoE_list),
}



outfile=open(outfilename,'w')
cPickle.dump(info,outfile)
outfile.close()


