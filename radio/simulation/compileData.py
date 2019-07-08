import numpy as np
import os
import cPickle


data_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/radio/simulation/output/'
#data_dir='output/'

file_list=[]

for file in os.listdir(data_dir):
    if file.startswith('proton'):
        file_list.append(data_dir+file)

print file_list

outfilename='compiled_sims_proton.dat'

em_energy=np.empty([0])
other_energy=np.empty([0])
total_energy=np.empty([0])
zenith=np.empty([0])
azimuth=np.empty([0])
energy=np.empty([0])
xmax=np.empty([0])
alpha=np.empty([0])
Erad=np.empty([0])
clip=np.empty([0])
charge_excess_ratio=np.empty([0])
S_basic=np.empty([0])
Srd_1=np.empty([0])
Srd_2=np.empty([0])
StoEM=np.empty([0])





for i in np.arange(len(file_list)):
    infile=open(file_list[i],'r')
    info=cPickle.load(infile)
    infile.close()

    em_energy_hold=info['em_energy']
    other_energy_hold=info['other_energy']
    total_energy_hold=info['total_energy']
    zenith_hold=info['zenith']
    azimuth_hold=info['azimuth']
    energy_hold=info['energy']
    xmax_hold=info['xmax']
    alpha_hold=info['alpha']
    Erad_hold=info['Erad']
    clip_hold=info['clip']
    charge_excess_ratio_hold=info['charge_excess_ratio']
    S_basic_hold=info['S_basic']
    Srd_1_hold=info['Srd_1']
    Srd_2_hold=info['Srd_2']
    StoEM_hold=info['StoEM']
    print xmax_hold

    
    try:    
    	energy=np.concatenate((energy,energy_hold))
        em_energy=np.concatenate((em_energy,em_energy_hold))
        other_energy=np.concatenate((other_energy,other_energy_hold))
        total_energy=np.concatenate((total_energy,total_energy_hold))
        zenith=np.concatenate((zenith,zenith_hold))
        azimuth=np.concatenate((azimuth,azimuth_hold))
        xmax=np.concatenate((xmax,xmax_hold))
        alpha=np.concatenate((alpha,alpha_hold))
        Erad=np.concatenate((Erad,Erad_hold))
        clip=np.concatenate((clip,clip_hold))
        charge_excess_ratio=np.concatenate((charge_excess_ratio,charge_excess_ratio_hold))
        S_basic=np.concatenate((S_basic,S_basic_hold))
        Srd_1=np.concatenate((Srd_1,Srd_1_hold))
        Srd_2=np.concatenate((Srd_2,Srd_2_hold))
        StoEM=np.concatenate((StoEM,StoEM_hold))

   
    except:
	print 'oops'



info={'em_energy':np.asarray(em_energy),
    'other_energy':np.asarray(other_energy),
    'total_energy':np.asarray(total_energy),
    'zenith':np.asarray(zenith),
    'azimuth':np.asarray(azimuth),
    'energy':np.asarray(energy),
    'xmax':np.asarray(xmax),
    'alpha':np.asarray(alpha),
    'Erad':np.asarray(Erad),
    'clip':np.asarray(clip),
    'charge_excess_ratio':np.asarray(charge_excess_ratio),
    'S_basic':np.asarray(S_basic),
    'Srd_1':np.asarray(Srd_1),
    'Srd_2':np.asarray(Srd_2),
    'StoEM':np.asarray(StoEM),
}



outfile=open(outfilename,'w')
cPickle.dump(info,outfile)
outfile.close()









