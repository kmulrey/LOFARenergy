import numpy as np
import os
import cPickle


data_dir='/vol/astro3/lofar/sim/kmulrey/energy/LOFARenergy/radio/lofar_events/output/'
#data_dir='output/'

file_list=[]

for file in os.listdir(data_dir):
    file_list.append(data_dir+file)


outfilename='compiled_sims_all.dat'

em_energy=np.empty([0])
other_energy=np.empty([0])
total_energy=np.empty([0])
zenith=np.empty([0])
azimuth=np.empty([0])
energy=np.empty([0])
xmax=np.empty([0])
xmax_fit=np.empty([0])
core_x=np.empty([0])
core_y=np.empty([0])
p_ratio=np.empty([0])
d_ratio=np.empty([0])
combchi2=np.empty([0])
radiochi2=np.empty([0])
event=np.empty([0])
alpha=np.empty([0])
Erad=np.empty([0])
clip=np.empty([0])
charge_excess_ratio=np.empty([0])
S_basic=np.empty([0])
Srd_1=np.empty([0])
Srd_2=np.empty([0])
density=np.empty([0])
dmax=np.empty([0])
type=np.empty([0])


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
    xmax_fit_hold=info['xmax_fit']
    core_x_hold=info['core_x']
    core_y_hold=info['core_y']
    p_ratio_hold=info['p_ratio']
    d_ratio_hold=info['d_ratio']
    combchi2_hold=info['combchi2']
    radiochi2_hold=info['radiochi2']
    event_hold=info['event']
    density_hold=info['density']
    dmax_hold=info['dmax']
    type_hold=info['type']

    
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

        xmax_fit=np.concatenate((xmax_fit,xmax_fit_hold))
        core_x=np.concatenate((core_x,core_x_hold))
        core_y=np.concatenate((core_y,core_y_hold))
        p_ratio=np.concatenate((p_ratio,p_ratio_hold))
        d_ratio=np.concatenate((d_ratio,d_ratio_hold))
        combchi2=np.concatenate((combchi2,combchi2_hold))
        radiochi2=np.concatenate((radiochi2,radiochi2_hold))
        event=np.concatenate((event,event_hold))
        density=np.concatenate((density,density_hold))
        dmax=np.concatenate((dmax,dmax_hold))
        type=np.concatenate((type,type_hold))


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
    'xmax_fit':np.asarray(xmax_fit),
    'core_x':np.asarray(core_x),
    'core_y':np.asarray(core_y),
    'p_ratio':np.asarray(p_ratio),
    'd_ratio':np.asarray(d_ratio),
    'combchi2':np.asarray(combchi2),
    'radiochi2':np.asarray(radiochi2),
    'event':np.asarray(event),
    'density':np.asarray(density),
    'dmax':np.asarray(dmax),
    'type':np.asarray(type),

}



outfile=open(outfilename,'w')
cPickle.dump(info,outfile)
outfile.close()









