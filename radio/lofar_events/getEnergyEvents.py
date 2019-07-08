import numpy as np
from optparse import OptionParser
import os
import cPickle
import os.path
import matplotlib.pyplot as plt
import scipy.interpolate as intp
from scipy.interpolate import interp1d
from matplotlib import cm
import re
import sys
sys.path.insert(0, '../simulation')

import fluence as flu
import radiation_energy as rad
import sim_functions as sim
import helper as helper

#import slope_functions as slope
#import density

freq_a=30.0
freq_b=77.0
nAnt=160
go_ahead=1

n_spoke=8
n_rad=20

dist=np.arange(25,n_rad*25+25,25)




parser = OptionParser()
parser.add_option('-e', '--event',type='int',help='line number of event',default=0)

(options, args) = parser.parse_args()
e=int(options.event)


reco_dir='/vol/astro7/lofar/sim/pipeline/production_analysis_oct2018/'
event_path='/vol/astro7/lofar/sim/pipeline/events/'

def integrate(r,flu0,flu1):
    n=len(r)
    dr=r[1]-r[0]
    integral=0
    for i in np.arange(n-1):
        r0=r[i]
        r1=r[i+1]
        val0=r0*(flu0[i]+flu1[i])
        val1=r1*(flu0[i+1]+flu1[i+1])
        integral=integral+(val0+val1)*0.5*dr
    
    
    
    return 2*np.pi*integral*6.2415e18 # to eV


def return_file_list(dirIron,dirProton,filelist,xmax_list):
    
    for file in os.listdir(dirIron):
        if file.endswith(".long"):
            longfile=dirIron+file
            runnr=file.split(".")[0]
            runnr=runnr.split("T")[1]
            sim_path=dirIron+'SIM'+runnr+'_coreas/'
            hillas = np.genfromtxt(re.findall("PARAMETERS.*",open(longfile,'r').read()))[2:]
            #print 'iron: {0}      {1}   {2}'.format(file,hillas[2],sim_path)
            filelist.append(sim_path)
            xmax_list.append(hillas[2])
    
    
    for file in os.listdir(dirProton):
        if file.endswith(".long"):
            longfile=dirProton+file
            runnr=file.split(".")[0]
            runnr=runnr.split("T")[1]
            sim_path=dirProton+'SIM'+runnr+'_coreas/'
            #print sim_path
            hillas = np.genfromtxt(re.findall("PARAMETERS.*",open(longfile,'r').read()))[2:]
            #print 'proton: {0}      {1}  {2}'.format(file,hillas[2],sim_path)
            filelist.append(sim_path)
            xmax_list.append(hillas[2])

    return filelist,xmax_list

def find_SIM_dir(event_no,xmax):
    
    event=str(int(event_no))
    filelist=[]
    xmax_list=[]
    
    try:
        dirIron=event_path+event+'/2/coreas/iron/'
        dirProton=event_path+event+'/2/coreas/proton/'
        filelist,xmax_list=return_file_list(dirIron,dirProton,filelist,xmax_list)
    except:
        #print 'no 2 directory'
        pass
    try:
        dirIron=event_path+event+'/1/coreas/iron/'
        dirProton=event_path+event+'/1/coreas/proton/'
        filelist,xmax_list=return_file_list(dirIron,dirProton,filelist,xmax_list)
    except:
        #print 'no 1 directory'
        pass
    
    try:
        dirIron=event_path+event+'/0/coreas/iron/'
        dirProton=event_path+event+'/0/coreas/proton/'
        #print dirProton
        filelist,xmax_list=return_file_list(dirIron,dirProton,filelist,xmax_list)
    except:
        #print 'no 0 directory'
        pass

    closest_xmax=100000
    closest_index=100000
    min_diff=10000
    

    for i in np.arange(len(xmax_list)):
        min_diff_here=np.abs(xmax_list[i]-xmax)
        if min_diff_here<min_diff:
            closest_xmax=xmax_list[i]
            closest_index=i
            min_diff=min_diff_here

    try:
        return_file=filelist[closest_index]
        return_xmax=xmax_list[closest_index]
    except:
        return_file='null'
        return_xmax=0

    return_dir=(return_file.split('_')[0]).split('SIM')[0]
    return_sim=(return_file.split('_')[0]).split('SIM')[1]

    return return_dir, return_sim,return_xmax



event_file=open('energy_events.txt','r')
event_list=np.genfromtxt(event_file,dtype='int')
event_file.close()

nEvents=len(event_list)

############# variables that i want to save ################




core_x_list=[]
type_list=[]

core_y_list=[]
events_list=[]
p_ratio_list=[]
d_ratio_list=[]
combchi2_list=[]
radiochi2_list=[]

clip_ratio_list=[]
charge_excess_ratio_list=[]
S_basic_list=[]
Srd_1_list=[]
Srd_2_list=[]
StoE_list=[]
rho_list=[]
a_list=[]

em_energy_list=[]
other_energy_list=[]
total_energy_list=[]
zenith_list=[]
azimuth_list=[]
energy_list=[]
xmax_list=[]
dmax_list=[]

xmax_fit_list=[]

alpha_list=[]
Erad_list=[]
clip_ratio_list=[]
charge_excess_ratio_list=[]
S_basic_list=[]
Srd_1_list=[]
Srd_2_list=[]
StoE_list=[]
rho_list=[]
a_list=[]








#for i in np.arange(nEvents):
for i in np.arange(e,e+1):

    print '______________________'
    print event_list[i]

    my_file3 = reco_dir+'reco'+str(int(event_list[i]))+'_3.dat'
    my_file2 = reco_dir+'reco'+str(int(event_list[i]))+'_2.dat'
    my_file1 = reco_dir+'reco'+str(int(event_list[i]))+'_1.dat'
    my_file = reco_dir+'reco'+str(int(event_list[i]))+'.dat'



    if os.path.exists(my_file3)==True:
        use_file=my_file3
    elif os.path.exists(my_file2)==True:
        use_file=my_file2
    elif os.path.exists(my_file1)==True:
        use_file=my_file1
    elif os.path.exists(my_file)==True:
        use_file=my_file
    else:
        print 'no file for {0}'.format(int(event_list[i]))
        continue


#try:
    for g in np.arange(1):
        analysisinfo = cPickle.load(open(use_file,'r'))
        #print analysisinfo.keys()


        xmaxreco_temp=analysisinfo['xmaxreco']
        xmax_temp=analysisinfo['realxmax']
        xmaxbest_temp=analysisinfo['xmaxbest']

        core_x_temp=analysisinfo['core_x']+analysisinfo['xoff']
        core_y_temp=analysisinfo['core_y']+analysisinfo['xoff']
        p_ratio_temp=analysisinfo['p_ratio']
        d_ratio_temp=analysisinfo['d_ratio']
        combchi2_temp=analysisinfo['combchi2']
        radiochi2_temp=analysisinfo['radiochi2']
        energy_temp=analysisinfo['energy']
        zenith_temp=analysisinfo['zenith']
        azimuth_temp=analysisinfo['azimuth']








        sim_dir,runnr,closest_xmax=find_SIM_dir(int(event_list[i]),xmaxreco_temp)
      
        event= sim_dir.split('events/')[1].split('/')[0]
        
        if 'proton' in sim_dir:
            type=0
        elif 'iron' in sim_dir:
            type=1
        else:
            type=2

        print 'type: {0}'.format(type)

        sim_zenith, sim_azimuth, alpha,sim_energy, hillas, antenna_positions, ant_pos_uvw,fluence11,fluence21,fluence41=sim.ProcessSim(sim_dir, runnr)# fluence in J/m2
        xmax=hillas[2]
        em_energy_hold,other_energy_hold,total_energy_hold=sim.getEM(sim_dir, runnr)
        atm=sim.get_atm(event)

        pos_uvw_vxb=ant_pos_uvw[0::8]
        pos_uvw_vxvxb=ant_pos_uvw[2::8]
        neg_uvw_vxvxb=ant_pos_uvw[6::8]


        fluence_vxb=fluence21[0::8]
        fluence_pos_vxvxb=fluence21[2::8]
        fluence_neg_vxvxb=fluence21[6::8]

        pos_vxvxb_all=np.concatenate([neg_uvw_vxvxb.T[1],pos_uvw_vxvxb.T[1]])
        fluence_vxvxb_0=np.concatenate([fluence_pos_vxvxb.T[0],fluence_neg_vxvxb.T[0]])
        fluence_vxvxb_1=np.concatenate([fluence_pos_vxvxb.T[1],fluence_neg_vxvxb.T[1]])


        inds = pos_vxvxb_all.argsort()

        sorted_pos=pos_vxvxb_all[inds]
        sorted_fluence_vxb=fluence_vxvxb_0[inds]
        sorted_fluence_vxvxb=fluence_vxvxb_1[inds]

        f0 = interp1d(sorted_pos, sorted_fluence_vxb, kind='cubic')
        f1 = interp1d(sorted_pos, sorted_fluence_vxvxb, kind='cubic')

        xnew = np.linspace(0, 400, num=1000, endpoint=True)

        Erad=integrate(xnew,f0(xnew),f1(xnew))




        hi=helper.get_vertical_height(xmax,atm)   # height in cm
    
        at=helper.get_atmosphere(hi,atm)
        rho=helper.get_density(hi, atm, allow_negative_heights=True)
        dmax=helper.get_distance_xmax_geometric(sim_zenith, xmax, atm)
    
        clip_ratio=rad.get_clipping(dmax)
        
        a=rad.get_a(rho)
        #print '_________________________________________'
        #print 'xmax: {0}'.format(xmax)
        #print 'zenith: {0}'.format(zenith*180/np.pi)
        
        #print 'relative charge excess: {0}'.format(a)
        
        S=rad.get_S_basic(Erad, np.sin(alpha))
        Srd=rad.get_S(Erad, np.sin(alpha), rho)
        Srd_final=rad.get_radiation_energy(Erad, np.sin(alpha), rho)
        
        StoE=rad.StoEm(Srd_final)



        xmax_fit_list.append(xmaxreco_temp)
        xmax_list.append(xmax)
        energy_list.append(sim_energy*1e9)
        core_x_list.append(core_x_temp)
        core_y_list.append(core_y_temp)
        p_ratio_list.append(p_ratio_temp)
        d_ratio_list.append(d_ratio_temp)
        combchi2_list.append(combchi2_temp)
        radiochi2_list.append(radiochi2_temp)
        events_list.append(event_list[i])
        zenith_list.append(sim_zenith)
        azimuth_list.append(sim_azimuth)
        Erad_list.append(Erad)
        alpha_list.append(alpha)
        em_energy_list.append(em_energy_hold*1e9)
        other_energy_list.append(other_energy_hold*1e9)
        total_energy_list.append(total_energy_hold*1e9)
        clip_ratio_list.append(clip_ratio)
        charge_excess_ratio_list.append(a)
        dmax_list.append(dmax)
        type_list.append(type)

        S_basic_list.append(S)
        Srd_1_list.append(Srd)
        Srd_2_list.append(Srd_final)
        rho_list.append(rho)

#   except:
#       print'can\'t open file'


info={'em_energy':np.asarray(em_energy_list),
    'other_energy':np.asarray(other_energy_list),
    'total_energy':np.asarray(total_energy_list),
    'zenith':np.asarray(zenith_list),
    'azimuth':np.asarray(azimuth_list),
    'energy':np.asarray(energy_list),
    'xmax':np.asarray(xmax_list),
    'xmax_fit':np.asarray(xmax_fit_list),
    'core_x':np.asarray(core_x_list),
    'core_y':np.asarray(core_y_list),
    'p_ratio':np.asarray(p_ratio_list),
    'd_ratio':np.asarray(d_ratio_list),
    'combchi2':np.asarray(combchi2_list),
    'radiochi2':np.asarray(radiochi2_list),
    'event':np.asarray(events_list),
    'alpha':np.asarray(alpha_list),
    'Erad':np.asarray(Erad_list),
    'clip':np.asarray(clip_ratio_list),
    'charge_excess_ratio':np.asarray(charge_excess_ratio_list),
    'S_basic':np.asarray(S_basic_list),
    'Srd_1':np.asarray(Srd_1_list),
    'Srd_2':np.asarray(Srd_2_list),
    'StoEM':np.asarray(StoE_list),
    'density':np.asarray(rho_list),
    'dmax':np.asarray(dmax_list),
    'type':np.asarray(type_list),

}

outfile=open('output/event_'+str(int(event_list[e]))+'_information.dat','w')
cPickle.dump(info,outfile)
outfile.close()






