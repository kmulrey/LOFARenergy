import numpy as np
from optparse import OptionParser
import os
import cPickle
import os.path


eventlist=np.genfromtxt(open('sim_event_list.txt'))
nEvents=len(eventlist)

print 'n starting events: {0}'.format(nEvents)


#reco_dir='/vol/astro3/lofar/sim/pipeline/run/analysis/'
reco_dir='/vol/astro7/lofar/sim/pipeline/production_analysis_radio_only/'

xmax=[]
energy=[]
core_x=[]
core_y=[]
event=[]
p_ratio=[]
d_ratio=[]
combchi2=[]
radiochi2=[]
event=[]
zenith=[]
azimuth=[]

for i in np.arange(nEvents):

    my_file3 = reco_dir+'reco'+str(int(eventlist[i]))+'_3.dat'
    my_file2 = reco_dir+'reco'+str(int(eventlist[i]))+'_2.dat'
    my_file1 = reco_dir+'reco'+str(int(eventlist[i]))+'_1.dat'
    my_file = reco_dir+'reco'+str(int(eventlist[i]))+'.dat'


    if os.path.exists(my_file3)==True:
        use_file=my_file3
    elif os.path.exists(my_file2)==True:
        use_file=my_file2
    elif os.path.exists(my_file1)==True:
        use_file=my_file1
    elif os.path.exists(my_file)==True:
        use_file=my_file
    else:
        print 'no file for {0}'.format(int(eventlist[i]))
        continue

    try:
        #for g in np.arange(1):
        analysisinfo = cPickle.load(open(use_file,'r'))
        #print analysisinfo.keys()
        xmax_temp=analysisinfo['xmaxreco']
        core_x_temp=analysisinfo['core_x']+analysisinfo['xoff']
        core_y_temp=analysisinfo['core_y']+analysisinfo['xoff']
        p_ratio_temp=analysisinfo['p_ratio']
        d_ratio_temp=analysisinfo['d_ratio']
        combchi2_temp=analysisinfo['combchi2']
        radiochi2_temp=analysisinfo['radiochi2']
        energy_temp=analysisinfo['energy']
        zenith_temp=analysisinfo['zenith']
        azimuth_temp=analysisinfo['azimuth']
        
        event.append(int(eventlist[i]))
        xmax.append(xmax_temp)
        core_x.append(core_x_temp)
        core_y.append(core_y_temp)
        combchi2.append(combchi2_temp)
        radiochi2.append(radiochi2_temp)
        energy.append(energy_temp)
        zenith.append(zenith_temp)
        azimuth.append(azimuth_temp)
        p_ratio.append(p_ratio_temp)
        d_ratio.append(d_ratio_temp)

    except:
        print 'missing data for {0}'.format(int(eventlist[i]))
        continue


outfile=open('CR_event_info.txt','w')

for i in np.arange(len(event)):
    outfile.write('{0} {1} {2} {3} {4} {5} {6} {7} {8}\n'.format(int(event[i]),xmax[i],core_x[i],core_y[i],energy[i],zenith[i],azimuth[i],p_ratio[i],d_ratio[i]))

outfile.close()

'''
info={'event_number': np.asarray(event), 'xmax':np.asarray(xmax), 'core_x':np.asarray(core_x), 'core_y':np.asarray(core_y), 'energy':np.asarray(energy), 'zenith':np.asarray(zenith), 'azimuth':np.asarray(azimuth), 'p_ratio':np.asarray(p_ratio), 'd_ratio':np.asarray(d_ratio), 'radiochi2': np.asarray(radiochi2), 'combchi2': np.asarray(combchi2)}
        
fout=open('fit_info.p','w')
            
cPickle.dump(info,fout)
fout.close()

'''