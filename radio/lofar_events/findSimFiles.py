#from pycrtools import crdatabase as crdb
import numpy as np
from optparse import OptionParser
import os
import cPickle
import re




CR_events=np.genfromtxt(open('CR_event_info.txt','r'))
event_path='/vol/astro7/lofar/sim/pipeline/events/'

print CR_events.shape

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
        dirProton=event_path+'/0/coreas/proton/'
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

    return return_file,return_xmax




i=10

event_list=[]
file_list=[]

for i in np.arange(len(CR_events)):
    #for i in np.arange(3):

    #print '............{0}......{1}.............'.format(CR_events[i][0],CR_events[i][8])
    file,closest_xmax=find_SIM_dir(CR_events[i][0],CR_events[i][1])
    print '{0}  {1}'.format(file,closest_xmax)
    if closest_xmax!=0 and file!='null':
        file_list.append(file)
        event_list.append(int(CR_events[i][0]))
    #print '{0}  {1}'.format(file,closest_xmax)


outputfile=open('best_sim_files.txt','w')

for i in np.arange(len(event_list)):
    print '{0}  {1}'.format(event_list[i],file_list[i])
    directory=file_list[i].split('SIM')[0]
    sim_n=file_list[i].split('SIM')[1].split('_')[0]

    #print directory
    #print sim_n
    outputfile.write('{0}  {1}  {2}\n'.format(event_list[i],directory, sim_n))

outputfile.close()









