import os
import numpy as np
import cPickle as pickle
from optparse import OptionParser
import scipy.fftpack as fftp
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import signal
import ROOT


nTrace=4000
nDet=20   # total number of LORA detectors
nLasa=5   # number of LORA stations (5 stations of 4 scintillators)

data_dir='/Users/kmulrey/LOFAR/LORA/LORAdata/data/'




def readfile(filename,e):
    root_file=ROOT.TFile.Open(data_dir+filename)


    #ROOT files contain:

    #Tree_event: 20 branches, one for each scintillator.  Contains an enetry for  each event, including GPS time stamp, trigger condition, pulse height, width, total integrated ADC counts, and traces
    #Tree_sec: contains branches with information saved every second about each LORA station (4 scintillators) to monitor behavior.
    #Tree_log: information about settings, including constants for electronics calibration, trigger settings on the individual detector level
    #Tree_noise: information about average noise levels in each detector, including mean and sigma for ADC counts

    tree_event = root_file.Get("Tree_event")
    tree_sec = root_file.Get("Tree_sec")
    tree_log = root_file.Get("Tree_log")
    tree_noise = root_file.Get("Tree_noise")

# check to see how many events are saved
    det1=tree_event.GetBranch('Det1')
    nE= det1.GetEntries()
    #print 'number of events in this file: {0}'.format(nE)




# function to read root data
    def getData(det, entry):
    
        det.GetEntry(entry)
        detector=det.GetLeaf('detector').GetValue()
        ymd=det.GetLeaf('YMD').GetValue()
        gps=det.GetLeaf('GPS_time_stamp').GetValue()
        ctd=det.GetLeaf('CTD').GetValue()
        nsec=det.GetLeaf('nsec').GetValue()
        trigg_condition=det.GetLeaf('Trigg_condition').GetValue()
        try:
            trigg_pattern=det.GetLeaf('Trigg_pattern').GetValue()
        except:
            trigg_pattern=-1
        total_counts=det.GetLeaf('Total_counts').GetValue()
        pulse_height=det.GetLeaf('Pulse_height').GetValue()
        pulse_width=det.GetLeaf('Pulse_width').GetValue()
        counts=det.GetLeaf('counts')
        hold=np.zeros([nTrace])
        for i in np.arange(nTrace):
            hold[i]=counts.GetValue(i)
    
    
        return detector,ymd,gps,nsec,trigg_condition,total_counts,pulse_height,pulse_width,hold






            #e=1000  # pick an event to look at


    detector=np.zeros([nDet])
    ymd=np.zeros([nDet])   # 0 if no trigger
    gps=np.zeros([nDet])   # UTC time of event, 1 if no trigger
    nsec=np.zeros([nDet])  # ns timing info, 0 if no trigger
    trigg_condition=np.zeros([nDet])  # how many scintillators out of 4 for a LORA station trigger

    total_counts=np.zeros([nDet]) # integrated counts in time window (-70ns before peak to 300ns after peak), background subtracted
    pulse_height=np.zeros([nDet]) # peak ADC count (background subtrated)
    pulse_width=np.zeros([nDet]) # RMS of root histogram containing pulse
    counts=np.zeros([nDet,nTrace])  # time traces for each detector (delta t = 2.5 ns)

# if the LORA station (4 scintillators) has triggered, traces for each scintillator are saved,  if not, trace contains 0s



# loop over 20 detectors to get event info
    for d in np.arange(nDet):
        detname='Det'+str(d+1)
        det=tree_event.GetBranch(detname)
        
        detector[d],ymd[d],gps[d],nsec[d],trigg_condition[d],total_counts[d],pulse_height[d],pulse_width[d],counts[d]=getData(det,e)


    return counts



























