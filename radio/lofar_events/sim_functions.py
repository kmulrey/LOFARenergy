import numpy as np
from optparse import OptionParser
import cPickle
import re
from scipy.signal import hilbert
from scipy.signal import resample
import scipy.fftpack as fftp
import os
from scipy import signal
from scipy.optimize import curve_fit


atm_dir='/Users/kmulrey/radio/atmosphere_files/'
#atm_dir='/vol/astro7/lofar/sim/pipeline/atmosphere_files/'

conversion_factor_integrated_signal = 2.65441729e-3 * 6.24150934e18  # to convert V**2/m**2 * s -> J/m**2 -> eV/m**2

mag=2.04
#import pycrtools as cr
#import process_func as prf

def GetUVW(pos, cx, cy, cz, zen, az, Binc):
    relpos = pos-np.array([cx,cy,cz])
    B = np.array([0,np.cos(Binc),-np.sin(Binc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    return np.array([np.inner(vxB,relpos),np.inner(vxvxB,relpos),np.inner(v,relpos)]).T

def GetAlpha(zen,az,Binc):
    B = np.array([0,np.cos(Binc),-np.sin(Binc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    return np.arccos(np.inner(np.asarray(B), np.asarray(v)) / (np.linalg.norm(B) * np.linalg.norm(v)))




def GetXYZ(pos, zen, az):
    inc=1.1837
    B = np.array([0,np.cos(inc),-np.sin(inc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    #print v
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    return pos[0]*vxB+pos[1]*vxvxB+pos[2]*v
'''
def GetUVW(pos, cx, cy, zen, az):
    relpos = pos-np.array([cx,cy,7.6])
    inc=1.1837
    B = np.array([0,np.cos(inc),-np.sin(inc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    return np.array([np.inner(vxB,relpos),np.inner(vxvxB,relpos),np.inner(v,relpos)]).T
'''
def GetUVW_efield(efield, cx, cy, zen, az):
    #relpos = pos-np.array([cx,cy,7.6])
    inc=1.1837
    B = np.array([0,np.cos(inc),-np.sin(inc)])
    v = np.array([-np.cos(az)*np.sin(zen),-np.sin(az)*np.sin(zen),-np.cos(zen)])
    vxB = np.array([v[1]*B[2]-v[2]*B[1],v[2]*B[0]-v[0]*B[2],v[0]*B[1]-v[1]*B[0]])
    vxB = vxB/np.linalg.norm(vxB)
    vxvxB = np.array([v[1]*vxB[2]-v[2]*vxB[1],v[2]*vxB[0]-v[0]*vxB[2],v[0]*vxB[1]-v[1]*vxB[0]])
    print '_______________ in getUVW_efield_______________'
    print vxB.shape
    print vxvxB.shape
    print efield.shape
#return np.array([np.inner(vxB,relpos),np.inner(vxvxB,relpos),np.inner(v,relpos)]).T
    return 0


def ProcessSim(datadir,fileno):
    
    lSample=128*2
    lFFT=lSample/2+1
    
    lFFTkeep=20
    
    
    longfile = '{0}/DAT{1}.long'.format(datadir,str(fileno).zfill(6))
    steerfile = '{0}/steering/RUN{1}.inp'.format(datadir,str(fileno).zfill(6))
    listfile = open('{0}/steering/SIM{1}.list'.format(datadir,str(fileno).zfill(6)))
    lines = listfile.readlines()
    nTotalAnt=len(lines)
    
    antenna_positions=np.zeros([0,3])
    antenna_files=[]
    
    for l in np.arange(nTotalAnt):
        antenna_position_hold=np.asarray([float(lines[l].split(" ")[2]),float(lines[l].split(" ")[3]),float(lines[l].split(" ")[4])])#read antenna position...
        antenna_file_hold=(lines[l].split(" ")[5].split()[0])   #... and output filename from the antenna list file
        antenna_files.append(antenna_file_hold)
        antenna_positions=np.concatenate((antenna_positions,[antenna_position_hold]))

    lowco=30.0
    hico=80.0
    nantennas=len(antenna_files)
   
    onskypower=np.zeros([nantennas,2])
    #antenna_position=np.zeros([nantennas,3])
    filteredpower=np.zeros([nantennas,2])
    power=np.zeros([nantennas,2])
    power11=np.zeros([nantennas,2])
    power21=np.zeros([nantennas,2])
    power41=np.zeros([nantennas,2])
    peak_time=np.zeros([nantennas,2])
    peak_bin=np.zeros([nantennas,2])
    peak_amplitude=np.zeros([nantennas,2])
    pol_angle=np.zeros([nantennas])
    pol_angle_filt=np.zeros([nantennas])


    hillas = np.genfromtxt(re.findall("PARAMETERS.*",open(longfile,'r').read()))[2:]
    zenith=(np.genfromtxt(re.findall("THETAP.*",open(steerfile,'r').read()))[1])*np.pi/180. #rad; CORSIKA coordinates
    azimuth=np.mod(np.genfromtxt(re.findall("PHIP.*",open(steerfile,'r').read()))[1],360)*np.pi/180.  #rad; CORSIKA coordinates
    energy=np.genfromtxt(re.findall("ERANGE.*",open(steerfile,'r').read()))[1] #GeV


    fullFFT=np.zeros([nantennas,lFFTkeep])
    fullPower=np.zeros([nantennas])
    fullFrequencies=np.zeros([lFFTkeep])
    wfX=np.zeros([nantennas,81])
    wfY=np.zeros([nantennas,81])
    wfZ=np.zeros([nantennas,81])
    time_all=np.zeros([81])
    wfX=np.zeros([nantennas,81])
    wfY=np.zeros([nantennas,81])
    wfZ=np.zeros([nantennas,81])

    for j in np.arange(nantennas):
        #for j in np.arange(1):

        antenna_file = lines[j].split(" ")[5]
        coreasfile = '{0}/SIM{1}_coreas/raw_{2}.dat'.format(datadir,str(fileno).zfill(6),antenna_files[j])

        data=np.genfromtxt(coreasfile)
        data[:,1:]*=2.99792458e4 # convert Ex, Ey and Ez (not time!) to Volt/meter
        dlength=data.shape[0]
        poldata=np.ndarray([dlength,2])
        XYZdata=np.ndarray([dlength,2])
        az_rot=3*np.pi/2+azimuth    #conversion from CORSIKA coordinates to 0=east, pi/2=north
        zen_rot=zenith
        XYZ=np.zeros([dlength,3])
        XYZ[:,0]=-data[:,2] #conversion from CORSIKA coordinates to 0=east, pi/2=north
        XYZ[:,1]=data[:,1]
        XYZ[:,2]=data[:,3]
        
        #x_data=XYZ.T[0]
        #y_data=XYZ.T[1]
        #z_data=XYZ.T[2]
        
        

        # Convert to, v, vxB, vxvxB coordinates to compute Stokes parameters and polarization angle
        UVW=GetUVW(XYZ,0,0,0,zen_rot,az_rot,1.1837)
        alpha= GetAlpha(zen_rot,az_rot,1.1837)
        #Stokes=prf.stokes_parameters(UVW[:,0],UVW[:,1],fftp.hilbert(UVW[:,0]),fftp.hilbert(UVW[:,1]))
        #pol_angle[j]=prf.polarization_angle(Stokes)
        #UVWfilt=prf.FreqFilter(UVW,30,80,data[1,0]-data[0,0])
        #Stokesfilt=prf.stokes_parameters(UVWfilt[:,0],UVWfilt[:,1],fftp.hilbert(UVWfilt[:,0]),fftp.hilbert(UVWfilt[:,1]))
        
        #pol_angle_filt[j]=prf.polarization_angle(Stokesfilt)
        # Convert to on-sky coordinates (n, theta, phi) to prepare for application of antenna model
        #poldata[:,0] = -1.0/np.sin(zen_rot)*data[:,3] # -1/sin(theta) *z
        #poldata[:,1] = np.sin(az_rot)*data[:,2] + np.cos(az_rot)*data[:,1] # -sin(phi) *x + cos(phi)*y in coREAS 0=positive y, 1=negative x
        #poldata[:,0] = -1.0/np.sin(zen_rot)*XYZ[:,2] # -1/sin(theta) *z
        #poldata[:,1] = np.sin(az_rot)*(-1)*XYZ[:,0] + np.cos(az_rot)*XYZ[:,1] # -sin(phi) *x + cos(phi)*y in coREAS 0=positive y, 1=negative x
        
        poldata[:,0] = UVW[:,0]
        poldata[:,1] = UVW[:,1]

        spec=np.fft.rfft(poldata, axis=-2)
        #print spec.shape
        # Apply antenna model
        tstep = data[1,0]-data[0,0]
        onskypower[j]=np.array([np.sum(poldata[:,0]*poldata[:,0]),np.sum(poldata[:,1]*poldata[:,1])])*tstep


        freqhi = 0.5/tstep/1e6 # MHz
        freqstep = freqhi/(dlength/2+1) # MHz
        frequencies = np.arange(0,freqhi,freqstep)*1e6 # Hz
        frequencies = np.arange(0,dlength/2+1)*freqstep*1e6

        #Apply window and reduce maximum frequency to acquire downsampled signal
        fb = int(np.floor(lowco/freqstep))
        lb = int(np.floor(hico/freqstep)+1)
        window = np.zeros([1,dlength/2+1,1])
        window[0,fb:lb+1,0]=1
        
        pow0=np.abs(spec[:,0])*np.abs(spec[:,0])
        pow1=np.abs(spec[:,1])*np.abs(spec[:,1])
        
        
        ospow0=np.abs(spec[:,0])*np.abs(spec[:,0])
        ospow1=np.abs(spec[:,1])*np.abs(spec[:,1])
        power[j]=np.array([np.sum(pow0[fb:lb+1]),np.sum(pow1[fb:lb+1])])/(dlength/2.)*tstep
        filteredpower[j]=np.array([np.sum(ospow0[fb:lb+1]),np.sum(ospow1[fb:lb+1])])/(dlength/2.)*tstep
        
        # assume that simulated time resolution is higher than LOFAR time resolution (t_step=5 ns)
        maxfreqbin= int(np.floor(tstep/5e-9 * dlength/2.)+1)
        shortspec=np.array([spec[0:maxfreqbin,0]*window[0,0:maxfreqbin,0],spec[0:maxfreqbin,1]*window[0,0:maxfreqbin,0]])
        filt=np.fft.irfft(shortspec, axis=-1)
        
        # after downsampling, renormalize the signal!
        dlength_new=filt.shape[1]
        filt=filt*1.0*dlength_new/dlength
        # to calculate the time of arrival upsample with a factor 5
        filt_upsampled=resample(filt,5*dlength_new,axis=-1)
        # compute hilbert enevelope
        hilbenv=np.abs(hilbert(filt,axis=-1))
        hilbenv_upsampled=np.abs(hilbert(filt_upsampled,axis=-1))
        
        # peak_time is the bin where the maximum is located; NOT the actual time of the peak!
        peak_bin[j]=np.argmax(hilbenv,axis=-1)
        peak_time[j]=np.argmax(hilbenv_upsampled,axis=-1)*1e-9 #in seconds
        peak_amplitude[j]=np.max(hilbenv_upsampled,axis=-1)
        if (peak_amplitude[j,0]>peak_amplitude[j,1]):
            pt=peak_bin[j,0]
        else:
            pt=peak_bin[j,1]
        
        # for 3 different window size, the total power is calculated. The window is allowed to `wrap around', so some voodoo is needed to determine the range:
        d=filt.shape[1]
        rng=5
        a=int(np.max([0,pt-rng]))
        b=int(pt+rng+1)
        c=int(np.min([d,pt+d-rng]))
        power11[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))*5e-9
        rng=10
        a=int(np.max([0,pt-rng]))
        b=int(pt+rng+1)
        c=int(np.min([d,pt+d-rng]))
        power21[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))*5e-9
        rng=20
        a=int(np.max([0,pt-rng]))
        b=int(pt+rng+1)
        c=int(np.min([d,pt+d-rng]))
        power41[j]=(np.sum(np.square(filt[:,a:b]),axis=-1)+np.sum(np.square(filt[:,c:d]),axis=-1))*5e-9

    temp=np.copy(antenna_positions)
    antenna_positions[:,0], antenna_positions[:,1], antenna_positions[:,2] = -1*(temp[:,1])/100.,(temp[:,0])/100., temp[:,2]/100.

    azimuth=3*np.pi/2+azimuth #  +x = east (phi=0), +y = north (phi=90)
    
    ant_pos_uvw=GetUVW(antenna_positions, 0, 0, 0, zenith, az_rot,1.1837)

    return zenith, azimuth, alpha, energy, hillas, antenna_positions, ant_pos_uvw, power11/(377.0),power21/(377.0),power41/(377.0)
#return zenith, azimuth, energy, hillas, antenna_positions, ant_pos_uvw, fullPower, fullFFT, fullFrequencies,  wfX,  wfY,  wfZ, time_all





# read coreas and return sim information and traces in V/m

def get_efield(datadir,fileno):

    dlength=4082
    longfile = '{0}/DAT{1}.long'.format(datadir,str(fileno).zfill(6))
    steerfile = '{0}/steering/RUN{1}.inp'.format(datadir,str(fileno).zfill(6))
    listfile = open('{0}/steering/SIM{1}.list'.format(datadir,str(fileno).zfill(6)))
    lines = listfile.readlines()
    nTotalAnt=len(lines)
    
    antenna_positions=np.zeros([0,3])
    antenna_files=[]
    
    efield=np.zeros([nTotalAnt,dlength,3])
    time=np.zeros([nTotalAnt,dlength])
    
    for l in np.arange(nTotalAnt):
        antenna_position_hold=np.asarray([float(lines[l].split(" ")[2]),float(lines[l].split(" ")[3]),float(lines[l].split(" ")[4])])#read antenna position...
        antenna_file_hold=(lines[l].split(" ")[5].split()[0])   #... and output filename from the antenna list file
        antenna_files.append(antenna_file_hold)
        antenna_positions=np.concatenate((antenna_positions,[antenna_position_hold]))


    nantennas=len(antenna_files)
    

    hillas = np.genfromtxt(re.findall("PARAMETERS.*",open(longfile,'r').read()))[2:]
    zenith=(np.genfromtxt(re.findall("THETAP.*",open(steerfile,'r').read()))[1])*np.pi/180. #rad; CORSIKA coordinates
    azimuth=np.mod(np.genfromtxt(re.findall("PHIP.*",open(steerfile,'r').read()))[1],360)*np.pi/180.  #rad; CORSIKA coordinates
    az_rot=3*np.pi/2+azimuth    #conversion from CORSIKA coordinates to 0=east, pi/2=north

    energy=np.genfromtxt(re.findall("ERANGE.*",open(steerfile,'r').read()))[1] #GeV

    for j in np.arange(nantennas):
            #for j in np.arange(1):
            
        antenna_file = lines[j].split(" ")[5]
        coreasfile = '{0}/SIM{1}_coreas/raw_{2}.dat'.format(datadir,str(fileno).zfill(6),antenna_files[j])

        data=np.genfromtxt(coreasfile)
        data[:,1:]*=2.99792458e4 # convert Ex, Ey and Ez (not time!) to Volt/meter
        dlength=data.shape[0]

        XYZ=np.zeros([dlength,3])
        XYZ[:,0]=-data[:,2] #conversion from CORSIKA coordinates to 0=east, pi/2=north
        XYZ[:,1]=data[:,1]
        XYZ[:,2]=data[:,3]
        
        UVW=GetUVW(XYZ,0,0,0,zenith,az_rot,1.1837)
        #poldata[:,0] = UVW[:,0]
        #poldata[:,1] = UVW[:,1]

        
        
        
        efield[j]=UVW#data[:,1:]#UVW#
        time[j]=data.T[0]
    
    temp=np.copy(antenna_positions)
    antenna_positions[:,0], antenna_positions[:,1], antenna_positions[:,2] = -1*(temp[:,1])/100.,(temp[:,0])/100., temp[:,2]/100.
    
    return antenna_positions,time,efield,zenith,az_rot,energy,hillas[2]


def ground_to_shower(array,zenith,az,B=1.1837):
    
    array_uvw=GetUVW(array, 0, 0, 0, zenith, az,B)
    return array_uvw


def filter(times,traces, fmin, fmax): # 2d-> time,data
    if traces.ndim == 2:
        traces = np.expand_dims(traces, axis=0)

    filetered_traces = np.full_like(traces, 0)

    nTraces=len(traces)
    tstep=times[0][1]-times[0][0]
    ldata=len(traces[1][0])

    rate=1/tstep
    nyq=0.5*rate
    lowcut=fmin*1e6
    highcut=fmax*1e6
    low = lowcut/nyq
    high = highcut/nyq
    order=3
    b, a = signal.butter(order, [low, high], btype='band')

    for i in np.arange(nTraces):
        filetered_traces[i].T[0]= signal.filtfilt(b, a, traces[i].T[0])  # this is data in the time domain
        filetered_traces[i].T[1]= signal.filtfilt(b, a, traces[i].T[1])  # this is data in the time domain
        filetered_traces[i].T[2]= signal.filtfilt(b, a, traces[i].T[2])  # this is data in the time domain

    return np.squeeze(filetered_traces)

def get_atm(event):
    
    # read atmospheric parameters ATMLAY, A, B, C respectively for event

    filename=atm_dir+'ATMOSPHERE_'+event+'.DAT'
    file=open(filename,'r')
    atm=np.genfromtxt(file,skip_header=1,max_rows=4)
    file.close()


    return atm




def getEM(datadir,fileno):
    
    longfile = '{0}/DAT{1}.long'.format(datadir,str(fileno).zfill(6))
    lookup1='LONGITUDINAL ENERGY DEPOSIT'
    lookup2='FIT OF THE HILLAS CURVE'
    start=-1
    stop=-1
    
    with open(longfile) as myFile:
        for num, line in enumerate(myFile, 1):
            if lookup1 in line:
                start=num+2
            if lookup2 in line:
                stop=num-1

    myFile.close()
    myFile=open(longfile,'r')
    longinfo=np.genfromtxt(myFile,skip_header=start-1,max_rows=(stop-start+1))

    myFile.close()

    
    total_dep= np.sum(longinfo.T[9])
    em_dep=np.sum(longinfo.T[1]+longinfo.T[2]+longinfo.T[3])
    other_dep=np.sum(longinfo.T[4]+longinfo.T[5]+longinfo.T[6]+longinfo.T[7]+longinfo.T[8])
    
    return em_dep,other_dep,total_dep

def test_lin(x, a,b):
    return b*x+a+7


def getFit(x,y):
    
    y_use=np.log10(y)
    x_use=np.log10(x/1.0e18)
    
    
    x_use=x_use[(y_use>-10)*(y_use<10)]
    y_use=y_use[(y_use>-10)*(y_use<10)]
    
    param1,param_cov1=curve_fit(test_lin,x_use,y_use,p0=[1,1])
    a=param1[0]
    b=param1[1]
    perr=np.sqrt(np.diag(param_cov1))
    a_err=np.power(10,a+perr[0])-np.power(10,a)
    
    
    return np.power(10,a),b,a_err,perr[1]

def energy_to_rad(x,a,b):
    return a*1e7*(x/1e18)**b#*(mag)**b

