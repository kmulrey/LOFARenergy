import numpy as np
import cPickle as pickle
from optparse import OptionParser
import scipy.fftpack as fftp
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import cmath
import math
import scipy.interpolate as intp
from scipy.interpolate import interp1d
from scipy import signal
from matplotlib.colors import LogNorm
plt.ion()


def theta_phi( THETA, PHI, PSI, X0, Y0, Z0):
    # to put ground plane into shower plane
    #1st ROTATION: counterclockwise from X-axis(N) for angle 'phi' about Z-axis.
    X1= X0*np.cos(PHI)+Y0*np.sin(PHI)
    Y1=-X0*np.sin(PHI)+Y0*np.cos(PHI)
    Z1= Z0
    #-------------xxx----------------
    #2nd ROTATION: clockwise from Z-axis for angle 'theta' about Y-axis(W).
    X2= X1*np.cos(THETA)-Z1*np.sin(THETA)
    Y2= Y1
    Z2= X1*np.sin(THETA)+Z1*np.cos(THETA)
    #-------------xxx----------------
    #3rd ROTATION: counterclockwise from X-axis(N) for angle 'psi' about Z-axis.
    a= X2*np.cos(PSI)+Y2*np.sin(PSI)
    b=-X2*np.sin(PSI)+Y2*np.cos(PSI)
    c= Z2
    return a,b,c
#-------------xxx----------------


def back_theta_phi(THETA, PHI, PSI, X0, Y0, Z0):
#1st ROTATION:clockwise from X-axis(N) for angle 'psi' about Z-axis.

    X1= X0*np.cos(PSI)-Y0*np.sin(PSI)
    Y1= X0*np.sin(PSI)+Y0*np.cos(PSI)
    Z1= Z0
    #-------------xxx------------------
    #2nd ROTATION: counterclockwise from Z-axis for angle 'theta' about Y-axis(W).
    X2= X1*np.cos(THETA)+Z1*np.sin(THETA)
    Y2= Y1
    Z2= -X1*np.sin(THETA)+Z1*np.cos(THETA)
    #--------------xxx-------------------
    #3rd ROTATION:clockwise from X-axis(N) for angle 'phi' about Z-axis.
    a= X2*np.cos(PHI)-Y2*np.sin(PHI)
    b= X2*np.sin(PHI)+Y2*np.cos(PHI)
    c= Z2
    return a,b,c

#--------------xxx------------------



file='geant/proton/LORA000059.geant'
set='20'
dir='/Users/kmulrey/LOFAR/energy/LOFARenergy/metadata/'
#info=deposits = np.genfromtxt(file,max_rows=1)
info = np.genfromtxt(file,max_rows=1)
deposits = np.genfromtxt(file,skip_header=1)


file1=open(dir+'Detector_Cord_All.dat')
file2=open(dir+'LOFARstation_centers.txt')
infile=open(dir+'detectors'+set+'.txt')

LORAdetectors=np.genfromtxt(file1,comments="//",usecols=(1,2,3))
LOFARstations=np.genfromtxt(file2,comments="//", usecols=(1,2,3))
detector_sets=np.genfromtxt(infile)

plot_det=LORAdetectors[detector_sets==1]
nDetTotal=len(LORAdetectors)
nDetPlot=len(plot_det)


energy = info[1]
zenith = info[2]
phi = info[3]
psi=(2*np.pi-phi)

print 'energy: {0:.2f}'.format(np.log10(energy)+9)
print 'theta: {0:.2f}'.format(zenith*180/np.pi)
print 'phi: {0:.2f}'.format(phi*180/np.pi)

core_x = 50
core_y = -150


nRad=200
nAng=24

nPos=nRad*nAng


xy=np.zeros([3,nRad*nAng])
z=np.zeros([nRad*nAng])
count=0

closest_rad=1000
closest_dep=1000

muon_dep=6.3


for i in np.arange(nAng):
    for j in np.arange(nRad):
        theta = i*(360.0/nAng)*(np.pi/180.0)
        r=0.5+5.0*j
        xy[0][count] = np.cos(theta)*r + core_x
        xy[1][count] = np.sin(theta)*r + core_y
        xy[2][count] = 0

        z[count] = deposits[j][1]
        count =count+1
        if i==0:
            if  np.abs(deposits[j][1]-muon_dep)<np.abs(closest_dep-muon_dep):
                closest_dep=deposits[j][1]
                closest_rad=r=0.5+5.0*j
detection_ring=np.zeros([3,3600])
detection_ring_plane=np.zeros([3,3600])

for i in np.arange(3600):
    th=360.0/3600.0*i
    detection_ring[0][i]=closest_rad*np.cos(th)+core_x
    detection_ring[1][i]=closest_rad*np.sin(th)+core_y
    detection_ring[2][i]=0.0


xy_plane=np.zeros([3,nRad*nAng])
LORA_plane=np.zeros([3,len(plot_det)])
LOFAR_plane=np.zeros([3,len(LOFARstations)])


xy_plane[0],xy_plane[1],xy_plane[2] = back_theta_phi(zenith, phi, psi, xy[0], xy[1], xy[2])
LORA_plane[0],LORA_plane[1],LORA_plane[2] = theta_phi(zenith, phi, psi, plot_det.T[0],plot_det.T[1],plot_det.T[2])
detection_ring_plane[0],detection_ring_plane[1],detection_ring_plane[2] = back_theta_phi(zenith, phi, psi, detection_ring[0], detection_ring[1], detection_ring[2])

#for i in np.arange(len(LORA_plane.T)):

#print '{0}    {1}'.format(plot_det.T[0][i], LORA_plane[0][i])

rbf0= np.ndarray([nPos],dtype=object)
rbf0 = intp.Rbf(xy_plane[0], xy_plane[1], z, function='linear')
#rbf0 = intp.Rbf(pos_sim_UVW.T[0], pos_sim_UVW.T[1], SNR.T[0],smooth =0)#,function='quintic')

dist_scale=1000.0
ti = np.linspace(-dist_scale, dist_scale, 100)

XI, YI = np.meshgrid(ti, ti)
ZI0 = rbf0(XI, YI)
maxp = np.max([np.max(ZI0)])

circle1 = plt.Circle((0, 0), closest_rad, color='r', fill=False,lw=3)

fig = plt.figure(facecolor='white')
ax = plt.gca()

plt.pcolor(XI, YI, ZI0,vmax=1000, vmin=0.5,cmap=cm.jet,norm = LogNorm())#,norm = LogNorm()
ax.set_facecolor(cm.jet(0))

#plt.plot(LOFARstations.T[0],LOFARstations.T[1],'.',markersize=8, mew=2, marker='x',color='white')
plt.plot(plot_det.T[0],plot_det.T[1],'.',markersize=4,marker='s',color='white')
#plt.plot(LORA_plane[0],LORA_plane[1],'.',markersize=7,marker='s',color='white')
#plt.plot(detection_ring_plane[0],detection_ring_plane[1],'.',markersize=2,color='red')

#plt.gcf().gca().add_artist(circle1)
plt.colorbar(label='energy deposited MeV/m$^2$')

plt.axes().set_aspect('equal')

plt.xlim((-1200,1200))
plt.ylim((-1200,1200))

ax.set_xlabel('X (east, m)',fontsize=12)
ax.set_ylabel('Y (north, m)',fontsize=12)


#plt.xlim((-800,800))
#plt.ylim((-800,800))

plt.show()
raw_input()
plt.close()

'''
fig=plt.figure(facecolor='white',figsize=(8,8))

ax1=fig.add_subplot(1,1,1)
colors.Normalize(clip=False)

ax1.pcolor(XI, YI, ZI0,vmax=maxp, vmin=1e-1,cmap=cm.jet,norm = LogNorm())
#ax1.plot(xy[0],xy[1],'o')

ax1.axis([-1200, 1200, -1200, 1200])
plt.show()
'''
