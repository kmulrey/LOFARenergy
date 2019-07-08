import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'simulation/')
import tools as tools

import fluence as flu
import radiation_energy as rad
import sim_functions as sim
import helper as helper

plt.ion()


mag=2.03
average_density = 6.5e-4 #in g/cm^3#atmc.get_density(average_zenith, average_xmax) * 1e-3  # in kg/m^3
p0=0.250524463912
p1=-2.95290494

p0_=0.239
p1_=-3.13

def read_file(filename):
    
    infile=open(filename,'r')
    info=cPickle.load(infile)
    infile.close()
    print info.keys()
    
    em_energy=info['em_energy']
    other_energy=info['other_energy']
    total_energy=info['total_energy']
    zenith=info['zenith']
    azimuth=info['azimuth']
    energy=info['energy']
    xmax=info['xmax']
    alpha=info['alpha']
    Erad=info['Erad']
    clip=info['clip']
    charge_excess_ratio=info['charge_excess_ratio']
    S_basic=info['S_basic']
    Srd_1=info['Srd_1']
    Srd_2=info['Srd_2']
    density=info['density']
    xmax_fit=info['xmax_fit']
    core_x=info['core_x']
    core_y=info['core_y']
    p_ratio=info['p_ratio']
    d_ratio=info['d_ratio']
    combchi2=info['combchi2']
    radiochi2=info['radiochi2']
    event=info['event']
    density=info['density']
    dmax=info['dmax']
    type=info['type']
    
    Erad=Erad+Erad*.11-Erad*0.0336
    S_basic=S_basic+S_basic*.11-S_basic*0.0336
    Srd_1=Srd_1+Srd_1*.11-Srd_1*0.0336
    Srd_2=(Srd_2+Srd_2*.11-Srd_2*0.0336)
    
    corr=((1 - p0 + p0 * np.exp(p1 * (density - average_density)*1e3)) )**2
    corr_=((1 - p0_ + p0_ * np.exp(p1_ * (density - average_density)*1e3)) )**2
    
    
    return em_energy,energy,zenith,azimuth,xmax,alpha,S_basic,Srd_1,Srd_2*corr/corr_,Erad,charge_excess_ratio,event,p_ratio,d_ratio,type


#em_energy_int,energy_int,zenith_int,azimuth_int,xmax_int,alpha_int,S_basic_int,Srd_1_int,Srd_2_int,Erad_int,a_int,event_int,p_ratio,d_ratio,type=read_file('lofar_events/compiled_sims_all.dat')


coreas_file='lofar_events/compiled_sims_all.dat'
ldf_file='lofar_events/parameterization_radiation_energy.p'
uncertainty_file='lofar_events/uncertainty_info.dat'


event,radiation_i,radiation_p,azimuth,zenith,xmax,em_energy,sim_energy,p_ratio,d_ratio,type,alpha_i,alpha_p,sigma_e,sigma_e_radio,lora_energy,ldf_particle_err,ldf_radio_err,S_basic,Srd_1,Srd_2=tools.combine_files(coreas_file,ldf_file,uncertainty_file)

### ldf already has div, sin2alpha

particle_scale=(6.3/6.7)#(1.1)*(6.1/6.7)

S_use=Srd_1*p_ratio
E_use=sim_energy*d_ratio*particle_scale

#S_use=S_use[E_use>7e16]
#E_use=E_use[E_use>7e16]

a_lofar,b_lofar,a_sig_lofar,b_sig_lofar=sim.getFit(E_use[(E_use>1e17)*(E_use<1e18)],S_use[(E_use>1e17)*(E_use<1e18)])
print '\n\nEM energy fit'
print '___________________'
print '{0:.3f} +/- {1:.3f}       {2:.3f} +/- {3:.3f}'.format(a_lofar,a_sig_lofar,b_lofar,b_sig_lofar)

a_em=1.76
b_em=2.00

a_tot_p=1.34
b_tot_p=2.03

a_tot_fe=1.214
b_tot_fe=2.06

print type

for i in np.arange(len(S_use)):
    if type[i]==0:
        energy_recon=tools.rad_to_energy(S_use,a_tot_p,b_tot_p)
    if type[i]==1:
        energy_recon=tools.rad_to_energy(S_use,a_tot_fe,b_tot_fe)

line=np.arange(1e16,1e19,1e16)
#line=np.arange(4e4,4e8,1e8)

x_fit = np.linspace(10**16,10**18.5,100)
#x_data = np.linspace(10**17,10**18,100)



rad_p_fit=sim.energy_to_rad(x_fit,a_tot_p,b_tot_p)
rad_fe_fit=sim.energy_to_rad(x_fit,a_tot_fe,b_tot_fe)

lofar_fit=sim.energy_to_rad(x_fit,a_lofar,b_lofar)
#x=sim_energy*np.sqrt(p_ratio)+sim_energy*np.sqrt(p_ratio)*0.05


y=S_use
x=sim_energy*d_ratio*particle_scale

#x=(radiation_i-radiation_i*0.11)*p_ratio
#y=radiation_p

x_err=x*sigma_e*d_ratio
y_err=y*sigma_e_radio*np.sqrt(p_ratio)



y_label='LOFAR radiation energy (eV)'
x_label='LORA particle energy (eV)'
hist_label='2 (E$_{part}$-E$_{rad}$)/(E$_{part}$+E$_{rad}$)'

'''
fig = plt.figure(figsize=(5,5))
ax1 = fig.add_subplot(1,1,1)#,aspect=1)
ax1.plot(x_fit,rad_p_fit,color='blue',linestyle=':',label='CORSIKA/CoREAS prediction, p')
ax1.plot(x_fit,rad_fe_fit,color='red',linestyle=':',label='CORSIKA/CoREAS prediction, fe')

#ax1.plot(x_fit,lofar_fit,label='LOFAR fit')

ax1.fill_between(x_fit,tools.pred_sys_up(x_fit),tools.pred_sys_down(x_fit),alpha=0.5,color='g')
ax1.plot(x_fit,tools.prediction(x_fit),label="AERA PRL prediction",c='g',linewidth=3)
#x = np.linspace(10**16,10**19.5,100)
#ax1.plot(x, simulation_pred(x), "--",label="CoREAS prediction", c='black',linewidth=3)#

#ax1.plot(x,x,color='black')

#ax1.plot(line,line,color='black')

#ax1.plot(x,y,'.',color='black')
ax1.errorbar(x, y, xerr=x_err,yerr=y_err, marker='o',elinewidth=0.7,linestyle='none',markeredgecolor='black',ecolor='black',label='LOFAR events')

ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel(x_label,fontsize=12)
ax1.set_ylabel(y_label,fontsize=12)

ax1.axis([2e16,3e18,1e4,2e8])

#ax1.axis([3e16,3e18,2e16,3e18])
#ax1.axis([2e16,3e18,7e3,1e8])

plt.tight_layout()


ax1.grid()
plt.show()
raw_input()
plt.close()



######## hist


z=energy_recon

diff=2*(x-z)/(x+z)
bin_start=-2.0
bin_stop=2.0
binWidth=0.1
#nBins=50.0
nBins=int((bin_stop-bin_start)/binWidth)
bins=np.arange(bin_start,bin_stop,binWidth)

p0,p1,p2, x_hist, y_hist=tools.fit_hist(diff,bin_start,bin_stop,int(nBins),-100)



fig = plt.figure(figsize=(4,4))
ax1 = fig.add_subplot(1,1,1)

ax1.hist(diff,bins=bins,alpha=0.5)
ax1.plot(x_hist,y_hist,color='green',label='mean {0:.2f}\n spread {1:.2f}'.format(p1,p2))

ax1.set_xlabel(hist_label,fontsize=12)

ax1.legend(numpoints=1,loc='upper left')
ax1.set_xlim([-1.0, 1.0])

plt.tight_layout()
plt.show()
raw_input()
plt.close()
'''


bin_start=0.0
bin_stop=1.0

binWidthR=0.007
nBinsR=int((bin_stop-bin_start)/binWidthR)
binsR=np.arange(bin_start,bin_stop,binWidthR)

binWidthP=0.02
nBinsP=int((bin_stop-bin_start)/binWidthP)
binsP=np.arange(bin_start,bin_stop,binWidthP)

p0_r,p1_r,p2_r, x_hist_r, y_hist_r=tools.fit_hist(sigma_e_radio,bin_start,bin_stop,int(nBinsR),-0)
p0_p,p1_p,p2_p, x_hist_p, y_hist_p=tools.fit_hist_land(sigma_e,bin_start,bin_stop,int(nBinsP),0)

#bins=np.arange(0,1,0.01)

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1)
ax2 = fig.add_subplot(1,2,2)

ax1.hist(sigma_e_radio,alpha=0.7,bins=binsR)
ax1.plot(x_hist_r, y_hist_r,color='green',linewidth=2,label='mean {0:.2f}\n spread {1:.2f}'.format(p1_r,p2_r))

ax2.hist(sigma_e,alpha=0.7,bins=binsP)
ax2.plot(x_hist_p, y_hist_p,color='green',linewidth=2,label='mean {0:.2f}\n spread {1:.2f}'.format(p1_p,p2_p))

ax1.set_xlabel('std($f_r$)',fontsize=12)
ax2.set_xlabel('std($f_p$)',fontsize=12)


ax1.legend(numpoints=1,loc='upper right')
ax2.legend(numpoints=1,loc='upper right')

ax1.set_xlim([0.0,0.2])
ax2.set_xlim([0.0,0.5])

#ax1.set_yscale('log')
#ax2.set_yscale('log')


plt.show()
raw_input()
plt.close()

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
line=np.arange(0,5,1)
ax1.plot(line,line)
ax1.plot(np.sqrt(p_ratio),d_ratio*(particle_scale),'.')
plt.show()
raw_input()
plt.close()

