import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sim_functions as sim


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
    
    Erad=Erad+Erad*.11-Erad*0.0336
    S_basic=S_basic+S_basic*.11-S_basic*0.0336
    Srd_1=Srd_1+Srd_1*.11-Srd_1*0.0336
    Srd_2=(Srd_2+Srd_2*.11-Srd_2*0.0336)
    
    corr=((1 - p0 + p0 * np.exp(p1 * (density - average_density)*1e3)) )**2
    corr_=((1 - p0_ + p0_ * np.exp(p1_ * (density - average_density)*1e3)) )**2

    print charge_excess_ratio

    return em_energy,energy,zenith,azimuth,xmax,alpha,S_basic,Srd_1,Srd_2*corr/corr_,Erad,charge_excess_ratio

em_energy_P,energy_P,zenith_P,azimuth_P,xmax_P,alpha_P,S_basic_P,Srd_1_P,Srd_2_P,Erad_P,a_P=read_file('compiled_sims_proton.dat')
em_energy_Fe,energy_Fe,zenith_Fe,azimuth_Fe,xmax_Fe,alpha_Fe,S_basic_Fe,Srd_1_Fe,Srd_2_Fe,Erad_Fe,a_Fe=read_file('compiled_sims_iron.dat')

print len(em_energy_P)
print len(em_energy_Fe)



#a_em_p,b_em_p,a_sig_em_p,b_sig_em_p=sim.getFit(em_energy_P,Erad_P/(np.sin(alpha_P)**2))
#a_em_fe,b_em_fe,a_sig_em_fe,b_sig_em_fe=sim.getFit(em_energy_Fe,Erad_Fe/(np.sin(alpha_Fe)**2))

a_P_corr=a_P/(np.power(mag,0.9))
a_Fe_corr=a_Fe/(np.power(mag,0.9))

#E_use_P=Erad_P/(a_P**2+(1-a_P**2)*(np.sin(alpha_P)**2))
#E_use_Fe=Erad_Fe/(a_Fe**2+(1-a_Fe**2)*(np.sin(alpha_Fe)**2))

E_use_P=Erad_P/(a_P_corr**2+(1-a_P_corr**2)*(np.sin(alpha_P)**2)*np.power(mag,1.8))
E_use_Fe=Erad_Fe/(a_Fe_corr**2+(1-a_Fe_corr**2)*(np.sin(alpha_Fe)**2)*np.power(mag,1.8))


#a_em_p,b_em_p,a_sig_em_p,b_sig_em_p=sim.getFit(em_energy_P,E_use_P)
#a_em_fe,b_em_fe,a_sig_em_fe,b_sig_em_fe=sim.getFit(em_energy_Fe,E_use_Fe)

S_P_use=Srd_1_P
S_Fe_use=Srd_1_Fe

a_em_p,b_em_p,a_sig_em_p,b_sig_em_p=sim.getFit(em_energy_P,S_P_use)
a_em_fe,b_em_fe,a_sig_em_fe,b_sig_em_fe=sim.getFit(em_energy_Fe,S_Fe_use)

a_tot_p,b_tot_p,a_sig_tot_p,b_sig_tot_p=sim.getFit(energy_P,S_P_use)
a_tot_fe,b_tot_fe,a_sig_tot_fe,b_sig_tot_fe=sim.getFit(energy_Fe,S_Fe_use)


print '\n\nEM energy fit'
print '___________________'
print '{0:.3f} +/- {1:.3f}       {2:.3f} +/- {3:.3f}'.format(a_em_p,a_sig_em_p,b_em_p,b_sig_em_p)
print '{0:.3f} +/- {1:.3f}       {2:.3f} +/- {3:.3f}'.format(a_em_fe,a_sig_em_fe,b_em_fe,b_sig_em_fe)


print '\n\nTotal energy fit'
print '___________________'
print '{0:.3f} +/- {1:.3f}       {2:.3f} +/- {3:.3f}'.format(a_tot_p,a_sig_tot_p,b_tot_p,b_sig_tot_p)
print '{0:.3f} +/- {1:.3f}       {2:.3f} +/- {3:.3f}'.format(a_tot_fe,a_sig_tot_fe,b_tot_fe,b_sig_tot_fe)


x = np.linspace(10**16,10**18.5,100)
x_data = np.linspace(10**17,10**18,100)

rad_em_fit_p=sim.energy_to_rad(x,a_em_p,b_em_p)
rad_em_fit_fe=sim.energy_to_rad(x,a_em_fe,b_em_fe)

rad_tot_fit_p=sim.energy_to_rad(x,a_tot_p,b_tot_p)
rad_tot_fit_fe=sim.energy_to_rad(x,a_tot_fe,b_tot_fe)
#rad_em_fit_p=sim.energy_to_rad(x,1.64,2)
#rad_em_fit_fe=sim.energy_to_rad(x,1.64,2)


#x = np.linspace(10**4,10**8,100)

fig = plt.figure()
ax1 = fig.add_subplot(1,2,1,aspect=1)
ax2 = fig.add_subplot(1,2,2,aspect=1)


#ax1.plot(x,x,color='black')

ax1.plot(em_energy_P,S_P_use,'.',color='blue',alpha=0.2,label='proton')
ax1.plot(x,rad_em_fit_p,color='blue')

ax1.plot(em_energy_Fe,S_Fe_use,'.',color='red',alpha=0.2,label='iron')
ax1.plot(x,rad_em_fit_fe,color='red')

ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax1.set_xlabel('EM Energy (eV)',fontsize=12)
ax1.set_ylabel('S$_{RD} \ (eV)$',fontsize=12)

ax1.axis([2e16,3e18,7e3,1e8])


###
ax2.plot(energy_P,S_P_use,'.',color='blue',alpha=0.2,label='proton')
ax2.plot(x,rad_tot_fit_p,color='blue')


ax2.plot(energy_Fe,S_Fe_use,'.',color='red',alpha=0.2,label='iron')
ax2.plot(x,rad_tot_fit_fe,color='red')

ax2.set_yscale('log')
ax2.set_xscale('log')

ax2.legend(numpoints=1,loc='upper left',framealpha=1,facecolor='white')

ax2.set_xlabel('Total CR Energy (eV)',fontsize=12)
ax2.set_ylabel('S$_{RD} (eV)$',fontsize=12)
ax2.axis([2e16,3e18,7e3,1e8])

ax1.grid()
ax2.grid()
plt.show()
raw_input()
plt.close()
