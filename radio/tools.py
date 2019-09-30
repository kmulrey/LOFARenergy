
import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sys
plt.ion()

from ROOT import TH1F
from ROOT import TF1

sys.path.insert(0, 'simulation/')

import fluence as flu
import radiation_energy as rad
import sim_functions as sim
import helper as helper


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
    
    p_ratio=p_ratio
    
    Erad=Erad+Erad*.11-Erad*0.0336
    S_basic=S_basic+S_basic*.11-S_basic*0.0336
    Srd_1=Srd_1+Srd_1*.11-Srd_1*0.0336
    Srd_2=(Srd_2+Srd_2*.11-Srd_2*0.0336)
    
    corr=((1 - p0 + p0 * np.exp(p1 * (density - average_density)*1e3)) )**2
    corr_=((1 - p0_ + p0_ * np.exp(p1_ * (density - average_density)*1e3)) )**2
    
    
    return em_energy,energy,zenith,azimuth,xmax,alpha,S_basic,Srd_1,Srd_2*corr/corr_,Erad,charge_excess_ratio,event,p_ratio,d_ratio,type



def combine_files(coreas_file,ldf_file,uncertainty_file):

    em_energy_int,energy_int,zenith_int,azimuth_int,xmax_int,alpha_int,S_basic_int,Srd_1_int,Srd_2_int,Erad_int,a_int,event_int,p_ratio_i,d_ratio_i,type_int=read_file(coreas_file)
    
    

    parameter_file=open(ldf_file,'r')
    info=cPickle.load(parameter_file)
    parameter_file.close
    event_number=info['event_number']
    lora_energy_p=info['lora_energy']
    radiation_energy_p=np.asarray(info['radiation_energy'])
    alpha_param=np.asarray(info['alpha'])
    lofar_e_err=np.asarray(info['lofar_e_err'])
    lora_e_err=np.asarray(info['lora_e_err'])
    
    
    

    uncertainty_file=open(uncertainty_file,'r')
    info=cPickle.load(uncertainty_file)

    uncertainty_file.close

    event_u=info['event']
    sigma_e_u=info['sigma_e']
    sigma_e_radio_u=info['sigma_e_radip']
    sigma_logE_u=info['sigma_logE']
    sigma_logE_radio_u=info['sigma_logE_radio']
    sigma_r_u=info['sigma_r']

    
    lora_energy=np.empty([0])
    radiation_i=np.empty([0])
    radiation_p=np.empty([0])
    alpha_i=np.empty([0])
    alpha_p=np.empty([0])

    azimuth=np.empty([0])
    zenith=np.empty([0])
    event=np.empty([0])
    em_energy=np.empty([0])

    sim_energy=np.empty([0])
    xmax=np.empty([0])
    p_ratio=np.empty([0])
    d_ratio=np.empty([0])
    type=np.empty([0])
    lora_err=np.empty([0])
    lofar_err=np.empty([0])
    combchi2=np.empty([0])
    chi2=np.empty([0])
    S_basic=np.empty([0])
    Srd_1=np.empty([0])
    Srd_2=np.empty([0])
    
    ldf_radio_err=np.empty([0])
    ldf_particle_err=np.empty([0])

    sigma_e=np.empty([0])
    sigma_e_radio=np.empty([0])
    sigma_logE=np.empty([0])
    sigma_logE_radio=np.empty([0])
    sigma_r=np.empty([0])
    
    
    
    nEvents=len(event_int)

    for e in np.arange(nEvents):
        n=event_int[e]
        index_i=np.where(event_int == n)[0]
        index_p=np.where(event_number == n)[0]
        index_u=np.where(event_u == n)[0]
    
    
    
        if len(index_p)>0 and len(index_i>0) and len(index_u>0) and sigma_e_radio_u[index_u]<1.0 and p_ratio_i[index_i]<5 and p_ratio_i[index_i]>0.1 and sigma_e_u[index_u]<1.3:
        
            lora_energy=np.concatenate((lora_energy,np.asarray(lora_energy_p[index_p])),axis=0)
            ldf_particle_err=np.concatenate((ldf_particle_err,np.asarray(lora_e_err[index_p])),axis=0)
            ldf_radio_err=np.concatenate((ldf_radio_err,np.asarray(lofar_e_err[index_p])),axis=0)
        
            scale=np.sin(alpha_param[index_p])**2/np.sin(alpha_int[index_i])**2
        
            radiation_i=np.concatenate((radiation_i,(1/(np.sin(alpha_int[index_i])**2))*np.asarray(Erad_int[index_i])),axis=0)
            radiation_p=np.concatenate((radiation_p,scale*np.asarray(radiation_energy_p[index_p])),axis=0)
        
            azimuth=np.concatenate((azimuth,np.asarray(azimuth_int[index_i])),axis=0)
            zenith=np.concatenate((zenith,np.asarray(zenith_int[index_i])),axis=0)
        
            event=np.concatenate((event,np.asarray(event_int[index_i])),axis=0)
        
            em_energy=np.concatenate((em_energy,np.asarray(em_energy_int[index_i])),axis=0)
        
            xmax=np.concatenate((xmax,np.asarray(xmax_int[index_i])),axis=0)
        
            p_ratio=np.concatenate((p_ratio,np.asarray(p_ratio_i[index_i])),axis=0)
            d_ratio=np.concatenate((d_ratio,np.asarray(d_ratio_i[index_i])),axis=0)
        
            type=np.concatenate((type,np.asarray(type_int[index_i])),axis=0)
            #combchi2=np.concatenate((combchi2,np.asarray(combchi2_i[index_i])),axis=0)
            #chi2=np.concatenate((chi2,np.asarray(radiochi2_i[index_i])),axis=0)

            S_basic=np.concatenate((S_basic,np.asarray(S_basic_int[index_i])),axis=0)
            Srd_1=np.concatenate((Srd_1,np.asarray(Srd_1_int[index_i])),axis=0)
            Srd_2=np.concatenate((Srd_2,np.asarray(Srd_2_int[index_i])),axis=0)

        
            alpha_i=np.concatenate((alpha_i,np.asarray(alpha_int[index_i])),axis=0)
            alpha_p=np.concatenate((alpha_p,np.asarray(alpha_param[index_p])),axis=0)
            sim_energy=np.concatenate((sim_energy,np.asarray(energy_int[index_i])),axis=0)

        
            sigma_e=np.concatenate((sigma_e,np.asarray(sigma_e_u[index_u])),axis=0)
            sigma_e_radio=np.concatenate((sigma_e_radio,np.asarray(sigma_e_radio_u[index_u])),axis=0)
            sigma_logE=np.concatenate((sigma_logE,np.asarray(sigma_logE_u[index_u])),axis=0)
            sigma_logE_radio=np.concatenate((sigma_logE_radio,np.asarray(sigma_logE_radio_u[index_u])),axis=0)
            sigma_r=np.concatenate((sigma_r ,np.asarray(sigma_r_u[index_u])),axis=0)


    return event,radiation_i,radiation_p,azimuth,zenith,xmax,em_energy,sim_energy,p_ratio,d_ratio,type,alpha_i,alpha_p,sigma_e,sigma_e_radio,lora_energy,ldf_particle_err,ldf_radio_err,S_basic,Srd_1,Srd_2

def rad_to_energy(x,a,b):
    return np.power((x/(a*1e7)),1/b)*1e18

def energy_to_rad(x,a,b):
    return a*1e7*(x/1e18)**b#*(mag)**b


def fit_hist(counts,bin_start,bin_stop,nBins,min):
    counts=counts[counts>min]
    binWidth=(bin_stop-bin_start)/nBins
    bins=np.arange(bin_start,bin_stop,binWidth)
    
    print bin_stop
    print bin_start
    Hist = TH1F("Hist", "Hist_x;x-axis;Frequency",int(nBins),bin_start,bin_stop)
    #fit = TF1("fit", "gaus", bin_start-1,bin_stop+1)
    fit = TF1("fit", "gaus", -0.2,0.2)
    
    for i in np.arange(len(counts)):
        Hist.Fill(counts[i])
    
    Hist.Fit("fit",'0')

    p0=fit.GetParameter(0)
    p1=fit.GetParameter(1)
    p2=fit.GetParameter(2)

    x=np.arange(bin_start,bin_stop,binWidth/5)
    y=np.zeros([len(x)])
    for i in np.arange(len(x)):
        y[i]=fit.Eval(x[i])


    '''
    fig = plt.figure()
    ax1 = fig.add_subplot(1,1,1)



    ax1.hist(counts,alpha=0.7,bins=bins)
    ax1.plot(x, y,color='green',linewidth=2,label='mean {0:.2f}\n spread {1:.2f}'.format(p1,p2))

    plt.show()
    raw_input()
    plt.close()
    '''
    
    return p0, p1, p2, x, y

def fit_hist_land(counts,bin_start,bin_stop,nBins,min):
    
    counts=counts[counts>min]
    binWidth=(bin_stop-bin_start)/nBins
    bins=np.arange(bin_start,bin_stop,binWidth)
    Hist = TH1F("Hist", "Hist_x;x-axis;Frequency",int(nBins),int(bin_start),int(bin_stop))
    fit = TF1("fit", "landau", bin_start,bin_stop)
    #fit = TF1("fit", "gaus", -200,10000)
    
    for i in np.arange(len(counts)):
        Hist.Fill(counts[i])

    Hist.Fit("fit",'0')

    p0=fit.GetParameter(0)
    p1=fit.GetParameter(1)
    p2=fit.GetParameter(2)

    x=np.arange(bin_start,bin_stop,0.005)

    y=np.zeros([len(x)])
    for i in np.arange(len(x)):
        y[i]=fit.Eval(x[i])

    
    return p0, p1, p2, x, y

def prediction(E,index=1.8):
    # PRL prediction from AERA
    E_rad = 15.8e6*(E*1e-18)**2#*2.04**index
    return E_rad
def pred_sys_up(E,index=1.8):
    # Plus uncertainties
    E_rad = (15.8+6.8)*10**6*(E*1e-18)**2#*(2.04**index)
    return E_rad
def pred_sys_down(E,index=1.8):
    # Minus uncertainties
    E_rad = (15.8-6.8)*10**6*(E*1e-18)**2#*(2.04**index)
    return E_rad
def inverse_pred(E_rad):
    E = 1.348*1e14*np.sqrt(E_rad)
    return E

def simulation_pred(E,index=1.8):
    E_rad = 11.9e6*(E*1e-18)**2*2.04**index
    return E_rad


def chi2(y1,y2,sig1):
    chi2=0
    for i in np.arange(len(y1)):
        chi2=chi2+((y1[i]-y2[i])**2)/(sig1[i]**2)
    return chi2


def e(pars,energy,rad,std):
    
    
    #pars a,b
    a=pars[0]
    b=pars[1]
    
    y=a*1e7*np.power((energy*1e-18),b)

    
    #logR=np.log10(a)+7+b*np.log10(energy*1e-18)#+b*np.log10(mag)
    #log_r=np.log10(rad)
    
  
    #log_std=((b*std_use*1e-18)/((energy_use*1e-18)*np.log(10)))**2
    #X2=np.sum(((logR-log_r)**2)/std**2)/len(energy)
    X2=np.sum(((np.log10(rad)-np.log10(y))**2)/(np.log10(std))**2)
    #X2=np.sum(((np.log10(rad)-np.log10(y))**2)/(np.log10(energy)))/len(energy)

    return X2
