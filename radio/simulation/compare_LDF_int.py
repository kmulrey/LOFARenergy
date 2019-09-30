import numpy as np
import os
import cPickle
import matplotlib.pyplot as plt
import sim_functions as sim
import LDF_functions as LDF
from scipy import integrate
from scipy.optimize import minimize


plt.ion()

mag=2.03
average_density = 6.5e-4 #in g/cm^3#atmc.get_density(average_zenith, average_xmax) * 1e-3  # in kg/m^3
p0=0.250524463912
p1=-2.95290494

p0_=0.239
p1_=-3.13

def read_file(filename):
    #try:
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
    
    Erad=Erad+Erad*.11-Erad*0.0336
    S_basic=S_basic+S_basic*.11-S_basic*0.0336
    Srd_1=Srd_1+Srd_1*.11-Srd_1*0.0336
    Srd_2=(Srd_2+Srd_2*.11-Srd_2*0.0336)
    
    corr=((1 - p0 + p0 * np.exp(p1 * (density - average_density)*1e3)) )**2
    corr_=((1 - p0_ + p0_ * np.exp(p1_ * (density - average_density)*1e3)) )**2


    return em_energy,energy,zenith,azimuth,xmax,alpha,S_basic,Srd_1,Srd_2*corr/corr_,Erad,charge_excess_ratio

em_energy_P,energy_P,zenith_P,azimuth_P,xmax_P,alpha_P,S_basic_P,Srd_1_P,Srd_2_P,Erad_P,a_P=read_file('compiled_sims_proton.dat')
em_energy_Fe,energy_Fe,zenith_Fe,azimuth_Fe,xmax_Fe,alpha_Fe,S_basic_Fe,Srd_1_Fe,Srd_2_Fe,Erad_Fe,a_Fe=read_file('compiled_sims_iron.dat')




S_P_use=Srd_1_P
S_Fe_use=Srd_1_Fe

print len(Srd_1_P)

# as written in the LOFAR pipeline
f = lambda x,y,p0,p1,p2,p3,p4,p5: p3*np.exp(-((x-p0)**2 + (y-p1)**2 )/p2**2) - p5*np.exp(- ((x-(p0+p4))**2 + (y-p1)**2) /(np.exp(2.788+0.008*p2))**2)

# no core fitting
f0 = lambda x,y,p2,p3,p4,p5: p3*np.exp(-((x-0)**2 + (y-0)**2 )/p2**2) - p5*np.exp(- ((x-(0+p4))**2 + (y-0)**2) /(np.exp(2.788+0.008*p2))**2)

# as written in the AERA paper
f2 = lambda x,y,p0,p1,p2,p3,p4,p5: p3*np.exp(-((x-p0)**2 + (y-p1)**2 )/p2**2) - p5*np.exp(- ((x-(p0+p4))**2 + (y-p1)**2) /(np.exp(2.788+0.008*p2))**2)

def GetIntegral(A_plus,A_min,sigma,C2=0.0086,C1=2.788,C3=1):  # from LOFAR
    C0 = A_min/A_plus
    E = A_plus*np.pi*((sigma**2)- np.exp(2*sigma*C2+2*C1)*C0)  #integration from Auger paper
    E *= 6.2415e18 # to eV
    return E
'''
def e(pars,X,Y,P):
    ndata=len(P)
    chi2=0
    for l in np.arange(ndata):
        pow=f(X[l],Y[l],pars[0],pars[1],pars[2],pars[3],pars[4],pars[5])
        chi2=chi2+(P[l]-pow)**2/1e-40
    return chi2
'''
def e(pars,X,Y,P):
    ndata=len(P)
    chi2=0
    for l in np.arange(ndata):
        pow=f0(X[l],Y[l],pars[0],pars[1],pars[2],pars[3])
        chi2=chi2+(P[l]-pow)**2/1e-40
    return chi2


def int_by_hand(energy,zenith,azimuth,xmax,x0,xf,y0,yf,nbins):

    divX=(xf-x0)/nbins
    divY=(yf-y0)/nbins
    
    a=divX*divY
    int=0

    for i in np.arange(nbins):
        for j in np.arange(nbins):
            x_pos=x0+i*divX+divX/2
            y_pos=y0+j*divY+divY/2
            int=int+a*LDF.returnPower(energy,zenith,azimuth,xmax,0,0,x_pos,y_pos)

    return int

def int_by_hand_rad(energy,zenith,azimuth,xmax,r0,rf,phi0,phif,nbins):
    
    nbins_phi=36
    div_phi=360/nbins_phi
    
    divR=(rf-r0)/nbins

    int=0
    
    for i in np.arange(nbins):
        a=divR*i*divR*np.pi/180*div_phi
        
        for j in np.arange(nbins_phi):
            phi=j*div_phi
            
            x_pos=divR*i*np.cos(phi)
            y_pos=divR*i*np.sin(phi)

            int=int+a*LDF.returnPower(energy,zenith,azimuth,xmax,0,0,x_pos,y_pos)

    return int




nUse=len(em_energy_Fe)

erad_coreas=[]
erad_int=[]
erad_int_dbl=[]

erad_ldf=[]
res_all=[]
em_energy_keep=[]
energy_keep=[]
zenith_keep=[]
azimuth_keep=[]
alpha_keep=[]
xmax_keep=[]
S_basic=[]


for i in np.arange(nUse):
    try:
        event=i

        energy=energy_Fe[event]#np.power(10.,17.)
        zenith=zenith_Fe[event]#45*np.pi/180
        azimuth=azimuth_Fe[event]#0*np.pi/180
        psi=(2*np.pi-azimuth) ;
        Xcore=0 #shower core
        Ycore=0 #shower core
        Zcore=0 #shower core

        xmax=xmax_P[event]#500#575




    # first guess
        p0=3.53e+01
        p1= -6.5e+00
        p2=1.20e+02
        p3=4.15e-19
        p4=-4.74e+01
        p5= 1.89e-19
    #['X_{+}', 'Y_{+}', '\\sigma_{+}', 'A_{+},', 'X_{-}', 'A_{-}']
    #test=integrate.dblquad(f, 0, np.inf, lambda x: 0, lambda x: np.inf, args=(p0,p1,p2,p3,p4,p5))

        nR=25
        nDeg=8
        x_star=np.zeros([nR*nDeg])
        y_star=np.zeros([nR*nDeg])
        p_star=np.zeros([nR*nDeg])

        count=0

        for r in np.arange(nR):
            for d in np.arange(nDeg):
                rad=(r+1)*25
                deg=(np.pi/180.0)*d*(360.0/nDeg)
                x_star[count]=rad*np.cos(deg)
                y_star[count]=rad*np.sin(deg)
                p_star[count]=LDF.returnPower(energy,zenith,azimuth,xmax,Xcore,Ycore,x_star[count],y_star[count])
                count=count+1






    #['X_{+}', 'Y_{+}', '\\sigma_{+}', 'A_{+},', 'X_{-}', 'A_{-}']
        pars=[p0,p1,p2,p3,p4,p5]
        res=minimize(e,[p2,p3,p4,p5],args=(x_star,y_star,p_star),method='Nelder-Mead', options={'disp': True})


#test_dbl=integrate.dblquad(f0, 0, 7000, lambda x: 0, lambda x: 7000, args=(res['x'][0],res['x'][1],res['x'][2],res['x'][3]))
#test=integrate.dblquad(f0, 0, 700, lambda x: 0, lambda x: 700, args=(res['x'][0],(np.power(10,-52.8)*energy*energy),res['x'][2],res['x'][3]))
        nbins=40
        x0=-600
        xf=600
        y0=-600
        yf=600
        r0=0
        rf=800
        phi0=0
        phif=360
        nbins=50
    #test=int_by_hand(energy,zenith,azimuth,xmax,x0,xf,y0,yf,nbins)
    
        test_r=int_by_hand_rad(energy,zenith,azimuth,xmax,r0,rf,phi0,phif,nbins)
    
    
        alpha=sim.GetAlpha(zenith,azimuth,1.1837)
        int_e_rad= test_r*6.242e18/(np.sin(alpha)**2)
    #int_e_dbl= test_dbl[0]*6.242e18/(np.sin(alpha)**2)

        E_rad = GetIntegral(res['x'][1],res['x'][3],res['x'][0])/(np.sin(alpha)**2)

    # erad_coreas[i]=Srd_1_P[event]
    # erad_int[i]=int_e_rad
    # erad_int_dbl[i]=int_e_dbl

    # erad_ldf[i]=E_rad





        erad_coreas.append(Srd_1_Fe[event])
        erad_int.append(int_e_rad)
        erad_ldf.append(E_rad)
        res_all.append(res['fun'])

        em_energy_keep.append(em_energy_Fe[event])
        energy_keep.append(energy_Fe[event])
        zenith_keep.append(zenith)
        azimuth_keep.append(azimuth)
        alpha_keep.append(alpha)
        xmax_keep.append(xmax)
        S_basic.append(S_basic_Fe[event])


    except:
        print 'issue...'





info={'erad_coreas':np.asarray(erad_coreas),
    'erad_int':np.asarray(erad_int),
    'erad_ldf':np.asarray(erad_ldf),
    'res_all':np.asarray(res_all),
    'em_energy':np.asarray(em_energy_keep),
    'energy':np.asarray(energy_keep),
    'xmax':np.asarray(xmax_keep),
    'alpha':np.asarray(alpha_keep),
    'S_basic':np.asarray(S_basic),
    'azimuth':np.asarray(azimuth_keep),
    'zenith':np.asarray(zenith_keep)}


outfilename='compare_radiation_iron.dat'

outfile=open(outfilename,'w')
cPickle.dump(info,outfile)
outfile.close()









'''
fig=plt.figure(facecolor='white',figsize=(12,4))
ax1=fig.add_subplot(1,3,1)
ax2=fig.add_subplot(1,3,2)
ax3=fig.add_subplot(1,3,3)

line=np.arange(1e3,1e8,1e7)

ax1.plot(line,line,color='black')
ax2.plot(line,line,color='black')
ax3.plot(line,line,color='black')

ax1.plot(erad_coreas,erad_int,'.')
ax2.plot(erad_coreas,erad_ldf,'.')
ax3.plot(erad_int,erad_ldf,'.')

ax1.plot(erad_coreas[good_fit==1],erad_int[good_fit==1],'.',color='red')
ax2.plot(erad_coreas[good_fit==1],erad_ldf[good_fit==1],'.',color='red')
ax3.plot(erad_int[good_fit==1],erad_ldf[good_fit==1],'.',color='red')

ax1.set_yscale('log')
ax1.set_xscale('log')

ax2.set_yscale('log')
ax2.set_xscale('log')

ax3.set_yscale('log')
ax3.set_xscale('log')

ax1.set_xlabel('coreas')
ax1.set_ylabel('int')

ax2.set_xlabel('coreas')
ax2.set_ylabel('ldf')

ax3.set_xlabel('int')
ax3.set_ylabel('ldf')

plt.show()
raw_input()
plt.close()





'''



'''
# make star shaped pattern
print p_star
cmap='gnuplot2_r'
fig=plt.figure(facecolor='white')
ax1=fig.add_subplot(1,1,1, aspect=1)
#ax1.plot(x_star,y_star,'.')
ax1.scatter(x_star, y_star, c=p_star, s=35, cmap=cmap)

plt.show()
raw_input()
plt.close()

'''
