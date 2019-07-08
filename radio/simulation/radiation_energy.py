import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import helper as helper

#from radiotools.atmosphere import models as atm

# see Glaser et al., JCAP 09(2016)024 for the derivation of the formulas

#average_xmax = 669.40191244545326  # 1 EeV, 50% proton, 50% iron composition
#average_zenith = np.deg2rad(45)
#atmc = atm.Atmosphere(model=1)
average_density = 6.5e-4 #in g/cm^3#atmc.get_density(average_zenith, average_xmax) * 1e-3  # in kg/m^3


def get_clipping(dxmax):
    """ get clipping correction
        
        Parameters
        ----------
        dxmax : float
        distance to shower maximum in g/cm^2
        
        Returns
        -------
        float
        fraction of radiation energy that is radiated in the atmosphere
        
        """
    return 1 - np.exp(-8.7 * (dxmax * 1e-3 + 0.29) ** 1.89)


def get_a(rho):
    """ get relative charge excess strength
        
        Parameters
        ----------
        rho : float
        density at shower maximum in g/cm^3
        
        Returns
        -------
        float
        relative charge excess strength a
        """
    return -0.23604683 + 0.43426141 * np.exp(1.11141046 * ((rho - average_density)*1e3))

'''
def get_a_zenith(zenith):
    """ get relative charge excess strength wo Xmax information
        
        Parameters
        ----------
        zentith : float
        zenith angle in rad according to radiotools default coordinate system
        
        Returns
        --------
        float
        relative charge excess strength a
        """
    rho = atmc.get_density(zenith, average_xmax) * 1e-3
    return -0.24304254 + 0.4511355 * np.exp(1.1380946 * (rho - average_density))
'''

def get_S_basic(Erad, sinalpha,b_scale=2.03, b=1.8):
    """ get corrected radiation energy (S_RD) no charge excess or density
        
        Parameters
        ----------
        Erad : float
        radiation energy (in eV)
        sinalpha: float
        sine of angle between shower axis and geomagnetic field
       
        
        Returns
        --------
        float:
        corrected radiation energy (in eV)
        """
    return Erad / (sinalpha ** 2 * b_scale ** b)



def get_S(Erad, sinalpha, density, p0=0.250524463912, p1=-2.95290494,
          b_scale=2.03, b=1.8):
    """ get corrected radiation energy (S_RD)
        
        Parameters
        ----------
        Erad : float
        radiation energy (in eV)
        sinalpha: float
        sine of angle between shower axis and geomagnetic field
        density : float
        density at shower maximum in g/cm^3
        
        Returns
        --------
        float:
        corrected radiation energy (in eV)
        """
    a = get_a(density) /( b_scale ** (0.9))
    
    Erad_corr=(Erad/(a**2+(1-a**2)*sinalpha**2*b_scale**b))#*((1)/(1-p0+p0*np.exp(p1*(density-average_density)*1000)**2))
    
    return Erad_corr #Erad / (a ** 2 + (1 - a ** 2) * sinalpha ** 2 * b_scale ** b) / \(1 - p0 + p0 * np.exp(p1 * (density - average_density)*1e3)) ** 2



def get_radiation_energy(Srd, sinalpha, density, p0=0.239,
                         p1=-3.13, b_scale=2.03, b=1.8):
    """ get radiation energy (S_RD)
        
        Parameters
        ----------
        Srd : float
        corrected radiation energy (in eV)
        sinalpha: float
        sine of angle between shower axis and geomagnetic field
        density : float
        density at shower maximum in g/cm^3
        
        Returns
        --------
        float:
        radiation energy (in eV)
        """
    return Srd/((1 - p0 + p0 * np.exp(p1 * (density - average_density)*1e3)) )**2





def integral(fluence,pos):

    pos_uvw_vxb=pos[0::8]
    pos_uvw_vxvxb=pos[2::8]
    neg_uvw_vxb=pos[4::8]
    neg_uvw_vxvxb=pos[6::8]
    
    fluence_pos_vxb=fluence[0::8]
    fluence_pos_vxvxb=fluence[2::8]
    fluence_neg_vxb=fluence[4::8]
    fluence_neg_vxvxb=fluence[6::8]
        
        
    
    pos_vxvxb_all=np.concatenate([neg_uvw_vxvxb.T[1],pos_uvw_vxvxb.T[1]])
    fluence_vxvxb=np.concatenate([fluence_pos_vxvxb,fluence_neg_vxvxb])
    inds = pos_vxvxb_all.argsort()
        
    sorted_pos=pos_vxvxb_all[inds]
    sorted_fluence_vxvxb=fluence_vxvxb[inds]
    f0 = interp1d(sorted_pos, sorted_fluence_vxvxb, kind='cubic')
    r= np.linspace(0, 500, num=1000, endpoint=True)
    
    
    # integrate positive vxvxB arm
    n=len(r)
    dr=r[1]-r[0]
    integral=0
    
    for i in np.arange(n-1):
        r0=r[i]
        r1=r[i+1]
        val0=r0*(f0(r0))
        val1=r1*(f0(r1))
        integral=integral+(val0+val1)*0.5*dr

    integral=integral*2*np.pi

    #print '{0:.2f}'.format(integral)
    '''
    plt.ion()
    fig=plt.figure(facecolor='white')
    ax1=fig.add_subplot(1,1,1)
    #ax1.plot(pos_uvw_vxb.T[0],pos_uvw_vxb.T[1],'.',label='+ vxb')
    #ax1.plot(neg_uvw_vxb.T[0],neg_uvw_vxb.T[1],'.',label='- vxb')
    #ax1.plot(pos_uvw_vxvxb.T[0],pos_uvw_vxvxb.T[1],'.',label='+ vxvxb')
    #ax1.plot(neg_uvw_vxvxb.T[0],neg_uvw_vxvxb.T[1],'.',label='- vxvxb')

    
    #ax1.plot(ant_pos_shower.T[0],ant_pos_shower.T[1],'.')
    #ax1.plot(pos_uvw_vxb.T[0],fluence_pos_vxb,'.',label='pos vxB')
    #ax1.plot(pos_uvw_vxvxb.T[1],fluence_pos_vxvxb,'.',label='pos vxvxB')
    #ax1.plot(neg_uvw_vxb.T[0],fluence_neg_vxb,'.',label='neg vxB')
    #x1.plot(neg_uvw_vxvxb.T[1],fluence_neg_vxvxb,'.',label='neg vxvxB')
    ax1.plot(sorted_pos,sorted_fluence_vxvxb,'.',label='vxvxB')
    ax1.plot(xnew,f0(xnew),'.',label='vxvxB')

    #ax1.plot(efield[a][1])
    #ax1.plot(efield[a][2])
    ax1.legend(loc='upper left', shadow=False,frameon=False)
    ax1.set_xlabel('vxvxB position (m)')
    ax1.set_ylabel('fluence')

    plt.show()
    raw_input()
    plt.close()
    '''

    return integral


def StoEm(S,A=1.683,B=2.006):
    
    Em=np.power((S/(A*1e7)),1/B)*1e18

    return Em

