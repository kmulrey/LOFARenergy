import numpy as np
import cPickle as pickle
from optparse import OptionParser
import scipy.fftpack as fftp
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy import integrate

C0 = 0.24;
C1 = np.power(10,-52.8);  #J m^-2 eV ^-2
C2 = -7.88; #m
C3 = 28.58; #m
C4 = 1.98; #m
C5 = -2.57; #m
C6 = -54.9; #m
C7 = 0.44; #m g^-1 cm^2
C8 = -1.27e-4; #m *(g^-1 cm^2)^2
C9 = 20.4; #m
C10 = 0.006; # m g^-1 cm^2
C11  = 9.7e-5; # m (g^-1 cm^2)^2
C12_0 = 107; #m
C12_1 = -0.94; # m g^-1 cm^2
C12_2 = 1.94e-3; #m (g^-1 cm^2)^2
C12_3 = -1.5e-6; #m (g^-1 cm^2)^3
C12_4 = 4.1e-10; #m (g^-1 cm^2)^4
Xatm = 1050; #vertically integrated column depth of atmosphere, check # for LOFAR position




def F1(E):
    return C1*E*E

def F2(phi, X, Y, xPrime, yPrime):
    f=(xPrime-(X+C2*np.sin(phi)+C3))*(xPrime-(X+C2*np.sin(phi)+C3))+(yPrime-(Y+C4*np.sin(phi)+C5))*(yPrime-(Y+C4*np.sin(phi)+C5))
    return f

def F3(theta, Xmax):
    f=(C6+C7*(Xatm/np.cos(theta)-Xmax)+C8*(Xatm/np.cos(theta)-Xmax)*(Xatm/np.cos(theta)-Xmax))*(C6+C7*(Xatm/np.cos(theta)-Xmax)+C8*(Xatm/np.cos(theta)-Xmax)*(Xatm/np.cos(theta)-Xmax))
    return f

def F4(theta, Xmax, X, Y, xPrime, yPrime):
    n0=C12_0*np.power((Xatm/np.cos(theta)-Xmax),0)
    n1=C12_1*np.power((Xatm/np.cos(theta)-Xmax),1)
    n2=C12_2*np.power((Xatm/np.cos(theta)-Xmax),2)
    n3=C12_3*np.power((Xatm/np.cos(theta)-Xmax),3)
    n4=C12_4*np.power((Xatm/np.cos(theta)-Xmax),4)
    f=(xPrime-(X+(n0+n1+n2+n3+n4)))*(xPrime-(X+(n0+n1+n2+n3+n4)))+(yPrime-Y)*(yPrime-Y)
    return f

def F5(theta, Xmax):
    f=(C9+C10*(Xatm/np.cos(theta)-Xmax)+C11*(Xatm/np.cos(theta)-Xmax)*(Xatm/np.cos(theta)-Xmax))*(C9+C10*(Xatm/np.cos(theta)-Xmax)+C11*(Xatm/np.cos(theta)-Xmax)*(Xatm/np.cos(theta)-Xmax))
    return f


def returnPower(E, theta, phi, Xmax, X, Y, xPrime, yPrime):
    power=0;
    f1=F1(E)
    f2=F2(phi,X,Y,xPrime,yPrime)
    f3=F3(theta,Xmax)
    f4=F4(theta,Xmax,X,Y,xPrime,yPrime)
    f5=F5(theta,Xmax)
    
    power=f1*np.exp(-1*f2/f3)-C0*f1*np.exp(-1*f4/f5)
    return power



'''
def param_LDF()


def LDF_integrand(x, a, b):
    return a*x**2 + b
'''



















