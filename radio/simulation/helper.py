import numpy as np
from optparse import OptionParser
import cPickle
import re
from scipy.signal import hilbert
from scipy.signal import resample
import scipy.fftpack as fftp
import os
from scipy import signal

r_e = 6.371 * 1e6  # radius of Earth
h_max = 11282900.2  # height above sea level where the mass overburden vanishes, cm

def get_atmosphere(h,atm):

    #h and layers in cm
    
    layers=atm[0]
    a=atm[1]
    b=atm[2]
    c=atm[3]
    
    y = np.where(h < layers[0], a[0] + b[0] * np.exp(-1 * h / c[0]), a[1] + b[1] * np.exp(-1 * h / c[1]))
    y = np.where(h < layers[1], y, a[2] + b[2] * np.exp(-1 * h / c[2]))
    y = np.where(h < layers[2], y, a[3] + b[3] * np.exp(-1 * h / c[3]))
    y = np.where(h < layers[3], y, a[4] - b[4] * h / c[4])
    y = np.where(h < h_max, y, 0)
    return y



def get_i_at(at,atm):
    layers=atm[0]
    a=atm[1]
    b=atm[2]
    c=atm[3]
    
    if at > get_atmosphere(layers[0], atm):
        i = 0
    elif at > get_atmosphere(layers[1], atm):
        i = 1
    elif at > get_atmosphere(layers[2], atm):
        i = 2
    elif at > get_atmosphere(layers[3], atm):
        i = 3
    else:
        i = 4
    if i == 4:
        h = -1. * c[i] * (at - a[i]) / b[i]
    else:
        h = -1. * c[i] * np.log((at - a[i]) / b[i])
    return h


    
def get_vertical_height(at,atm):
    return get_i_at(at,atm)




def get_density(h, atm, allow_negative_heights=True):
    """ returns the atmospheric density [g/m^3] for the height h above see level"""
    
    layers=atm[0]
    a=atm[1]
    b=atm[2]
    c=atm[3]
    
    y = 0#np.zeros_like(h, dtype=np.float)
    if not allow_negative_heights:
        y *= np.nan  # set all requested densities for h < 0 to nan
        y = np.where(h < 0, y, b[0] * np.exp(-1 * h / c[0]) / c[0])
    else:
        y = b[0] * np.exp(-1 * h / c[0]) / c[0]
    y = np.where(h < layers[0], y, b[1] * np.exp(-1 * h / c[1]) / c[1])
    y = np.where(h < layers[1], y, b[2] * np.exp(-1 * h / c[2]) / c[2])
    y = np.where(h < layers[2], y, b[3] * np.exp(-1 * h / c[3]) / c[3])
    y = np.where(h < layers[3], y, b[4] / c[4])
    y = np.where(h < h_max, y, 0)
    return y

def get_distance_xmax(zenith, xmax, observation_level=760.): # obs cm
    '''
    input:
        - xmax in g/cm^2
        - zenith in radians
        output: distance to xmax in g/cm^2
    '''
    return get_atmosphere(observation_level) - xmax

def get_distance_for_height_above_ground(h, zenith, observation_level=760.):
    """ inverse of get_height_above_ground() """
    r = r_e + observation_level
    return (h ** 2 + 2 * r * h + r ** 2 * np.cos(zenith) ** 2) ** 0.5 - r * np.cos(zenith)

    
def get_distance_xmax_geometric(zenith, xmax, atm,observation_level=760.):
    '''
    input:
        - xmax in g/cm^2
        - zenith in radians
        - atm params
        output: distance to xmax in cm
    '''

    h = get_vertical_height(xmax,atm) - observation_level
    return get_distance_for_height_above_ground(h, zenith, observation_level)

