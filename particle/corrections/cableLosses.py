import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

plt.ion()


def apply_cable_loss(time,data):
    
    cable_length=60.0  # meters->guess

    atten_freq=np.asarray([0.0,1.8,10.0,28.0,50.0,100.0,144.0,200.0,400.0])*1e6
    atten_dB=np.asarray([0.0,2.1,4.7,7.6,10.8,15.8,19.3,22.1,33.3])*cable_length/100.0
    
    
    f2 = interp1d(atten_freq, atten_dB)
    xnew = np.linspace(0, 400, 100, endpoint=True)
    
    spec  = np.fft.rfft(data)
    freq = np.fft.rfftfreq(len(time), (time[1] - time[0])*1e-9)
    atten_interp=f2(freq)
    ratio=1./(np.power(10.0,(atten_interp/20.0)))
    spec1=spec*ratio
    
    
    data_loss=np.fft.irfft(spec1)
    
    '''
    fig=plt.figure(facecolor='white')
    ax1=fig.add_subplot(1,2,1)
    ax2=fig.add_subplot(1,2,2)
    
    ax1.plot(time,data)
    ax1.plot(time,data_loss)


    ax2.plot(freq*1e-6,np.abs(spec)**2*1e-6)
    ax2.plot(freq*1e-6,np.abs(spec1)**2*1e-6)

    
    ax2.set_yscale('log', nonposy='clip')


    ax2.set_xlabel('MHz')
    ax2.set_ylabel('|fft|^2/MHz')
    
    


    plt.show()
    raw_input()
    plt.close()
    '''

    return data_loss







