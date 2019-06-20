import numpy as np
import matplotlib.pyplot as plt
import readLORA as read
import processFiles as process
import constants
import cableLosses as cables

plt.ion()




filename='20141009_101855.root'
#e=1000
for e in np.arange(40,60):
    
    counts=read.readfile(filename,e)

    ldata=len(counts[0])
    dt=2.5
    time=np.arange(0,ldata*dt,dt)
    data=counts[12] # using this one as an example
    if np.sum(data)!=0:
        data_bg_removal=process. avg_baseline_correction(time,data)

        data_clean,sine_fit=process.fix_ringing(time,data)
        #data_clean,sine_fit=process.rolling_correction(time,data)
        data_rolling=process.rolling_correction(np.asarray([time,data_clean]))

        data_losses=cables.apply_cable_loss(time,data_rolling[1])

        counts_bg_removal=process.counts(np.asarray([time,data_bg_removal]))
        counts_cleaned=process.counts(np.asarray([time,data_clean]))
        counts_rolling=process.counts(data_rolling)

        counts_data_losses=process.counts(np.asarray([time,data_losses]))

        # print '{0:.0f},  {1:.0f},  {2:.0f}  --> {3:.2f}'.format(counts_cleaned,counts_rolling,counts_data_losses,counts_data_losses/counts_rolling)

        
        print '{0:.0f},  {1:.0f} --> {2:.2f}'.format(counts_bg_removal,counts_cleaned,counts_cleaned/counts_bg_removal)

        
        
        peak_index=np.argmax(np.abs(data_losses))
        if ((peak_index-constants.WINDOW_OPEN/dt)<=0) or ((peak_index+constants.WINDOW_CLOSE/dt)>(ldata-1)):
            peak_index=int(ldata/2)
    
        start_ind=peak_index-int(constants.WINDOW_OPEN/dt)
        stop_ind=peak_index+int(constants.WINDOW_CLOSE/dt)
        

        fig=plt.figure(facecolor='white')
        ax1=fig.add_subplot(1,1,1)

        # ax1.plot(time,data,label='raw')
        ax1.axvspan(time[start_ind], time[stop_ind], alpha=0.25, color='red')
        ax1.axhline(y=0,color='black')

        #ax1.plot(time,data_bg_removal,label='bg removed')
        ax1.plot(time,data_clean,label='cleaned')
        ax1.plot(time,data_rolling[1],label='rolling')

        ax1.plot(time,data_losses,label='cleaned')

        ax1.legend(loc='upper right', shadow=False,frameon=False,fontsize='x-small')

        ax1.set_xlim([time[start_ind-100],time[stop_ind+200]])

        ax1.set_ylabel('ADC counts')

        ax1.set_xlabel('time (ns)')


        plt.show()
        raw_input()
        plt.close()
