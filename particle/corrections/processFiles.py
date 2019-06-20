import numpy as np
import math
import matplotlib.pyplot as plt
from ROOT import TH1F
from ROOT import TF1
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.interpolate import UnivariateSpline
from scipy.signal import find_peaks
import constants
plt.ion()


def process_list(list):
    nLines=len(list)
    
    
    ch1_start=-1
    ch2_start=-1
    ch3_start=-1
    ch4_start=-1

    for i in np.arange(nLines):
        if 'Channel 1' in list[i]:
            ch1_start=i
    
        if 'Channel 2' in list[i]:
            ch2_start=i

        if 'Channel 3' in list[i]:
            ch3_start=i
            
        if 'Channel 4' in list[i]:
            ch4_start=i
                
    ch1=[]
    ch2=[]
    ch3=[]
    ch4=[]
  
    for i in np.arange(ch1_start+1,ch2_start):
        hold=list[i].strip().split(' ')
        str_list = filter(None, hold) # fastest
        num_array = [ float(s) for s in str_list]
        ch1.append(num_array)
    flattened_ch1 = np.asarray([val for sublist in ch1 for val in sublist])

    for i in np.arange(ch2_start+1,ch3_start):
        hold=list[i].strip().split(' ')
        str_list = filter(None, hold) # fastest
        num_array = [ float(s) for s in str_list]
        ch2.append(num_array)
    flattened_ch2 = np.asarray([val for sublist in ch2 for val in sublist])

    for i in np.arange(ch3_start+1,ch4_start):
        hold=list[i].strip().split(' ')
        str_list = filter(None, hold) # fastest
        num_array = [ float(s) for s in str_list]
        ch3.append(num_array)
    flattened_ch3 = np.asarray([val for sublist in ch3 for val in sublist])

    for i in np.arange(ch4_start+1,len(list)):
        hold=list[i].strip().split(' ')
        str_list = filter(None, hold) # fastest
        num_array = [ float(s) for s in str_list]
        ch4.append(num_array)
    flattened_ch4 = np.asarray([val for sublist in ch4 for val in sublist])


    return flattened_ch1,flattened_ch2,flattened_ch3,flattened_ch4

def process_file(filename):
    
  
    data_flag=0

    event_line_count=0
    event_count=0
    file=open(filename)
    
    # count event size
    with file as fp:
        list=[]
        c=0
        new_c=0
        event_lines=[]
        
        diff=[]
        
        for cnt, line in enumerate(fp):
            if 'Event record' in line:
                if cnt>1:
                    diff.append((cnt-new_c))
                new_c=cnt
    file.close()

    if len(set(diff)) <= 1 == False:
        print 'diff'
        data_flag=1
    #print diff
    block_size=diff[0]

    if data_flag==0:
        file=open(filename)
        with file as fp:
            list=[]
            c=0
            new_c=0
            event_lines=[]
            listch1=[]
            listch2=[]
            listch3=[]
            listch4=[]

            for cnt, line in enumerate(fp):
                event_line_count=event_line_count+1
                #if 'Event record' in line:
                #print (cnt-new_c),line
                #new_c=cnt
                if cnt%block_size==0:
                #print 'new event', cnt
                    event_count=event_count+1

                    event_lines=[]
                if cnt%block_size<(block_size-1):
                    event_lines.append(line)
                elif cnt%block_size==(block_size-1):
                #print 'end',len(event_lines), cnt
                #if event_count==1:
                    test1,test2,test3,test4=process_list(event_lines)
                    listch1.append(test1)
                    listch2.append(test2)
                    listch3.append(test3)
                    listch4.append(test4)

    file.close()
    try:
        data_ch1= np.asarray(listch1)
        data_ch2= np.asarray(listch2)
        data_ch3= np.asarray(listch3)
        data_ch4= np.asarray(listch4)
    except:
        data_ch1=np.empty([])
        data_ch2=np.empty([])
        data_ch3=np.empty([])
        data_ch4=np.empty([])
        print 'issue with file'

    return data_ch1,data_ch2,data_ch3,data_ch4


def findCounts(data):
    

    nTraces=len(data)
    l=len(data[0])
    bg=np.zeros([nTraces])
    data_adjusted=np.zeros([nTraces,l])
    peaks=[]#np.zeros([nTraces])
    counts=[]#np.zeros([nTraces])
    
    # findbackground to subtract
    #try:
    for t in np.arange(1):
        for i in np.arange(nTraces):
        
            rms=nanrms(data[i][int(bg_min/dt):int(bg_max/dt)])
            
    
            bg[i]=np.average(data[i][int(bg_min/dt):int(bg_max/dt)])
            data_adjusted[i]=data[i]-bg[i]
            peak_index=np.argmax(np.abs(data_adjusted[i]))
            if (peak_index<(window_open/dt+10)) or (peak_index>(l-window_close/dt-10)):
                peak_index=int(l/2)
            start_ind=peak_index-int(window_open/dt)
            stop_ind=peak_index+int(window_close/dt)
            #data_use=data_adjusted[i]
            #data_use=data_use[peaks[i]>20]
            #counts[i]=np.abs(np.sum(data_adjusted[i][start_ind:stop_ind]))
            #peaks[i]=np.max(np.abs(data_adjusted[i][start_ind:stop_ind]))
            
            peak=np.max(np.abs(data_adjusted[i][start_ind:stop_ind]))
            #if peak>2*rms:
            counts.append(np.abs(np.sum(data_adjusted[i][start_ind:stop_ind])))
            peaks.append(np.max(np.abs(data_adjusted[i][start_ind:stop_ind])))

                #except:
                # counts[i]=0

    return np.asarray(peaks),np.asarray(counts)




def fit_peak(counts,bin_start,bin_stop,nBins,min):
    
    counts=counts[counts>min]
    binWidth=(bin_stop-bin_start)/nBins
    bins=np.arange(bin_start,bin_stop,binWidth)
    Hist = TH1F("Hist", "Hist_x;x-axis;Frequency",int(nBins),int(bin_start),int(bin_stop))
    fit = TF1("fit", "landau", -2,10000)
    #fit = TF1("fit", "gaus", -200,10000)

    for i in np.arange(len(counts)):
        Hist.Fill(counts[i])
    
    Hist.Fit("fit",'0')

    p0=fit.GetParameter(0)
    p1=fit.GetParameter(1)
    p2=fit.GetParameter(2)
    
    x=np.arange(-100,10000,1)
    y=np.zeros([len(x)])
    for i in np.arange(len(x)):
        y[i]=fit.Eval(x[i])

    return p0, p1, p2, x, y





def nanrms(x, axis=None):
    return np.sqrt(np.nanmean(x**2, axis=axis))



def my_sin(x, freq, amplitude, phase, offset):
    return (np.sin(2*np.pi*(freq*x) + phase) * amplitude + offset)




# to remove the 100 MHZ ringing from existing LORA data

def fix_ringing(time,data):

    ldata=len(time)
    dt=time[1]-time[0]
    guess_freq = (1/ldata)
    guess_amplitude = 5
    guess_phase = np.pi
    guess_offset = -20
    guess_freq = 1/(time[4]-time[0])
    guess_amplitude = 5
    guess_phase = np.pi
    guess_offset = 0
    
    peak_index=np.argmax(np.abs(data))
    
    # check that we aren't too close to beginning or end
    
    if (peak_index-int(constants.WINDOW_OPEN/dt)<=0) or (peak_index+int(constants.WINDOW_CLOSE/dt))>=(ldata-1):
        peak_index=int(ldata/2)
    
    
    x=np.concatenate((time[0:(peak_index-int(constants.WINDOW_OPEN/dt))],time[(peak_index+int(constants.WINDOW_CLOSE/dt)):]))
    y=np.concatenate((data[0:(peak_index-int(constants.WINDOW_OPEN/dt))],data[(peak_index+int(constants.WINDOW_CLOSE/dt)):]))
    
    
    p0=[guess_freq, guess_amplitude,guess_phase, guess_offset]
    
    try:
        fit = curve_fit(my_sin, x, y, p0=p0)
    #time2=np.arange(fitting_time[0],fitting_time[0]+500,0.1)
    #data_first_guess = my_sin(x, *p0)
    #data_fit = my_sin(time2, *fit[0])
    #counts_fit_use= y-data_first_guess

    #fit = curve_fit(my_sin, x, counts_fit_use, p0=fit[0])
        counts_fit= data-my_sin(time, *fit[0])
        return counts_fit,my_sin(time, *fit[0])

    except:
        counts_fit=data
        return counts_fit,my_sin(time, guess_freq, guess_amplitude,guess_phase, guess_offset)


def rolling_average(x,y):


    peak_index=np.argmax(np.abs(y))
    
    if (peak_index<(window_open/dt+10)) or (peak_index>(len(x)-window_close/dt-10)):
        peak_index=int(len(x)/2)
    
    
    #print peak_index
    start_ind=peak_index-int(window_open/dt)
    stop_ind=peak_index+int(window_close/dt)
    start_avg=np.average(y[:n_avg])
    stop_avg=np.average(y[(len(y)-n_avg-1):])

    x_fit=np.concatenate((x[0:(peak_index-int(window_open/dt))], x[(peak_index+int(window_close/dt)):]))
    y_fit=np.concatenate((y[0:(peak_index-int(window_open/dt))], y[(peak_index+int(window_close/dt)):]))
    
    roll=np.zeros(len(y))
    total_window=len(x)
    pre=len(x[0:(peak_index-int(window_open/dt))])
    sig=len(x[(peak_index-int(window_open/dt)):(peak_index+int(window_close/dt))])
    post=len(x[(peak_index+int(window_close/dt)):len(x)])

    for i in np.arange(len(y)):
        if i<n_avg:
            roll[i]=start_avg
        elif i>(len(y)-n_avg-1):
            roll[i]=stop_avg
        else:
            roll[i]=np.average(y[(i-int(n_avg/2)):(i+int(n_avg/2))])

    first_sig=np.average(y[(peak_index-int(window_open/dt)-n_avg-5):(peak_index-int(window_open/dt)-5)])
    last_sig=np.average(y[(peak_index+int(window_close/dt)+5):(peak_index+int(window_close/dt)+n_avg)+5])
    #print first_sig,last_sig
    

    for i in np.arange((peak_index-int(window_open/dt)-10),(peak_index+int(window_close/dt)+10)):
        roll[i]=np.average([first_sig,first_sig])




    return x,roll


def rolling_average_var(x,y,dt,window_open,window_close,n_avg):
    window_open=70
    window_close=1000
    
    peak_index=np.argmax(np.abs(y))
    
    if (peak_index<(window_open/dt+20)) or (peak_index>(len(x)-window_close/dt-20)):
        peak_index=int(len(x)/2)
    

    #print peak_index
    start_ind=peak_index-int(window_open/dt)
    stop_ind=peak_index+int(window_close/dt)
    start_avg=np.average(y[:n_avg])
    stop_avg=np.average(y[(len(y)-n_avg-1):])

    x_fit=np.concatenate((x[0:(peak_index-int(window_open/dt))], x[(peak_index+int(window_close/dt)):]))
    y_fit=np.concatenate((y[0:(peak_index-int(window_open/dt))], y[(peak_index+int(window_close/dt)):]))

    roll=np.zeros(len(y))
    #print peak_index
    total_window=len(x)
    pre=len(x[0:(peak_index-int(window_open/dt))])
    sig=len(x[(peak_index-int(window_open/dt)):(peak_index+int(window_close/dt))])
    post=len(x[(peak_index+int(window_close/dt)):len(x)])
    
    
    #print '{0}:  {1} + {2} + {3} = {4}'.format(total_window,pre,sig,post,pre+sig+post)
    
    # print '{1} : {2} : {3}'.format(x[0:(peak_index-int(window_open/dt))][0])
    
    for i in np.arange(len(y)):
        if i<n_avg:
            roll[i]=start_avg
        elif i>(len(y)-n_avg-1):
            roll[i]=stop_avg
        else:
            roll[i]=np.average(y[(i-int(n_avg/2)):(i+int(n_avg/2))])

    first_sig=np.average(y[(peak_index-int(window_open/dt)-n_avg-10):(peak_index-int(window_open/dt)-10)])
    last_sig=np.average(y[(peak_index+int(window_close/dt)+10):(peak_index+int(window_close/dt)+n_avg)+10])
    #print first_sig,last_sig


    for i in np.arange((peak_index-int(window_open/dt)-20),(peak_index+int(window_close/dt)+20)):
        roll[i]=np.average([first_sig,last_sig])
    
    
    
    
    return x,roll



def findPeaks(trace):
    # peaks, x=find_peaks(trace,distance=3.,height=6.0,prominence=10.)
    peaks, x=find_peaks(trace,distance=3.,height=6.0,prominence=30.)
    return peaks



def fit_gauss(counts,bin_start,bin_stop,nBins,min=0):
    
    counts=counts[counts>min]
    binWidth=(bin_stop-bin_start)/nBins
    bins=np.arange(bin_start,bin_stop,binWidth)
    Hist = TH1F("Hist", "Hist_x;x-axis;Frequency",int(nBins),int(bin_start),int(bin_stop))
    fit = TF1("fit", "gaus", bin_start,bin_stop)
    #fit = TF1("fit", "gaus", -200,10000)
    
    for i in np.arange(len(counts)):
        Hist.Fill(counts[i])

    Hist.Fit("fit",'0')

    p0=fit.GetParameter(0)
    p1=fit.GetParameter(1)
    p2=fit.GetParameter(2)

    x=np.arange(bin_start,bin_stop,1)
    y=np.zeros([len(x)])
    for i in np.arange(len(x)):
        y[i]=fit.Eval(x[i])

    return p0, p1, p2, x, y




def avg_baseline_correction(time,data):
    
    dt=time[1]-time[0]
    ldata=len(data)
    peak_index=np.argmax(np.abs(data))

    if (peak_index-int(constants.BG_MAX/dt)<=0) or peak_index+int(constants.BG_MAX/dt)>=(ldata-1):
        peak_index=int(ldata/2)
    
    
    avg=np.average(data[constants.BG_MIN:constants.BG_MAX])


    return data-avg


def rolling_correction(traces): # 2d-> time,data
    
    if traces.ndim == 2:
        traces = np.expand_dims(traces, axis=0)
    
    nTraces=len(traces)
    dt=traces[0][0][1]-traces[0][0][0]
    ldata=len(traces[0][0])

    for i in np.arange(nTraces):
            
        peak_index=np.argmax(np.abs(traces[i][1]))
        
        # check that we aren't too close to beginning or end
        
        if (peak_index-int(constants.WINDOW_OPEN/dt)<=0) or (peak_index+int(constants.WINDOW_CLOSE/dt))>=(ldata-1):
            peak_index=int(ldata/2)
        
        start_ind=peak_index-int(constants.WINDOW_OPEN/dt)-int(constants.WINDOW_PRE/dt)
        stop_ind=peak_index+-int(constants.WINDOW_OPEN/dt)

        avg=np.average(traces[i][1][start_ind:stop_ind])
        traces[i][1]=traces[i][1]-avg




    return np.squeeze(traces)


def counts(traces):

    if traces.ndim == 2:
        traces = np.expand_dims(traces, axis=0)
    
    nTraces=len(traces)
    counts=np.zeros([len(traces)])
    dt=traces[0][0][1]-traces[0][0][0]
    ldata=len(traces[0][0])
    for i in np.arange(nTraces):

        peak_index=np.argmax(np.abs(traces[i][1]))
    
    # check that we aren't too close to beginning or end

        if (peak_index-int(constants.WINDOW_OPEN/dt)<=0) or (peak_index+int(constants.WINDOW_CLOSE/dt))>=(ldata-1):
            peak_index=int(ldata/2)

        start_ind=peak_index-int(constants.WINDOW_OPEN/dt)
        stop_ind=peak_index+int(constants.WINDOW_CLOSE/dt)
        counts[i]=np.abs(np.sum(traces[i][1][start_ind:stop_ind]))
                

    return np.squeeze(counts)
