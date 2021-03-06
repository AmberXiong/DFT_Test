#!/usr/bin/env python
# -*- coding: utf-8 -*-

## @package DFT
# Analysis data with DFT.
#

import numpy
import numpy as np
import math
import pylab as pl
import matplotlib.pyplot as plt
import scipy.signal as signal
import seaborn as sns

## Calculate the power spectrum of array x.
# @param[in] xf is FFT result in frequency domain.
# @param[in] fft_size is the lenth of xf.
# @return xps_dBm is power value in frequency domain.
def Simple_DFT(x, fft_size):
    xs = x[:fft_size]              # reshape input array with fft_size
    xf = np.fft.rfft(xs)/fft_size  # DFT
    xps = 2*abs(xf)*abs(xf)/50

    if fft_size%2:                 # fft_size is even
        xps[0] /= 2
    else:                          # fft_size is odd
        xps[0] /= 2
        xps[int(fft_size/2)] /= 2
    xps_dBm = 10*np.log10(xps*1000)
    return [xps_dBm, numpy.sqrt(xps)]
    #return xf

## Multiply window function to array x.
# @param[in] x is a finite sequence of equally-spaced samples of a signal.
# @param[in] window is an array of window function.
# @param[in] fft_size is the lenth of window function and the effective lenth of x. 
# @return xw is the multiply result of x and window.
def Apply_Window(x, window, fft_size):
    xs = x[:fft_size]
    xw = xs*window
    return xw

## Add wihte noise to sampled signal x.
# @param[in] x is sampled signal array in time domain.
# @param[in] mean is mean value of white noise.
# @param[in] std is std value of white noise.
# @param[in] fft_size is the effective lenth of x.
# @return rtn is the add result of x and white noise.
def Add_White_Noise(x, mean, std, fft_size):
    white_noise = numpy.random.normal(mean, std, size=fft_size)
    xs = x[:fft_size]
    rtn = xs + white_noise
    return rtn

## Add quantization noise(uniform distributed) to sampled signal x.
# @param[in] x is sampled signal array in time domain.
# @param[in] v_fullscale is the fullscale voltage of ADC.
# @param[in] width is the bit number of ADC.
# @param[in] fft_size is the effective lenth of x.
# @return rtn is the output of ADC.
def Add_Qnt_Noise(x, v_fullscale, width, fft_size):
    Q = v_fullscale/numpy.power(2, width)
    xs = x[:fft_size]
    rtn = []
    for i in xs:
        n = i//Q
        p = (i/Q)-n
        if n > (numpy.power(2, width-1)):
            j = v_fullscale/2
        elif n < -(numpy.power(2, width-1)):
            j = -(v_fullscale/2)
        else:
            if p >= 0.5:
                j = (n+1)*Q
            else:
                j = n*Q
        rtn.append(j)
    return rtn

## Add phase noise to sample process.
# @param[in] w is value of phase item at time t.
# @param[in] mean is the mean value of phase noise.
# @param[in] std is the std value of phase noise.
# @param[in] fft_size is the sampling frequency.
# @return rtn is the add result of w and phase noise.
def Add_Phase_Noise(w, mean, std, fft_size):
    phase_noise = numpy.random.normal(mean, std, size=fft_size)
    ws = w[:fft_size]
    rtn = ws + phase_noise
    #print(phase_noise)
    #print(rtn)
    return rtn

## Divide a long array x into several fft_size arrays, calculate the average value of seperated FFT results.
# @param[in] x is long data array.
# @param[in] m is the number of divided array, m*fft_size <= len(x).
# @param[in] fft_size is the lenth of divided 
# @return avg_xps_dBm is the average value of power spectrum in frequency domain.
def Average_FFT(x, m, fft_size):
    xs  = x[:m].reshape(-1, fft_size)
    xf  = np.fft.rfft(xs)/fft_size
    xps = 2*abs(xf)*abs(xf)
    avg_xps = np.average(xps, axis=0)

    if fft_size%2:                 # fft_size is even
        avg_xps[0] /= 2
    else:                          # fft_size is odd
        avg_xps[0] /= 2
        avg_xps[int(fft_size/2)] /= 2
    avg_xps_dBm = 10*np.log10(avg_xps)
    return [avg_xps_dBm, avg_xps]

## Extract data from each frequency bin, count the numbers of the histogram bins, save data to RF_power_spectrum.dat.
# @param[in] x is a 2_D array, it saves multiple times FFT results.
# @param[in] freqs is an array of frequency bins.
# @param[in] loops is the number of FFT times.
# @param[in] lenth is the data lenth of each FFT result.
# @param[in] bins is number of bins of histogram.
# @param[in] output_file is a string, it's the name of output file.
def Data_Dist(x, freqs, loops, lenth, num_bins, output_file):
    D    = []
    for j in range(lenth):
        d = []
        for i in range(loops):
            d.append(x[i][j])
        D.append(d)

    cnt     = []
    value   = []
    for i in range(lenth):
        n, val, patches = plt.hist(D[i], bins=num_bins)
        cnt.append(n)
        value.append(val)
    
    FREQ = []
    VAL  = []
    CNT  = []
    DATA = []
    for i in range(lenth):
        for j in range(num_bins):
            FREQ.append(freqs[i])
            VAL.append(value[i][j])
            CNT.append(cnt[i][j])
            DATA.append([freqs[i],value[i][j],cnt[i][j]])
    numpy.savetxt(output_file ,DATA)    

## Calculate the average value of x every avrg_num times.
# @param[in] x is an array.
# @param[in] avrg_num is the data lenth for each averaging.
# @param[in] times is the number of averaging, len(x)>=avrg_num*times.
# @param[in] num_bins is the number of bins for histogram.
# @param[in] output_file is the name of output file.
# @return an array which stores the average value of x, len(X_avrg)=times.
def Avrg_Amp(x, avrg_num, times, num_bins, output_file):
    X_avrg = []
    for n in range(times):
        x_sum = 0
        for m in range(avrg_num):
            x_sum += x[n*avrg_num+m]
        x_avrg = x_sum/avrg_num
        X_avrg.append(x_avrg)
    n, val, patches = plt.hist(X_avrg, bins=num_bins)
    
    DATA = []
    for i in range(num_bins):
        DATA.append([val[i], n[i]])
    numpy.savetxt(output_file, DATA)
    return X_avrg

## Plot the power spectrum with/without noise.
# @param[in] xps1 is DFT result of given signal.
# @param[in] xps2 is DFT result with noise.
# @param[in] freqs is an array of frequency.
# @param[in] signal_amp is the amplitude Vpp of the given signal.
# @param[in] sampling_rate and fft_size determine the frequency resolution.
# @param[in] noise_type is a tring, indecate the tyoe of noise we applied to signal.
# @param[in] m1 is mean value or lower boundary of noise.
# @param[in] m2 is std value or higher boundary of noise.
# @return a figure of power spectrum.
def Plot_Multi_PS(xps1, xps2, freqs, signal_amp, sampling_rate, fft_size, noise_type, m1, m2): 
    xps1_dBm    = xps1[0]
    xps1_watt  = xps1[1]
    xps2_dBm    = xps2[0]
    xps2_watt  = xps2[1] 
    a1 = max(xps1_dBm)
    b1 = numpy.argmax(xps1_dBm)
    a2 = max(xps2_dBm)
    b2 = numpy.argmax(xps2_dBm)
    signal_power = signal_amp*signal_amp/2
    accuracy1 = abs(xps1_watt[b1]-signal_power)/signal_power
    accuracy2 = abs(xps2_watt[b2]-signal_power)/signal_power
    print(a2,b2)
    
    pl.figure(figsize=(8,6))
    pl.subplot(211)
    pl.plot(freqs, xps1_dBm, label=u"Without noise")
    pl.ylabel(u'Power(dBm)')
    pl.annotate('%s Hz, %s dBm'%(b1*sampling_rate/fft_size, a1),xy=(freqs[b1], xps1_dBm[b1]), xytext=(freqs[b1]+5, xps1_dBm[b1]-30))
    props = dict(boxstyle='round', facecolor='none', alpha=0.5)
    pl.text(1500, -200, 'accuracy:%s'%accuracy1, size=10, bbox=props)
    pl.legend(fontsize=8)
    
    pl.title(u"Power Spectrum(fs=%sHz, n=%s)"%(sampling_rate,fft_size))
    pl.subplot(212)
    pl.plot(freqs, xps2_dBm, color='green',label=u"With %s"%(noise_type))
    pl.xlabel(u'Frequency(Hz)')
    pl.ylabel(u'Power(dBm)')
    pl.annotate('%s Hz, %s dBm'%(b2*sampling_rate/fft_size, a2),xy=(freqs[b2], xps2_dBm[b2]), xytext=(freqs[b2]+5, xps2_dBm[b2]-10))
    props = dict(boxstyle='round', facecolor='none', alpha=0.5)
    pl.text(1500, -80, 'mean/lower-boundary:%s\nstd/higher-boundary:%s\naccuracy:%s'%(m1,m2,accuracy2), size=10, bbox=props)
    pl.legend(fontsize=8)
    pl.savefig('power_spectrum.pdf')

## Plot the power spectrum with/without window.
# @param[in] xps1 is DFT result of the signal without noise window.
# @param[in] xps1_w is DFT result of signal only with window.
# @param[in] xps2 is DFT result of signal only with noise.
# @param[in] xps2_w is DFT result of signal both with noise and window.
# @param[in] freqs is an array of frequency.
# @param[in] signal_amp is the amplitude of given signal.
# @parma[in] sampling_rate and fft_size determine the frequency resolution.
# @param[in] noise_type is a string, indecate the tyoe of noise we applied to signal.
# @param[in] window_type is a string, indecate the type of window we applied to signal.
# @param[in] m1 is mean value or lower boundary of noise.
# @param[in] m2 is std value or higher boundary of noise.
# @return a figure of power spectrum.
def Plot_PS_W(xps1, xps1_w, xps2, xps2_w, freqs, signal_amp, sampling_rate, fft_size, noise_type, window_type, m1, m2):
    xps1_dBm     = xps1[0]
    xps1_w_dBm   = xps1_w[0]
    xps2_dBm     = xps2[0]
    xps2_w_dBm   = xps2_w[0]
    
    xps1_watt   = xps1[1]
    xps1_w_watt = xps1_w[1]
    xps2_watt   = xps2[1]
    xps2_w_watt = xps2_w[1]

    a1   = max(xps1_dBm)
    a1_w = max(xps1_w_dBm)
    b1   = numpy.argmax(xps1_dBm)
    b1_w = numpy.argmax(xps1_w_dBm)
    a2   = max(xps2_dBm)
    a2_w = max(xps2_w_dBm)
    b2   = numpy.argmax(xps2_dBm)
    b2_w = numpy.argmax(xps2_w_dBm)

    signal_power = signal_amp*signal_amp/2
    accuracy1   = abs(xps1_watt[b1]-signal_power)/signal_power
    accuracy1_w = abs(xps1_w_watt[b1_w]-signal_power)/signal_power
    accuracy2   = abs(xps2_watt[b2]-signal_power)/signal_power
    accuracy2_w = abs(xps2_w_watt[b2_w]-signal_power)/signal_power
    print(a2,b2)

    pl.figure(figsize=(8,6))
    pl.subplot(211)
    pl.step(freqs, xps1_dBm, color='green', label=u"Without noise, without window")
    pl.step(freqs, xps1_w_dBm, label=u'Without noise, with %s window'%(window_type))
    pl.ylabel(u'Power(dBm)')
    pl.annotate('%s Hz, %s dBm'%(b1*sampling_rate/fft_size, a1),xy=(freqs[b1], xps1_dBm[b1]), xytext=(freqs[b1]+5, xps1_dBm[b1]-30))
    props = dict(boxstyle='round', facecolor='none', alpha=0.5)
    pl.text(1500, -120, 'accuracy:%s\naccuracy_w:%s'%(accuracy1,accuracy1_w), size=10, bbox=props)
    pl.legend(fontsize=8)

    pl.title(u"Power Spectrum(fs=%sHz, n=%s)"%(sampling_rate,fft_size))
    pl.subplot(212)
    pl.step(freqs, xps2_dBm, color='green',label=u"With %s, without window"%(noise_type))
    pl.step(freqs, xps2_w_dBm, label=u'With %s, with %s window'%(noise_type, window_type))
    pl.xlabel(u'Frequency(Hz)')
    pl.ylabel(u'Power(dBm)')
    pl.annotate('%s Hz, %s dBm'%(b2*sampling_rate/fft_size, a2),xy=(freqs[b2], xps2_dBm[b2]), xytext=(freqs[b2]+5, xps2_dBm[b2]-10))
    props = dict(boxstyle='round', facecolor='none', alpha=0.5)
    pl.text(1500, -100, 'mean/lower-boundary:%s\nstd/higher-boundary:%s\naccuracy:%s\naccuracy_w:%s'%(m1,m2,accuracy2,accuracy2_w), size=10, bbox=props)
    pl.legend(fontsize=8)
    pl.savefig('power_spectrum_wind.pdf')

if __name__ == "__main__":
    sampling_rate = 400
    fft_size      = 512
    signal_freq   = 130
    signal_amp    = 1
    wht_mean      = 0
    wht_std       = 0.00005
    phs_mean      = 0
    phs_std       = 0.00001
    v_fs          = 2
    width         = 10
    Q             = v_fs/numpy.power(2, width)
    #mean_num      = 20  # Take an average of every mean_num numbers
    #times         = 200  # number of results
    window_num    = 8000  # number of window
    loops         = 500
    lenth         = (fft_size+2)//2  # data lenth of FFT results 
    num_bins      = 50  # bin number for full frequency scale
    bins          = 10  # bin number for signal amplitude

    t     = np.arange(0, fft_size/sampling_rate, 1.0/sampling_rate)  #step is out of range
    print(len(t))
    freqs = np.linspace(0, sampling_rate/2, fft_size/2+1)
    w     = 2*np.pi*signal_freq*t

    Avrg_Y     = []
    Avrg_Y_wht = []
    Avrg_Y_phs = []
    Avrg_Y_qnt = []
    Avrg_Y_all = []

    Avrg_Y_dB     = []
    Avrg_Y_wht_dB = []
    Avrg_Y_phs_dB = []
    Avrg_Y_qnt_dB = []
    Avrg_Y_all_dB = []
    for n in range(loops):
        #Y          = []
        #Y_wht      = []
        #Y_phs      = []
        #Y_qnt      = []
        #Y_all      = []
        Y_w        = []
        Y_wht_w    = []
        Y_phs_w    = []
        Y_qnt_w    = []
        Y_all_w    = []

        #Y_dB       = []
        #Y_wht_dB   = []
        #Y_phs_dB   = []
        #Y_qnt_dB   = []
        #Y_all_dB   = []
        Y_w_dB     = []
        Y_wht_w_dB = []
        Y_phs_w_dB = []
        Y_qnt_w_dB = []
        Y_all_w_dB = []
        # repeat m times, for white noise and phase noise
        for i in range(window_num):
            x     = signal_amp*np.sin(w)  #sample the points near the +-1
            x_wht = Add_White_Noise(x, wht_mean, wht_std, fft_size)
            w_phs = Add_Phase_Noise(w, phs_mean, phs_std, fft_size)
            x_phs = signal_amp*np.sin(w_phs)

            x_p   = signal_amp*np.sin(w+i*np.pi/window_num+n*np.pi/loops)
            x_qnt = Add_Qnt_Noise(x_p, v_fs, width, fft_size)
            
            x_phs_wht = Add_White_Noise(signal_amp*np.sin(w_phs+i*np.pi/window_num+n*np.pi/loops), wht_mean, wht_std, fft_size)
            x_all = Add_Qnt_Noise(x_phs_wht, v_fs, width, fft_size) 
            
            #window  = signal.hamming(fft_size, sym=0)
            window  = signal.hann(fft_size, sym=0)
            #window  = np.array([1 for i in range(fft_size)])
            #window  = signal.bartlett(fft_size, sym=0)

            #x_w     = Apply_Window(x, window, fft_size)
            #x_wht_w = Apply_Window(x_wht, window, fft_size)
            #x_phs_w = Apply_Window(x_phs, window, fft_size)
            x_qnt_w = Apply_Window(x_qnt, window, fft_size)
            x_all_w = Apply_Window(x_all, window, fft_size)

            #y       = Simple_DFT(x, fft_size)
            #y_wht   = Simple_DFT(x_wht, fft_size)
            #y_phs   = Simple_DFT(x_phs, fft_size)
            #y_qnt   = Simple_DFT(x_qnt, fft_size)
            #y_all   = Simple_DFT(x_all, fft_size)
            #y_w     = Simple_DFT(x_w, fft_size)
            #y_wht_w = Simple_DFT(x_wht_w, fft_size)
            #y_phs_w = Simple_DFT(x_phs_w, fft_size)
            y_qnt_w = Simple_DFT(x_qnt_w, fft_size)
            y_all_w = Simple_DFT(x_all_w, fft_size)
            
            # power
            #Y.append(y[1])
            #Y_wht.append(y_wht[1])
            #Y_phs.append(y_phs[1])
            #Y_qnt.append(y_qnt[1])
            #Y_all.append(y_all[1])

            #Y_w.append(y_w[1])
            #Y_wht_w.append(y_wht_w[1])
            #Y_phs_w.append(y_phs_w[1])
            Y_qnt_w.append(y_qnt_w[1])
            Y_all_w.append(y_all_w[1])

            freq_y = numpy.argmax(y_all_w[1])*sampling_rate/fft_size
            
            #Y_dB.append(y[0])
            #Y_wht_dB.append(y_wht[0])
            #Y_phs_dB.append(y_phs[0])
            #Y_qnt_dB.append(y_qnt[0])
            #Y_all_dB.append(y_all[0])

            #Y_w_dB.append(y_w[0])
            #Y_wht_w_dB.append(y_wht_w[0])
            #Y_phs_w_dB.append(y_phs_w[0])
            Y_qnt_w_dB.append(y_qnt_w[0])
            Y_all_w_dB.append(y_all_w[0])
        
        #Y_w_trans = numpy.transpose(Y_w)
        #Y_wht_w_trans = numpy.transpose(Y_wht_w)
        #Y_phs_w_trans = numpy.transpose(Y_phs_w)
        Y_qnt_w_trans = numpy.transpose(Y_qnt_w)
        Y_all_w_trans = numpy.transpose(Y_all_w)

        #Y_w_dB_trans     = numpy.transpose(Y_w_dB)
        #Y_wht_w_dB_trans = numpy.transpose(Y_wht_w_dB)
        #Y_phs_w_dB_trans = numpy.transpose(Y_phs_w_dB)
        Y_qnt_w_dB_trans = numpy.transpose(Y_qnt_w_dB)
        Y_all_w_dB_trans = numpy.transpose(Y_all_w_dB)

        #single_Y_w = []
        #single_Y_wht_w = []
        #single_Y_phs_w = []
        single_Y_qnt_w = []
        single_Y_all_w = []

        #single_w_dB = []
        #single_wht_w_dB = []
        #single_phs_w_dB = []
        single_qnt_w_dB = []
        single_all_w_dB = []

        for k in range(lenth):
            #s_w     = numpy.mean(Y_w_trans[k])
            #s_wht_w = numpy.mean(Y_wht_w_trans[k])
            #s_phs_w = numpy.mean(Y_phs_w_trans[k])
            s_qnt_w = numpy.mean(Y_qnt_w_trans[k])
            s_all_w = numpy.mean(Y_all_w_trans[k])

            #s_w_dB     = numpy.mean(Y_w_dB_trans[k])
            #s_wht_w_dB = numpy.mean(Y_wht_w_dB_trans[k])
            #s_phs_w_dB = numpy.mean(Y_phs_w_dB_trans[k])
            s_qnt_w_dB = numpy.mean(Y_qnt_w_dB_trans[k])
            s_all_w_dB = numpy.mean(Y_all_w_dB_trans[k])

            #single_Y_w.append(s_w)
            #single_Y_wht_w.append(s_wht_w)
            #single_Y_phs_w.append(s_phs_w)
            single_Y_qnt_w.append(s_qnt_w)
            single_Y_all_w.append(s_all_w)

            #single_w_dB.append(s_w_dB)
            #single_wht_w_dB.append(s_wht_w_dB)
            #single_phs_w_dB.append(s_phs_w_dB)
            single_qnt_w_dB.append(s_qnt_w_dB)
            single_all_w_dB.append(s_all_w_dB)
        
        print(len(single_Y_qnt_w)) 
        #average amplitude for each frequency bin
        #Avrg_Y.append(max(single_Y_w))
        #Avrg_Y_wht.append(max(single_Y_wht_w))
        #Avrg_Y_phs.append(max(single_Y_phs_w))
        Avrg_Y_qnt.append(max(single_Y_qnt_w))
        Avrg_Y_all.append(max(single_Y_all_w))

        #Avrg_Y_dB.append(single_w_dB)
        #Avrg_Y_wht_dB.append(single_wht_w_dB)
        #Avrg_Y_phs_dB.append(single_phs_w_dB)
        Avrg_Y_qnt_dB.append(single_qnt_w_dB)
        Avrg_Y_all_dB.append(single_all_w_dB)
    #A       = numpy.sqrt(Y)  #Amplitude of signal
    #A_wht   = numpy.sqrt(Y_wht)
    #A_wht_w = numpy.sqrt(Y_wht_w)
    #A_qnt   = numpy.sqrt(Y_qnt)
    #A_qnt_w = numpy.sqrt(Y_qnt_w)
    #A_phs   = numpy.sqrt(Y_phs)
    #A_phs_w = numpy.sqrt(Y_phs_w)
    #A_all   = numpy.sqrt(Y_all)
    #A_all_w = numpy.sqrt(Y_all_w)
        
    #Y_wht_w_mean = Avrg_Amp(A_wht_w, mean_num, times, bins, "Average_amp_dist_y_wht_w.csv")
    #Y_qnt_w_mean = Avrg_Amp(A_qnt_w, mean_num, times, bins, "Average_amp_dist_y_qnt_w.csv")
    #Y_all_w_mean = Avrg_Amp(A_all_w, mean_num, times, bins, "all_hann_160.csv")
    #Y_phs_w_mean = Avrg_Amp(A_phs_w, mean_num, times, bins, "Average_amp_dist_y_phs_w.csv")

    #D_wht_w = Data_Dist(Y_wht_w_dB, freqs, loops, lenth, num_bins, "RF_wht_window.csv")
    D_qnt_w = Data_Dist(Y_qnt_w_dB, freqs, loops, lenth, num_bins, "RF_qnt_window_dB.csv")
    D_all_w = Data_Dist(Avrg_Y_all_dB, freqs, loops, lenth, num_bins, "RF_all_window_dB.csv")
    #D_phs_w = Data_Dist(Y_phs_w_dB, freqs, loops, lenth, num_bins, "RF_phs_window_dB.csv")

    #print("std_A_wht_w:",numpy.std(Avrg_Y_wht),"mean_A_wht_w:",numpy.mean(Avrg_Y_wht))
    #print("std_A_wht_w_mean:",numpy.std(Y_wht_w_mean),"mean_A_wht_w_mean:",numpy.mean(Y_wht_w_mean))
    print("std_A_qnt_w:",numpy.std(Avrg_Y_qnt),"mean_A_qnt_w:",numpy.mean(Avrg_Y_qnt))
    #print("std_A_qnt_w_mean:",numpy.std(Y_qnt_w_mean),"mean_A_qnt_w_mean:",numpy.mean(Y_qnt_w_mean))
    #print("std_A_phs_w:",numpy.std(Avrg_Y_phs),"mean_A_phs_w:",numpy.mean(Avrg_Y_phs))
    #print("std_A_phs_w_mean:",numpy.std(Y_phs_w_mean),"mean_A_phs_w_mean:",numpy.mean(Y_phs_w_mean))
    print("std_A_all_w:",numpy.std(Avrg_Y_all),"mean_A_all_w:",numpy.mean(Avrg_Y_all))
    #print("std_A_all_w_mean:",numpy.std(Y_all_w_mean),"mean_A_all_w_mean:",numpy.mean(Y_all_w_mean))
    
    #plt.hist(Y_qnt_mean, bins)
    #plt.show()
