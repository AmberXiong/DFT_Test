#!/usr/bin/env python
# -*- coding: utf-8 -*-

## @package DFT
# Analysis data with DFT.
#

import numpy
import numpy as np
import math
import pylab as pl
import scipy.signal as signal
import seaborn as sns

## Calculate the power spectrum of array x.
# @param[in] xf is FFT result in frequency domain.
# @param[in] fft_size is the lenth of xf.
# @return xps_dB is power value in frequency domain.
def Simple_DFT(x, fft_size):
    xs = x[:fft_size]              # reshape input array with fft_size
    xf = np.fft.rfft(xs)/fft_size  # DFT
    xps = 2*abs(xf)*abs(xf)

    if fft_size%2:                 # fft_size is even
        xps[0] /= 2
    else:                          # fft_size is odd
        xps[0] /= 2
        xps[int(fft_size/2)] /= 2
    xps_dB = 10*np.log10(xps)
    return [xps_dB, xps]
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
# @param[in] low is lower boundary of quantization noise.
# @param[in] high is higher boundary of quantization noise.
# @param[in] fft_size is the effective lenth of x.
# @return rtn is the add result of x and quantization noise.
def Add_Qnt_Noise(x, v_fullscale, width, fft_size):
    Q = v_fullscale/numpy.power(2, width)
    xs = x[:fft_size]
    rtn = []
    for i in xs:
        n = i//Q
        if n>=numpy.power(2,width-1):
            j = (numpy.power(2, width-1)-0.5)*Q
        elif n<(-numpy.power(2, width-1)):
            j = (0.5-numpy.power(2, width-1))*Q
        else:
            j = (n+0.5)*Q
        rtn.append(j)
    return rtn

## Add phase noise to sample process.
# @param[in] w is value of phase item at time t.
# @param[in] mean is the mean value of phase noise.
# @param[in] std is the std value of phase noise.
# @param[in] sampling_rate is the sampling frequency.
# @return rtn is the add result of w and phase noise.
def Add_Phase_Noise(w, mean, std, sampling_rate):
    phase_noise = numpy.random.normal(mean, std, size=sampling_rate)
    ws = w[:sampling_rate]
    rtn = ws + phase_noise
    #print(phase_noise)
    #print(rtn)
    return rtn

## Divide a long array x into several fft_size arrays, calculate the average value of seperated FFT results.
# @param[in] x is long data array.
# @param[in] m is the number of divided array, m*fft_size <= len(x).
# @param[in] fft_size is the lenth of divided 
# @return avg_xps_dB is the average value of power spectrum in frequency domain.
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
    avg_xps_dB = 10*np.log10(avg_xps)
    return [avg_xps_dB, avg_xps]

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
    xps1_dB    = xps1[0]
    xps1_watt  = xps1[1]
    xps2_dB    = xps2[0]
    xps2_watt  = xps2[1] 
    a1 = max(xps1_dB)
    b1 = numpy.argmax(xps1_dB)
    a2 = max(xps2_dB)
    b2 = numpy.argmax(xps2_dB)
    signal_power = signal_amp*signal_amp/2
    accuracy1 = abs(xps1_watt[b1]-signal_power)/signal_power
    accuracy2 = abs(xps2_watt[b2]-signal_power)/signal_power
    print(a2,b2)
    
    pl.figure(figsize=(8,6))
    pl.subplot(211)
    pl.plot(freqs, xps1_dB, label=u"Without noise")
    pl.ylabel(u'Power(dB)')
    pl.annotate('%s Hz, %s dB'%(b1*sampling_rate/fft_size, a1),xy=(freqs[b1], xps1_dB[b1]), xytext=(freqs[b1]+5, xps1_dB[b1]-30))
    props = dict(boxstyle='round', facecolor='none', alpha=0.5)
    pl.text(1500, -200, 'accuracy:%s'%accuracy1, size=10, bbox=props)
    pl.legend(fontsize=8)
    
    pl.title(u"Power Spectrum(fs=%sHz, n=%s)"%(sampling_rate,fft_size))
    pl.subplot(212)
    pl.plot(freqs, xps2_dB, color='green',label=u"With %s"%(noise_type))
    pl.xlabel(u'Frequency(Hz)')
    pl.ylabel(u'Power(dB)')
    pl.annotate('%s Hz, %s dB'%(b2*sampling_rate/fft_size, a2),xy=(freqs[b2], xps2_dB[b2]), xytext=(freqs[b2]+5, xps2_dB[b2]-10))
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
    xps1_dB     = xps1[0]
    xps1_w_dB   = xps1_w[0]
    xps2_dB     = xps2[0]
    xps2_w_dB   = xps2_w[0]
    
    xps1_watt   = xps1[1]
    xps1_w_watt = xps1_w[1]
    xps2_watt   = xps2[1]
    xps2_w_watt = xps2_w[1]

    a1   = max(xps1_dB)
    a1_w = max(xps1_w_dB)
    b1   = numpy.argmax(xps1_dB)
    b1_w = numpy.argmax(xps1_w_dB)
    a2   = max(xps2_dB)
    a2_w = max(xps2_w_dB)
    b2   = numpy.argmax(xps2_dB)
    b2_w = numpy.argmax(xps2_w_dB)

    signal_power = signal_amp*signal_amp/2
    accuracy1   = abs(xps1_watt[b1]-signal_power)/signal_power
    accuracy1_w = abs(xps1_w_watt[b1_w]-signal_power)/signal_power
    accuracy2   = abs(xps2_watt[b2]-signal_power)/signal_power
    accuracy2_w = abs(xps2_w_watt[b2_w]-signal_power)/signal_power
    print(a2,b2)

    pl.figure(figsize=(8,6))
    pl.subplot(211)
    pl.step(freqs, xps1_dB, color='green', label=u"Without noise, without window")
    pl.step(freqs, xps1_w_dB, label=u'Without noise, with %s window'%(window_type))
    pl.ylabel(u'Power(dB)')
    pl.annotate('%s Hz, %s dB'%(b1*sampling_rate/fft_size, a1),xy=(freqs[b1], xps1_dB[b1]), xytext=(freqs[b1]+5, xps1_dB[b1]-30))
    props = dict(boxstyle='round', facecolor='none', alpha=0.5)
    pl.text(1500, -120, 'accuracy:%s\naccuracy_w:%s'%(accuracy1,accuracy1_w), size=10, bbox=props)
    pl.legend(fontsize=8)

    pl.title(u"Power Spectrum(fs=%sHz, n=%s)"%(sampling_rate,fft_size))
    pl.subplot(212)
    pl.step(freqs, xps2_dB, color='green',label=u"With %s, without window"%(noise_type))
    pl.step(freqs, xps2_w_dB, label=u'With %s, with %s window'%(noise_type, window_type))
    pl.xlabel(u'Frequency(Hz)')
    pl.ylabel(u'Power(dB)')
    pl.annotate('%s Hz, %s dB'%(b2*sampling_rate/fft_size, a2),xy=(freqs[b2], xps2_dB[b2]), xytext=(freqs[b2]+5, xps2_dB[b2]-10))
    props = dict(boxstyle='round', facecolor='none', alpha=0.5)
    pl.text(1500, -100, 'mean/lower-boundary:%s\nstd/higher-boundary:%s\naccuracy:%s\naccuracy_w:%s'%(m1,m2,accuracy2,accuracy2_w), size=10, bbox=props)
    pl.legend(fontsize=8)
    pl.savefig('power_spectrum_wind.pdf')

if __name__ == "__main__":
    sampling_rate = 5120
    fft_size      = 512
    signal_freq   = 200
    signal_amp    = 1
    wht_mean      = 0
    wht_std       = 0.00005
    phs_mean      = 0
    phs_std       = 0.0001
    v_fs          = 2
    width         = 16
    Q             = v_fs/numpy.power(2, width)
    loops             = 500

    t     = np.arange(0, 1.0, 1.0/sampling_rate)
    freqs = np.linspace(0, sampling_rate/2, fft_size/2+1)
    w     = 2*np.pi*signal_freq*t

    Y       = []
    Y_wht   = []
    Y_phs   = []
    Y_qnt   = []
    Y_w     = []
    Y_wht_w = []
    Y_phs_w = []
    Y_qnt_w = []
# repeat m times, for white noise and phase noise
    for i in range(loops):
        x     = signal_amp*np.sin(w)
        x_wht = Add_White_Noise(x, wht_mean, wht_std, fft_size)
        w_phs = Add_Phase_Noise(w, phs_mean, phs_std, sampling_rate)
        x_phs = signal_amp*np.sin(w_phs)

        x_p   = signal_amp*np.sin(w+i*2*np.pi/loops+1)
        x_qnt = Add_Qnt_Noise(x_p, v_fs, width, fft_size)

        window  = signal.hann(fft_size, sym=0)
        x_w     = Apply_Window(x, window, fft_size)
        x_wht_w = Apply_Window(x_wht, window, fft_size)
        x_phs_w = Apply_Window(x_phs, window, fft_size)
        x_qnt_w = Apply_Window(x_qnt, window, fft_size)

        y       = Simple_DFT(x, fft_size)
        y_wht   = Simple_DFT(x_wht, fft_size)
        y_phs   = Simple_DFT(x_phs, fft_size)
        y_qnt   = Simple_DFT(x_qnt, fft_size)
        y_w     = Simple_DFT(x_w, fft_size)
        y_wht_w = Simple_DFT(x_wht_w, fft_size)
        y_phs_w = Simple_DFT(x_phs_w, fft_size)
        y_qnt_w = Simple_DFT(x_qnt_w, fft_size)

        Y.append(max(y[1]))
        Y_wht.append(max(y_wht[1]))
        Y_phs.append(max(y_phs[1]))
        Y_qnt.append(max(y_qnt[1]))
        Y_w.append(max(y_w[1]))
        Y_wht_w.append(max(y_wht_w[1]))
        Y_phs_w.append(max(y_phs_w[1]))
        Y_qnt_w.append(max(y_qnt_w[1]))

    print("std_Y_wht:",numpy.std(Y_wht),"len_Y_wht:",len(Y_wht))
    print("std_Y_wht_w:",numpy.std(Y_wht_w))
    print("std_Y_qnt:",numpy.std(Y_qnt))
    print("std_Y_qnt_w:",numpy.std(Y_qnt_w))
    print("std_Y_phs_w:",numpy.std(Y_phs_w),"std_Y_phs:",numpy.std(Y_phs))
    print("mean_Y_qnt:",numpy.mean(Y_qnt))
    print("mean_Y_wht:",numpy.mean(Y_wht),"mean_Y_wht_w:",numpy.mean(Y_wht_w))
 
    #sns.distplot(Y) 
    #sns.distplot(Y_wht,label="with white noise")
    #sns.distplot(Y_qnt,label="with quantization noise")
    #sns.distplot(Y_phs,label="with phase noise")
    #sns.distplot(Y-w)
    #sns.distplot(Y_wht_w,label="with white noise and window")
    #sns.distplot(Y_qnt_w)
    sns.distplot(Y_phs_w)
    sns.plt.show()
