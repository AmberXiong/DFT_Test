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
def Add_Qnt_Noise(x, low, high, fft_size):
    qnt_noise = numpy.random.uniform(low, high, fft_size)
    xs = x[:fft_size]
    rtn = xs + qnt_noise
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

## Plot the power spectrum.
# @param[in] xps1 is DFT result of given signal.
# @param[in] xps2 is DFT result with noise.
# @param[in] freqs is an array of frequency.
# @param[in] sampling_rate and fft_size determine the frequency resolution.
# @param[in] noise_type is a tring, indecate the tyoe of noise we applied to signal.
# @param[in] m1 is mean value or lower boundary of noise.
# @param[in] m2 is std value or higher boundary of noise.
# @return a figure of power spectrum.
def Plot_Multi_PS(xps1, xps2, freqs, sampling_rate, fft_size, noise_type, m1, m2): 
    xps1_dB = xps1[0]
    xps1_w  = xps1[1]
    xps2_dB = xps2[0]
    xps2_w  = xps2[1] 
    a1 = max(xps1_dB)
    b1 = numpy.argmax(xps1_dB)
    a2 = max(xps2_dB)
    b2 = numpy.argmax(xps2_dB)
    accuracy = abs(xps2_w[b2]-xps1_w[b2])/xps1_w[b2]
    print(a2,b2)
    
    pl.figure(figsize=(8,6))
    pl.subplot(211)
    pl.plot(freqs, xps1_dB, label=u"Without noise")
    pl.ylabel(u'Power(dB)')
    pl.annotate('%s Hz, %s dB'%(b1*sampling_rate/fft_size, a1),xy=(freqs[b1], xps1_dB[b1]), xytext=(freqs[b1]+5, xps1_dB[b1]-30))
    pl.legend()
    
    pl.title(u"Power Spectrum")
    pl.subplot(212)
    pl.plot(freqs, xps2_dB, color='green',label=u"With %s"%(noise_type))
    pl.xlabel(u'Frequency(Hz)')
    pl.ylabel(u'Power(dB)')
    pl.annotate('%s Hz, %s dB'%(b2*sampling_rate/fft_size, a2),xy=(freqs[b2], xps2_dB[b2]), xytext=(freqs[b2]+5, xps2_dB[b2]-10))
    props = dict(boxstyle='round', facecolor='none', alpha=0.5)
    pl.text(1500, -80, 'mean/lower-boundary:%s\nstd/higher-boundary:%s\naccuracy:%s'%(m1,m2,accuracy), size=10, bbox=props)
    pl.legend()
    pl.savefig('power_spectrum.pdf')

if __name__ == "__main__":


    sampling_rate = 5120
    fft_size      = 512
    signal_freq   = 200
    wht_mean      = 0.5
    wht_std       = 0.1
    phs_mean      = 0
    phs_std       = 0.001
    Q = 1/numpy.power(2, 16)

    t = np.arange(0, 1.0, 1.0/sampling_rate)
    freqs = np.linspace(0, sampling_rate/2, fft_size/2+1)
    w = 2*np.pi*signal_freq*t
    
    x = 1*np.sin(w)
    x_wht = Add_White_Noise(x, wht_mean, wht_std, fft_size)
    x_qnt = Add_Qnt_Noise(x, -Q/2, Q/2, fft_size)
    w_phs = Add_Phase_Noise(w, phs_mean, phs_std, sampling_rate)
    x_phs = 1*np.sin(w_phs)
    
    y = Simple_DFT(x, fft_size)
    y_wht = Simple_DFT(x_wht, fft_size)
    y_qnt = Simple_DFT(x_qnt, fft_size)
    y_phs = Simple_DFT(x_phs, fft_size)
    
    print(y[1][20], y_wht[1][20], y_qnt[1][20], y_phs[1][20])  
    print(y[0][20], y_wht[0][20], y_qnt[0][20], y_phs[0][20])
    
    #Plot_Multi_PS(y, y_wht, freqs, 5120, 512, 'white noise', wht_mean, wht_std)
    #Plot_Multi_PS(y, y_qnt, freqs, 5120, 512, 'quantization noise', -Q/2, Q/2)
    Plot_Multi_PS(y, y_phs, freqs, 5120, 512, 'phase noise', phs_mean, phs_std)


