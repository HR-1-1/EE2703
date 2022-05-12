"""
Title			: dft2.py
Author			: Harish R EE20B044
Created			: May 12 2022
Last Modified	: May 12 2022
Purpose			: Understand Spectra of non-periodic signals
Inputs			: (:
"""

from __future__ import division
import numpy.fft as nf
import numpy.random as nr
import numpy as np
import matplotlib.pyplot as plt

parameters = {'axes.labelsize': 10, 
              'axes.titlesize': 12, 
              'legend.fontsize': 8,  
             'mathtext.fontset':'cm',
             'figure.figsize': (10,15)} 

PATH = './plots/'

def f1(t,w=0.86):
	return np.cos(w*t)**3

def f2(t,w=1.5,d=0.5):
	return np.cos(w*t+d)

def f3(t):
	return (np.sin(t))**3

def f4(t):
	return (np.cos(t))**3

def f5(t):
	return np.cos(20*t+ 5*np.cos(t))

def gs(t):
	return np.exp(-1*t**2/2)

def gs_true(w):
	return np.exp(-(w**2)/2)/np.sqrt(2*np.pi)

def fft_acc(N):
	x = nr.rand(N)
	return np.abs(x-nf.ifft(nf.fft(x))).max()

def estim(w, Y, s=1e-4, wd=1):
	ii = np.where(w>0)
	w0 = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2))
	

def spectrum(xl, N, f, xlim, title, name,t0, tlim = False, window= False):
    
	if(tlim):
        t = t0
    else:
        t = np.linspace(-xl,xl,N+1)[:-1]
    
    fmax = 1/(t[1]-t[0])
	y = f(t)
    
	if (window):
        y=y*np.fftshift(0.54+0.46*cos(2*np.pi*np.arange(N)/N))
    
	y[0]=0 # the sample corresponding to -tmax should be set zeroo
    y = np.fftshift(y) # make y start with y(t=0)
    Y = np.fftshift(np.fft(y))/float(n)
    w = np.linspace(-np.pi*fmax, np.pi*fmax, n+1)[:-1]
	
	fig, ax = plt.subplots(2,1)
	ax[0].plot(w, abs(Y), lw=2)
	ax[0].set_xlim([-xlim, xlim])
	ax[0].set_ylabel(r'$|Y|$')
	ax[0].set_title(r'Spectrum of ' + title)
	ax[0].grid(True)
	
    phi = np.angle(Y)
	phi[np.where(abs(Y)<3e-3)] = 0
	
	ax[1].plot(w, phi, 'ro', lw=2)
	ax[1].set_xlim([-xlim, xlim])
	ax[1].set_ylabel(r'Phase of $Y$')
	ax[1].set_xlabel(r'$\omega$')
	ax[1].grid(True)
	
	fig.savefig(PATH+name)
    
	return w,Y


def main():
	
	print("FFT Max accuracy {}".format(fft_acc(128)))
	calc_fft(0, 2*np.pi, 128, f1, 1e-3, 10, '$\sin(5t)$', 'Fig1.png') 
	calc_fft(-4*np.pi, 4*np.pi, 512, f2, 1e-3, 15, '$(1+0.1\cos(t))\cos(10t)$', 'Fig2.png')
	calc_fft(-4*np.pi, 4*np.pi, 512, f3, 1e-3, 15, '$sin^3(t)$', 'Fig3.png')
	calc_fft(-4*np.pi, 4*np.pi, 512, f4, 1e-3, 15, '$cos^3(t)$',  'Fig4.png')
	calc_fft(-4*np.pi, 4*np.pi, 512, f5, 1e-3, 30, '$cos(20t + 5cos(t))$', 'Fig5.png')
	
	err, N, T = estim_fft()
	print(err, N, T/np.pi)
	calc_fft(-T/2, T/2, N, gs, 1e-3, 5, 'Gaussian - Calculated', 'Fig6.png')
	calc_fft(-T/2, T/2, N, gs_true, 1e-3, 5, 'Gaussian - True Plot', 'Fig7.png', fft=False)

main()
