"""
Title			: dft.py
Author			: Harish R EE20B044
Created			: Apr 15 2022
Last Modified	: May 12 2022
Purpose			: Understand Digital Fourier Transforms in python
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

def f1(t):
	return np.sin(5*t)

def f2(t):
	return (1+0.1*np.cos(t))*np.cos(10*t)

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

def calc_fft(xn, xp, N, f, thr, xlim, title, name, fft=True):
	t = np.linspace(xn, xp, N+1)[:-1]
	w = (N/(xp-xn))*np.linspace(-np.pi, np.pi, N+1)[:-1]

	if fft:
		Y = nf.fftshift(nf.fft(nf.ifftshift(f(t))))/float(N)
	else:
		Y = f(w)
	fig, ax = plt.subplots(2,1)
	ax[0].plot(w, abs(Y), lw=2)
	ax[0].set_xlim([-xlim, xlim])
	ax[0].set_ylabel(r'$|Y|$')
	ax[0].set_title(r'Spectrum of ' + title)
	ax[0].grid(True)
	if thr:
		ii = np.where(np.abs(Y)>thr)
	else:
		ii = np.where(Y!=Y)
   	
	phi = np.angle(Y)
	phi[np.where(np.abs(Y)<1e-6)]=0

	ax[1].plot(w, phi, 'ro', lw=2)
	ax[1].plot(w[ii], phi[ii], 'go', lw=2)
	ax[1].set_xlim([-xlim, xlim])
	ax[1].set_ylabel(r'Phase of $Y$')
	ax[1].set_xlabel(r'$\omega$')
	ax[1].grid(True)
	
	fig.savefig(PATH+name)

def estim_fft(fn = gs, true_fn = gs_true, tol=1e-6, N=128):

	T = 2*np.pi
	Y_o = 0
	err = tol+10
	iters = 0
		
	while err>tol:  
		x = np.linspace(-T/2,T/2,N+1)[:-1]
		w = np.linspace(-N*np.pi/T,N*np.pi/T,N+1)[:-1]
		y = fn(x)
		Y = nf.fftshift(nf.fft(nf.ifftshift(fn(x))))*T/(2*np.pi*N)
		err = sum(abs(Y[::2]-Y_o))
		Y_o = Y
		iters+=1
		T*=2
		N*=2

	return sum(abs(Y-true_fn(w))), N, T

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
