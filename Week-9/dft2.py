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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

parameters = {'axes.labelsize': 10, 
              'axes.titlesize': 12, 
              'legend.fontsize': 8,  
             'mathtext.fontset':'cm',
             'figure.figsize': (10,15),
			'figure.max_open_warning':0} 

PATH = './plots/'

def debug(var):
	for name, value in var.items():
		print("The value of " + name)
		print(value)
		print('*'*5)

def f1(t,w=0.86):
	return np.cos(w*t)**3

def f2(t,w=1.5, d=0.5):
	return np.cos(w*t+d)

def f3(t, w=1.5, d=0.5):
	return np.cos(w*t+d) + 0.1*np.random.randn(len(t))

def f4(t):
	return np.cos(16*(1.5 + t/(2*np.pi))*t) 

def estim(w, Y, s=1e-4, wd=1):
	"""
	Function to estimate time-signal parameters 
	given the frequency spectra
	"""
	ii = np.where(w>0)
	w0 = (sum(abs(Y[ii])**2*w[ii])/sum(abs(Y[ii])**2)) #Weighted average of w
    
	ii = np.sort(np.where(np.logical_and(np.abs(Y)>s, w>0))[0])
	pts = ii[1:wd+1] #Extracting wd points that satisy the criteria
     
	dlt = np.sum(np.angle(Y[pts]))/len(pts) #Sum of window points

	return w0,dlt

def calc_fft(xl, N, f, xlim, title, name, window=False, t0=None, tlim = False):
	"""
	A paramerized function to calculate the fft of 
	Aperiodic signals by performing windowing
	and hamming
	"""
	if(tlim):
		t = t0
	else:
		t = np.linspace(-xl,xl,N+1)[:-1]
    
	fmax = 1/(t[1]-t[0])
	y = f(t)
    
	if (window):
		y=y*nf.fftshift(0.54+0.46*np.cos(2*np.pi*np.arange(N)/N))
    
	y[0]=0 
	Y = nf.fftshift(nf.fft(nf.fftshift(y)))/float(N)
	w = np.linspace(-np.pi*fmax, np.pi*fmax, N+1)[:-1]
	
	fig, ax = plt.subplots(2,1)
	ax[0].plot(w, abs(Y), lw=2)
	ax[0].set_xlim([-xlim, xlim])
	ax[0].set_ylabel(r'$|Y| \longrightarrow$')
	ax[0].set_title(r'Spectrum of ' + title)
	ax[0].grid(True)
	
	phi = np.angle(Y)
	phi[np.where(abs(Y)<3e-3)] = 0
	
	ax[1].plot(w, phi, 'ro', lw=2)
	ax[1].set_xlim([-xlim, xlim])
	ax[1].set_ylabel(r'Phase of $Y \longrightarrow$')
	ax[1].set_xlabel(r'$\omega \longrightarrow$')
	ax[1].grid(True)
	
	if(name):
		fig.savefig(PATH+name)
    
	return w,Y

def freq_time(window=False,e=7):
	"""
	Get a frequency-time plot for chirp function
	"""
	t=np.linspace(-np.pi,np.pi,1025)[:-1]
	t_arrays=np.split(t,16)

	Y_mags=np.zeros((16,64))
	Y_angles=np.zeros((16,64))
	
	#splitting array and doing fft
	for i in range(len(t_arrays)):
		w,Y = calc_fft(10, 64, f4, 60, r'Spectrum of chirp function', t0=t_arrays[i], name=None, tlim=True, window=window)
		Y_mags[i] =  abs(Y)
		Y_angles[i] = np.angle(Y)
	
	fig1 = plt.figure()
	ax1 = fig1.add_subplot(111, projection='3d')

	fmax = 1/(t[1]-t[0])
	t=t[::64]
	w = np.linspace(-fmax*np.pi,fmax*np.pi,64+1);w=w[:-1]
	t,w = np.meshgrid(t,w)

	surf=ax1.plot_surface(w,t,Y_mags.T,cmap=cm.coolwarm,linewidth=0, antialiased=False)
	fig1.colorbar(surf, shrink=0.5, aspect=5)
	ax1.set_ylabel(r'Frequency$\longrightarrow$')
	ax1.set_xlabel(r'Time$\longrightarrow$')

	fig1.savefig(PATH+"Figure"+str(e)+".png")

	fig2 = plt.figure() 
	ax2 = fig2.add_subplot(111, projection='3d')
	
	surf=ax2.plot_surface(w,t,Y_angles.T,cmap=cm.coolwarm,linewidth=0, antialiased=False)
	fig2.colorbar(surf, shrink=0.5, aspect=5)
	ax2.set_ylabel(r'Frequency$\longrightarrow$')
	ax2.set_xlabel(r'Time$\longrightarrow$')

	fig2.savefig(PATH+"Figure"+str(e+1)+".png")

def main():
	
	#plt.rcParams.update(parameters) 
	calc_fft(4*np.pi, 64*4, f1, xlim= 3, title= r'$cos^3(w_0t)$',name= 'Figure1.png')
	calc_fft(4*np.pi, 64*4, f1, xlim= 3, title= r'$cos^3(w_0t)$',name= 'Figure2.png', window=True)
	
	w,Y = calc_fft(np.pi, 128, f2, xlim= 3, title = r'$cos(w_0t + \delta)$', name = 'Figure3.png', window=True)
	w0,delta = estim(w,Y) #Estimate the time-signal parameters
	debug({"Estimated w0":w0, "Estimated delta":delta})
	w,Y = calc_fft(np.pi, 128, f3, xlim= 3, title = r'$cos(w_0t + \delta) + Noise$', name = 'Figure4.png', window=True)
	w0,delta = estim(w,Y)
	debug({"Estimated w0":w0, "Estimated delta":delta})
	calc_fft(np.pi, 1024, f4, xlim= 60, title= r'Chirp function', name='Figure5.png')
	calc_fft(np.pi, 1024, f4, xlim= 60, title= r'Chirp function', name='Figure6.png', window=True)
	
	freq_time(True) #Calculate with Hamming function
	freq_time(False,9) #Calculate without Hamming function
main()
