##############################################################
#Title   : fourierFit.py
#Author  : Harish R EE20B044
#Date    : Feb 23 2022
#Purpose : 1. Linear fit a curve to its fourier series
#		   2. Study the difference in the best fit and integration approach
#Inputs  : Sit back and relax
################################################################

import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import clabel
from scipy import special as sp
from scipy.linalg import lstsq
from scipy.integrate import quad
PATH = "./plots/"

def exp(x):
	return np.exp(x)
def ccos(x):
	return np.cos(np.cos(x))

def cos_cf(x, k, f):
	return f(x)*np.cos(k*x)
def sin_cf(x, k, f):
	return f(x)*np.sin(k*x)

def calcFourierCfs (f, N):
	cf_a = np.zeros(int(N/2)+1)
	cf_b = np.zeros(int(N/2)+1)
	cf_a[0] = quad(cos_cf, 0, 2*np.pi, args=(0, f))[0]/(2*np.pi)
	for i in range(1, int(N/2)+1):
		cf_a[i] = quad(cos_cf, 0, 2*np.pi, args=(i, f))[0]/(np.pi)
		cf_b[i] = quad(sin_cf, 0, 2*np.pi, args=(i, f))[0]/(np.pi)
	cfs = np.zeros(N)
	cfs[0] = cf_a[0]
	cfs[1::2] = cf_b[1:]
	cfs[2::2] = cf_a[1:]
	return cfs

def plot_data(fn='line',x=None,y=None,z=None,
				xl=r'x-axis',yl=r'y-axis',title='title',fmt=None,
				legends=None,stdev=None,logtype=None,
				save=None):
		
	parameters = {'axes.labelsize': 12,
			  'axes.titlesize': 15,
			  'legend.fontsize': 10,
			  'mathtext.fontset':'cm'}
	
	plt.rcParams.update(parameters)
	fig, ax = plt.subplots(figsize =(8,8))
	
	if fn=='line':
		if fmt:
			ax.plot(x,y,fmt)
		else:
			ax.plot(x,y)
	elif fn=='err':
		ax.errorbar(x,y,stdev,fmt=fmt)
	elif fn=='contour':
		cont = ax.contour(x,y,z,20)
		clabel(cont, np.linspace(0.025,0.100,4), inline=True)
	elif fn=='dashed':
		ax.plot(x,y,linestyle='dashed',marker='o', markersize=5)
	elif fn=='scatter':
		ax.scatter(x,y,z)

	if logtype=='x' or logtype=='xy':
		ax.set_xscale('log')
	
	if logtype == 'y' or logtype=='xy':
		ax.set_yscale('log')
	
	if legends:
		ax.legend(tuple(legends), loc='upper right')
	
	ax.set(xlabel=xl, ylabel=yl, title=title)
	ax.grid()
	
	if save:
		fig.savefig(PATH+save)
	return fig

def linfit(t,f):
	a_fit = np.zeros(f.shape[1]).astype(float)
	b_fit = np.zeros(f.shape[1]).astype(float)

	for i in range(f.shape[1]):
		a_fit[i], b_fit[i] = lstsq(np.c_[sp.jn(2,t), t], f[:,i])[0]
	
	return a_fit, b_fit

def main():
		
	x = np.linspace(-2*np.pi, 4*np.pi, num=600)
	exp_x = plot_data(fn='line', x=x, y=np.c_[exp(x), exp(x%(2*np.pi))], 
						xl=r'$x$', yl=r'True Function',
						title=r'$\exp(x)$',
						legends= ['True Function', 'Fourier Expansion'],
						logtype='y',
						save='exp_x.png')
	
	ccos_x = plot_data(fn='line', x=x, y=np.c_[ccos(x), ccos(x%(2*np.pi))],
						xl=r'$x$', yl=r'True Function',
						title=r'$\cos(cos(x))$',
						legends=['True Function', 'Fourier Expansion'],
						logtype='y',
						save='ccos_x.png')
	
	ccos_cfs = calcFourierCfs(ccos,51)
	exp_cfs = calcFourierCfs(exp,51)
	
	exp_fr_semi = plot_data(fn='scatter', x=np.arange(25), y=abs(exp_cfs[1::2]),	z=abs(exp_cfs[2::2]),
							xl=r'$i$', yl=r'',
							title=r'First 51 Fourier Coefficients - $e^x$ ',
							legends=['$a_n$', '$b_n$'],
							logtype='y',
							save='exp_fr_semi.png')
	
	ccos_fr_semi = plot_data(fn='scatter', x=np.arange(25), y=abs(exp_cfs[1::2]),	z=abs(exp_cfs[2::2]),
							xl=r'$i$', yl=r'',
							title=r'First 51 Fourier Coefficients - $e^x$ ',
							legends=['$a_n$', '$b_n$'],
							logtype='y',
							save='exp_fr_semi.png')
	
	exp_fr_log = plot_data(fn='scatter', x=np.arange(25), y=abs(exp_cfs[1::2]),	z=abs(exp_cfs[2::2]),
						xl=r'$i$', yl=r'',
						title=r'First 51 Fourier Coefficients - $cos(cos(x)$',
						legends=['$a_n$', '$b_n$'],
						logtype='xy',
						save='exp_fr_log.png')

	ccos_fr_log = plot_data(fn='scatter', x=np.arange(25), y=abs(exp_cfs[1::2]),	z=abs(exp_cfs[2::2]),
						xl=r'$i$', yl=r'',
						title=r'First 51 Fourier Coefficients - $cos(cos(x)$',
						legends=['$a_n$', '$b_n$'],
						logtype='xy',
						save='ccos_fr_log.png')

	
main()
