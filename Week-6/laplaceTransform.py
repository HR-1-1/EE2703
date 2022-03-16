"""
Title			: laplaceTransform.py
Author			: Harish R EE20B044
Last Modified 	: Mar 15 2022
Purpose 		:
Inputs			: 
"""
import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt
parameters = {'axes.labelsize': 12, 
			  'axes.titlesize': 15, 
			  'legend.fontsize': 10, 
		 	 'mathtext.fontset':'cm'} 

PATH='./plots/'

def ltiSolve(num, den, time):
	"""
	
	"""
	H = sp.lti(num, den)
	t, x = sp.impulse(H, None, np.linspace(0,time,1000))
	return t,x

def main():
		
	plt.rcParams.update(parameters) 
	t1, x1 = ltiSolve(np.poly1d([1,0.5]), np.polymul([1,0,2.25],[1,1,2.5]), 50)
	fig1, ax1 = plt.subplots(figsize=(8,8))
	ax1.plot(t1,x1)
	ax1.set_title(r'Time response of a spring forced by $f(t) = \cos(1.5t)\exp(-0.5t)u_0(t)$')
	ax1.set_xlabel(r'Time$\longrightarrow$')
	ax1.set_ylabel(r'Position$\longrightarrow$')
	fig1.savefig(PATH+'Figure1.png')

	t2, x2 = ltiSolve(np.poly1d([1,0.05]), np.polymul([1,0,2.25],[1,1,2.5]), 50)
	fig2, ax2 = plt.subplots(figsize=(8,8))
	ax2.plot(t2,x2)
	ax2.set_title(r'Time response of a spring forced by $f(t) = \cos(1.5t)\exp(-0.05t)u_0(t)$')
	ax2.set_xlabel(r'Time$\longrightarrow$')
	ax2.set_ylabel(r'Position$\longrightarrow$')
	fig2.savefig(PATH+'Figure2.png')

	plt.show()
	
main()
