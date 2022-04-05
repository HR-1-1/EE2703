"""
Title			: laplaceTransform.py
Author			: Harish R EE20B044
Last Modified 	: Mar 15 2022
Purpose 		: Assert laplace dominance and leverage scipy.signal
Inputs			: Lite bro
"""

import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt

parameters = {'axes.labelsize': 12, 
			  'axes.titlesize': 15, 
			  'legend.fontsize': 10, 
		 	 'mathtext.fontset':'cm'} 

PATH='./plots/'

def spring(a, b, time):
	"""
	Inputs : 
	a -> driver frequency
	b -> decay factor
	"""
	num = np.poly1d([1, b])
	den = np.polymul([1,0,2.25],[1,2*b,a**2+b**2])
	H = sp.lti(num, den)
	t, x = sp.impulse(H, None, np.linspace(0,time,1000))
	return t,x

def spring_lsim(fi, n, time, ax):

	t = np.linspace(0, time, 1000)
	num = np.poly1d([1])
	den = np.poly1d([1,0,2.25])
	H = sp.lti(num, den)
	for i in range(n):
		f = round(fi+i*0.05,2)
		u = np.cos(f*t)*np.exp(-0.05*t)
		tout, yout, xout = sp.lsim(H, U=u, T=t)		
		ax.plot(tout, yout, linewidth=1.5, label=r'$Driver freq ='+str(f)+'$')
	return ax

def coupled_spring(time):
	
	H1 = sp.lti(np.poly1d([1,0,2]), np.poly1d([1,0,3,0]))
	H2 = sp.lti(np.poly1d([2]), np.poly1d([1,0,3,0]))
	t, x = sp.impulse(H1, None, np.linspace(0, time, 1000))
	t, y = sp.impulse(H2, None, np.linspace(0, time, 1000))
	return t, x, y

def two_port(R, C, L, ax):
	
	num = np.poly1d([1])
	den = np.poly1d([L*C,C*R,1])
	H = sp.lti(num, den)
	w,S,phi = H.bode()
	ax[0].semilogx(w,S)
	ax[1].semilogx(w,phi)
	return H, ax

def main():
		
	plt.rcParams.update(parameters) 
	t1, x1 = spring(1.5, 0.5 , 50)
	fig1, ax1 = plt.subplots(figsize=(8,8))
	ax1.plot(t1,x1, label='x')
	ax1.set_title(r'Time response of a spring forced by $f(t) = \cos(1.5t)\exp(-0.5t)u_0(t)$')
	ax1.set_xlabel(r'Time$\longrightarrow$')
	ax1.set_ylabel(r'Position$\longrightarrow$')
	ax1.legend(loc='best', shadow=True, framealpha=1)
	ax1.grid()
	fig1.savefig(PATH+'Figure1.png')


	t2, x2 = spring(1.5, 0.05, 50)
	fig2, ax2 = plt.subplots(figsize=(8,8))
	ax2.plot(t2,x2, label='x')
	ax2.set_title(r'Time response of a spring forced by $f(t) = \cos(1.5t)\exp(-0.05t)u_0(t)$')
	ax2.set_xlabel(r'Time$\longrightarrow$')
	ax2.set_ylabel(r'Position$\longrightarrow$')
	ax2.legend(loc='best', shadow=True, framealpha=1)
	ax2.grid()
	fig2.savefig(PATH+'Figure2.png')
	
	
	fig3, ax3 = plt.subplots(figsize=(8,8))
	ax3 = spring_lsim(1.4, 5, 50, ax3)
	ax3.set_title(r'Time response of a spring forced by $f(t) = \cos(f*t)\exp(-0.05t)u_0(t)$')
	ax3.set_xlabel(r'Time$\longrightarrow$')
	ax3.set_ylabel(r'Position$\longrightarrow$')
	ax3.legend(loc='best', shadow=True, framealpha=1)
	ax3.grid()
	fig3.savefig(PATH+'Figure3.png')
	
	fig4, ax4 = plt.subplots(figsize=(8,8))
	t, x, y = coupled_spring(20)
	ax4.plot(t,x,label='x')
	ax4.plot(t,y,label='y')
	ax4.set_title(r'Time evolution of a coupled spring system')
	ax4.set_xlabel(r'Time$\longrightarrow$')
	ax4.set_ylabel(r'Position$\longrightarrow$')
	ax4.legend(loc='best', shadow=True, framealpha=1)
	ax4.grid()
	fig4.savefig(PATH+'Figure4.png')

	fig5, ax5 = plt.subplots(2,1,figsize=(8,8))
	H, ax5 = two_port(100, 1e-6, 1e-6, ax5)
	ax5[0].set_title(r'Magnitude response of a second order LCR system')
	ax5[0].set_xlabel(r'Time$\longrightarrow$')
	ax5[0].set_ylabel(r'Magnitude$\longrightarrow$')
	ax5[0].grid()
	ax5[1].set_title(r'Phase response of a second order LCR system')
	ax5[1].set_xlabel(r'Time$\longrightarrow$')
	ax5[1].set_ylabel(r'Phase$\longrightarrow$')
	ax5[1].grid()
	fig5.savefig(PATH+'Figure5.png')
	
	fig6, ax6 = plt.subplots(2,1,figsize=(8,8))
	t = np.linspace(0, 0.1, 1000000)
	u = np.cos(1e3*t)-np.cos(1e6*t)
	tout, yout, xout = sp.lsim(H, U=u, T=t)
	ax6[0].plot(tout, yout,'r')
	ax6[0].set_title(r'Output Voltage of a Second order system excited with $V_i(t)=\cos(10^3*t)-\cos(10^6*t)$')
	ax6[0].set_xlabel(r'Time$\longrightarrow$')
	ax6[0].set_ylabel(r'Voltage$\longrightarrow$')
	ax6[0].grid()
	ax6[1].plot(tout, yout,'r')
	ax6[1].set_title(r'Zoomed out view for time less than 30 $\mu$ seconds')
	ax6[1].set_xlabel(r'Time$\longrightarrow$')
	ax6[1].set_ylabel(r'Voltage$\longrightarrow$')
	ax6[1].grid()
	ax6[1].set_ylim(0, 0.35)
	ax6[1].set_xlim(0, 3e-5)
	fig6.savefig(PATH+'Figure6.png')
	
	plt.show()
main()
