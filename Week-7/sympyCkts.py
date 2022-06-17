"""
Title			: sympyCkts.py
Author			: Harish R EE20B044
Created			: Mar 28 2022
Last Modified 	: Apr 6 2022
Purpose 		: Analyze circuits in python using sympy module
Inputs			: Peace Potiya
"""

from __future__ import division
import sympy as sy
import scipy.signal as sp
import pylab as p
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display

sy.init_printing()  # LaTeX like pretty printing for IPython
parameters = {'axes.labelsize': 10, 
			  'axes.titlesize': 12, 
			  'legend.fontsize': 8, 
		 	 'mathtext.fontset':'cm',
			 'figure.figsize': (10,15)} 

PATH='./plots/'

def sympy_to_lti(xpr, s=sy.Symbol('s'), display=False):
	""" Convert Sympy transfer function polynomial to Scipy LTI """
	if display:
		print("Sympy transfer function ", xpr)
	num, den = sy.simplify(xpr).as_numer_denom()  # expressions
	l_num = np.array(sy.Poly(num, s).all_coeffs(), dtype=float)
	l_den = np.array(sy.Poly(den, s).all_coeffs(), dtype=float)
	H_lti = sp.lti(l_num, l_den)
	if display:
		print("LTI Transfer function", H_lti)
	return H_lti

def lowpass(R1, R2, C1, C2, G, Vi, display=False):
	s = sy.Symbol('s')
	A = sy.Matrix([[0, 0, 1, -1/G], 
				[-1/(1+s*R2*C2), 1, 0, 0], 
				[0, -G, G, 1], 
				[-1/R1-1/R2-s*C1, 1/R2, 0, s*C1]])
	b = sy.Matrix([0, 0, 0, -1*Vi/R1]) #Doubt whether -1 needs to be put --> Yes, Required
	V = A.inv()*b
	Vo = V[3]
	if display:
		print("V0 :", Vo)
	hf = sy.lambdify(s, Vo, 'numpy')
	return s, Vo, hf

def highpass(R1, R3, C1, C2, G, Vi, display=False): #Are the Op-Amp signs right? No
	s = sy.Symbol('s')
	A = sy.Matrix([[0, -1, 0, 1/G], 
				[s*R3*C2/(1+s*R3*C2), 0, -1, 0], 
				[0, G, -G, 1], 
				[-1/R1-s*C2-s*C1, 0, s*C2, 1/R1]])
	b = sy.Matrix([0, 0, 0, -Vi*s*C1])
	V = A.inv()*b
	Vo = V[3]
	if display:
		print("V0 :", Vo)
	hf = sy.lambdify(s, Vo, 'numpy')
	return s, Vo, hf

def main():
		
	plt.rcParams.update(parameters) 
	s, Vo, hf = lowpass(1e5, 1e5, 1e-9, 1e-9, 1.586, 1, display=False)
	lti_H1 = sympy_to_lti(Vo, s)
	w, mag, phase = sp.bode(lti_H1, w=np.linspace(1, 1e6, 1000000))
	t = np.linspace(0, 1e-2, 1000000)
	fig1, ax1 = plt.subplots(2,1)
	ax1[0].semilogx(w, mag)
	ax1[0].set_title(r'Magnitude response of a Butterworth filter')
	ax1[0].set_xlabel(r'Frequency$\longrightarrow$')
	ax1[0].set_ylabel(r'$20log(|(H(j*\omega)|)\longrightarrow$')
	ax1[0].grid()
	
	ax1[1].semilogx(w, phase)
	ax1[1].set_title(r'Phase response of a Butterworth filter')
	ax1[1].set_xlabel(r'Frequency$\longrightarrow$')
	ax1[1].set_ylabel(r'$\angle (H(j\omega)|\longrightarrow$')
	ax1[1].grid()

	fig1.savefig(PATH+'Figure1.png')
	
	u = np.sin(2e3*np.pi*t) + np.cos(2e6*np.pi*t)
	tout, yout, xout = sp.lsim(lti_H1, u, t)    
	fig2, ax2 = plt.subplots()
	ax2.plot(t, u, label=r'Input Sum of sinusoids')
	ax2.plot(tout, yout, linewidth=1.5, label=r'Output Voltage')
	ax2.set_title(r'Output Voltage of a Butterworth filter driven by $V_i(t)=\sin(2000*\pi*t)+\cos(2*10^6*\pi*t)*u_0(t)$')
	ax2.set_xlabel(r'Time$\longrightarrow$')
	ax2.set_xlim(0, 1e-2)	
	ax2.set_ylabel(r'$Voltage (Volts)\longrightarrow$')
	ax2.legend(loc='best', shadow=True, framealpha=1)
	ax2.grid()
	fig2.savefig(PATH+'Figure2.png')

	s, Vo, hf = highpass(1e4, 1e4, 1e-9, 1e-9, 1.586, 1, display=True)
	lti_H2 = sympy_to_lti(Vo, s)
	w, mag, phase = sp.bode(lti_H2, w=np.linspace(1, 1e6, int(1e6)))

	fig3, ax3 = plt.subplots(2,1)
	ax3[0].semilogx(w, mag)
	ax3[0].set_title(r'Magnitude response of a Butterworth HighPass filter')
	ax3[0].set_xlabel(r'Frequency$\longrightarrow$')
	ax3[0].set_ylabel(r'$20log(|(H(j*\omega)|)\longrightarrow$')
	ax3[0].grid()
	
	ax3[1].semilogx(w, phase)
	ax3[1].set_title(r'Phase response of a Butterworth HighPass filter')
	ax3[1].set_xlabel(r'Frequency$\longrightarrow$')
	ax3[1].set_ylabel(r'$\angle (H(j\omega)|\longrightarrow$')
	ax3[1].grid()

	fig3.savefig(PATH+'Figure3.png')
	
	f = 1e2*np.pi
	b = 1000
	u_low = np.heaviside(t, 1)*(np.exp(-b*t)*np.sin(f*t))
	tout, yout, xout = sp.lsim(lti_H2, U=u_low, T=t)    
	fig4, ax4 = plt.subplots()
	ax4.plot(t, u_low, label=r'Sum of Low freq Dampened Sinusoids')
	ax4.plot(t, yout, linewidth=1.5, label=r'Output Voltage')
	ax4.set_title(r'Output Voltage of a BHPF driven by $V_i(t)=\sin(2*10^3*\pi*t)+\cos(2*10^6*\pi*t)*u_0(t)$')
	ax4.set_xlabel(r'Time$\longrightarrow$')
	ax4.set_ylabel(r'$Voltage (Volts)\longrightarrow$')
	ax4.set_xlim(0, 1e-3)
	ax4.legend(loc='best', shadow=True, framealpha=1)
	ax4.grid()
	fig4.savefig(PATH+'Figure4.png')
	
	f = 1e6*np.pi
	b = 1000
	u_high = np.heaviside(t, 1)*(np.exp(-b*t)*np.sin(f*t))
	tout, yout, xout = sp.lsim(lti_H2, U=u_high, T=t)    
	fig5, ax5 = plt.subplots()
	ax5.plot(t, u_high, label=r'Sum of High freq Dampened Sinusoids')
	ax5.plot(t, yout, linewidth=1.5, label=r'Output Voltage')
	ax5.set_title(r'Output Voltage of a BHPF driven by $V_i(t)=\sin(2*10^3*\pi*t)+\cos(2*10^6*\pi*t)*u_0(t)$')
	ax5.set_xlabel(r'Time$\longrightarrow$')
	ax5.set_ylabel(r'$Voltage (Volts)\longrightarrow$')
	ax5.legend(loc='best', shadow=True, framealpha=1)
	ax5.set_xlim(0, 1e-3)
	ax5.grid()
	fig5.savefig(PATH+'Figure5.png')

# Unit functions Responses

	tout, yout  = sp.step(lti_H2, None, T=t)    
	fig6, ax6 = plt.subplots()
	ax6.plot(t, yout, linewidth=1.5, label=r'Output Voltage')
	ax6.set_title(r'Response of a BHPF driven to Unit step function')
	ax6.set_xlabel(r'Time$\longrightarrow$')
	ax6.set_ylabel(r'$Voltage (Volts)\longrightarrow$')
	ax6.legend(loc='best', shadow=True, framealpha=1)
	ax6.set_xlim(0, 1e-3)
	ax6.grid()
	fig6.savefig(PATH+'Figure6.png')

	tout, yout  = sp.step(lti_H1, None, T=t)    
	fig7, ax7 = plt.subplots()
	ax7.plot(t, yout, linewidth=1.5, label=r'Output Voltage')
	ax7.set_title(r'Response of a BLPF driven to Unit step function')
	ax7.set_xlabel(r'Time$\longrightarrow$')
	ax7.set_ylabel(r'$Voltage (Volts)\longrightarrow$')
	ax7.legend(loc='best', shadow=True, framealpha=1)
	ax7.set_xlim(0, 1e-3)
	ax7.grid()
	fig7.savefig(PATH+'Figure6.png')
	plt.show()

main()
