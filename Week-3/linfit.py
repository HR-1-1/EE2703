##############################################################
#Title   : linfit.py
#Author  : Harish R EE20B044
#Date    : Feb 14 2022
#Purpose : 1. Linear Fit a Bessel function to the given data
#		   2. 
#Inputs  : Textfile with a 10x10 matrix first column being time [Given as a CommandLine input]
################################################################

import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import clabel
from scipy import special as sp
from scipy.linalg import lstsq

A_ = 1.05
B_ = -0.105 
PATH = "./plots/"
N = 100

def g(t,a,b):
	return a*sp.jn(2,t)+b*t

def generate_data():
	
	# script to generate data files for the least squares assignment
	
	N=101                           # no of data points
	k=9                             # no of sets of data with varying noise

	# generate the data points and add noise
	t= np.linspace(0,10,N)              # t vector
	y=1.05*sp.jn(2,t)-0.105*t       # f(t) vector
	Y=np.meshgrid(y,np.ones(k),indexing='ij')[0] # make k copies
	scl=np.logspace(-1,-3,k)           # noise stdev
	n=np.dot(np.random.randn(N,k),np.diag(scl))     # generate k vectors
	yy=Y+n                          # add noise to signal

	return t,yy
	

def plot_data(fn,x,y,xl,yl,title,legends=None,stdev=None,z=None,loglog=False):
		
	parameters = {'axes.labelsize': 12,
			  'axes.titlesize': 15,
			  'legend.fontsize': 10,
			  'mathtext.fontset':'cm'}
	
	plt.rcParams.update(parameters)
	fig, ax = plt.subplots(figsize =(8,8))
	if fn=='line':
		ax.plot(x,y)
	elif fn=='err':
		ax.errorbar(x,y,stdev,fmt='ro')
	elif fn=='contour':
		cont = ax.contour(x,y,z,20)
		clabel(cont, np.linspace(0.025,0.100,4), inline=True)
	elif fn=='dashed':
		ax.plot(x,y,linestyle='dashed',marker='o', markersize=5)
	
	if loglog:
		ax.set_xscale('log')
		ax.set_yscale('log')

	if legends:
		ax.legend(tuple(legends), loc='upper right')
	
	ax.set(xlabel=xl, ylabel=yl, title=title)
	ax.grid()

	return fig

def linfit(t,f):
	a_fit = np.zeros(f.shape[1]).astype(float)
	b_fit = np.zeros(f.shape[1]).astype(float)

	for i in range(f.shape[1]):
		a_fit[i], b_fit[i] = lstsq(np.c_[sp.jn(2,t), t], f[:,i])[0]
	
	return a_fit, b_fit

def main(f):
	
	x,y = np.split(np.loadtxt(f), [1,], axis=1)
	original = g(x, A_, B_)
	std = np.logspace(-1,-3,9)
	legends = [r'$\sigma_'+str(i+1)+' = $' + str(round(x,3)) for i,x in enumerate(std)]
	legends.append("True Value")
	in_data = plot_data('line', x, np.c_[y,original],
						r'$t	\longrightarrow$',
						r'$F(t) = f(t)+noise	\longrightarrow$',
						r'Q4 : Data to be fitted in theory',
						legends)
	
	in_data.axes[0].plot(x,original, color='skyblue', linewidth=5)
	in_data.savefig(PATH+"noisy_data.png")
	
	legends = [r'$f(t)$', r'$ErrorBar$']
	err_bar = plot_data('err', x[::5],y[::5, 0],
						r't	$\longrightarrow$',
						r'',
						r'Q5 : Data points for $\sigma = '+str(round(std[0],3))+'$along with exact function',
						legends,
						std[0])
	
	err_bar.axes[0].plot(x, original, color='skyblue', linewidth=5)
	err_bar.savefig(PATH+"error_bar.png")
	
	M = np.c_[sp.jn(2,x), x]
	p = np.array([A_, B_])
	
	if np.allclose(original.squeeze(), np.matmul(M,p), 1e-2):
		print("The arrays are close")
	else :
		print("The arrays aren't close enough")
	
	A = np.linspace(0,2,21).squeeze()
	B = np.linspace(-0.2,0,21).squeeze()
	mse = np.zeros((A.shape[0], B.shape[0]))
	for i,a in enumerate(A):
		for j,b in enumerate(B):
			mse[i][j] = np.square(np.subtract(y[:,0],g(x,a,b).squeeze())).mean()
			
	
	contour = plot_data('contour', A, B,
						r'A    $\longrightarrow$',
						r'B	   $\longrightarrow$',
						r'Q8: Contour plot of $\epsilon_{ij}$',
						z=mse)

	min_A, min_B = np.where(mse == mse.min())
	print(A[min_A], B[min_B])
	contour.axes[0].plot(A[min_A], B[min_B], 'ro')
	contour.axes[0].annotate("Minima Point", (A[min_A],B[min_B]))
	contour.axes[0].plot(A_,B_,'ro')
	contour.axes[0].annotate("Exact Location", (A_, B_))
	contour.savefig(PATH+"contour_plot.png")
	
	a_err = []
	b_err = []

	for i in range(N):
		x,y = generate_data()
		a, b = linfit(x,y)
		a_err.append(np.square(a-A_))
		b_err.append(np.square(b-B_))
	
	a_err_mse = np.stack(a_err).mean(axis=0)
	b_err_mse = np.stack(b_err).mean(axis=0)

	legends = ['Aerr', 'Berr']
	err_plot = plot_data('dashed', std, np.c_[a_err_mse, b_err_mse], 
						r'Noise standard deviation	$\longrightarrow$',
						r'MS Error	$\longrightarrow$',
						r'Q10: Variation of error with noise',
						legends)

	err_plot.savefig(PATH+"error_plot.png")
	
	err_loglog = plot_data('dashed', std, np.c_[a_err_mse, b_err_mse], 
						r'Noise standard deviation	$\longrightarrow$',
						r'MS Error	$\longrightarrow$',
						r'Q11: Variation of error with noise',
						legends,
						loglog=True)
	
	err_loglog.savefig(PATH+"error_plot_loglog.png")



if(len(sys.argv)!=2):
	sys.exit("Correct Usage : python3 {0} <input-data.txt>".format(sys.argv[0]))

main(sys.argv[1])
