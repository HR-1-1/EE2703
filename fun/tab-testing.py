"""
Title			: tab-testing.py
Author			: Harish R EE20B044
Created 		: Apr 15 2022
Last Modified 	: 
Purpose 		: Numerical solution to current flow distribution in a resistor
Inputs			: Number of iterations
"""

"""
Some comments on the plate : 
45 mm X 19.35 mm in size
X goes from -22.5 to 22.5
Y goes from -9.675 to 9.675
"""

import argparse
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.pyplot as plt

PLOTS = './plots/'

def debug(var):
	global args
	if args.debug:
		for name, value in var.items():
			print("The value of " + name)
			print(value)
			print('*'*5)

def laplaceSolver(phi, Niter):
	error = []
	for n in range(Niter):
		oldphi = phi.copy()

		# Update Potential
		phi[1:-1, 1:-1] = 0.25*(phi[1:-1, 0:-2]+phi[1:-1, 2:]+phi[0:-2, 1:-1]+phi[2:, 1:-1])

		# Boundary Conditions
		phi[1:-1, 0] = phi[1:-1, 1]
		phi[1:-1, -1] = phi[1:-1, -2]
		#phi[0, :] = phi[1, :]
		
		phi[:,0] = 0
		phi[:,-1] = 1.0

		error.append(np.abs(phi-oldphi).max())
	
	return np.array(error)

def cErr(N, A, B):
	return -(A/B)*np.exp(B*(N+0.5))

def getErrorTol(bestFit, Niter, error_tol):
	cErr = []
	for n in range(0, Niter):
		cErr.append(cErr(n, np.exp(bestFit[0]), bestFit[1]))
		if(cErr[n-1] <= error_tol):
			return cErr[n-1], n
	return cErr[-1], Niter

def main():
	parser = argparse.ArgumentParser(description="Program to calculate current flow distribution in a resistor", epilog="Feel free to checkout the source code")
	parser.add_argument('-l', '--length', nargs='?', type=int, help="Number of points for length of copper plate", default=25)
	parser.add_argument('-w', '--width', nargs='?', type=int, help="Number of points for width of copper plate", default=25)
	parser.add_argument('-r', '--radius', nargs='?',  type=int, help="radius of the central lead", default=8)
	parser.add_argument('-i', '--iter', nargs='?',  type=int, help="Number of iterations to perform", default=1500)
	parser.add_argument('-d', '--debug', help="Better debugging in code", action="store_true")
		
	global args 
	args = parser.parse_args()
	
	phi = np.zeros((args.width, args.length))
	xunit = np.linspace(-22.5, 22.5, args.length)
	yunit = np.flip(np.linspace(-9.675, 9.675, args.width))
	xx, yy = meshgrid(xunit, yunit)
	phi[:,-1] = 1.0
	phi[:,0] = 0
	
	debug({"xx":xx,"yy":yy,"phi":phi})

	fig1, ax = plt.subplots(figsize =(8,8))
	cont = ax.contour(xx, yy, phi, 20)
	#ax.scatter(xunit[ii[0]], yunit[ii[1]], c='r', label='1 V')
	ax.set_title('Contour Plot of potential function $\phi$')
	ax.set_xlabel(r'X axis $\longrightarrow$')
	ax.set_ylabel(r'Y axis $\longrightarrow$')
	ax.legend()
	ax.grid(True)
	fig1.savefig(PLOTS + "Figure1.png")
	
	error = laplaceSolver(phi, args.iter)

	fig2, ax = plt.subplots(figsize =(8,8))
	ax.semilogy(np.arange(args.iter), error)
	ax.set_title('Semilog plot of error')
	ax.set_xlabel(r'Number of iteration$\longrightarrow$')
	ax.set_ylabel(r'Error$\longrightarrow$')
	ax.grid(True)
	fig2.savefig(PLOTS + "Figure2.png")
	
#	Niter = args.iter
#
#	bestFit = np.polyfit(np.arange(Niter), np.log(error), 1)
#	bestFit500 = np.polyfit(np.arange(500, Niter), np.log(error[500:]), 1)
#	
#	bestFitValues = np.arange(Niter)*bestFit[0] + bestFit[1]
#	bestFit500Values = np.arange(500, Niter)*bestFit500[0] + bestFit500[1]
#	
#	fig3, ax = plt.subplots(figsize =(8,8))
#	ax.semilogy(np.arange(Niter)[::50], np.exp(bestFitValues[::50]), 'ro', label='Entire curve Linear Fit')
#	ax.semilogy(np.arange(500, Niter)[::50], np.exp(bestFit500Values[::50]), 'gx', label='Linear Fit after 500')
#	ax.semilogy(np.arange(Niter), error, 'b', label='Original Error')
#	ax.set_title('Linear Fit Estimate of Error')
#	ax.set_xlabel(r'Number of iteration$\longrightarrow$')
#	ax.set_ylabel(r'Error Estimates$\longrightarrow$')
#	ax.legend()
#	ax.grid(True)
#	fig3.savefig(PLOTS + "Figure3.png")
#
#
	#cErr, end = getErrorTol(bestFit, Niter, 1e-5)
	#print(end, cErr)

	fig4 = plt.figure(4)
	ax = p3.Axes3D(fig4)
	ax.set_title('The 3-D surface plot of $\phi$')
	surfacePlot = ax.plot_surface(xx, yy, phi, rstride=1, cstride=1, cmap=cm.jet)
	cax = fig4.add_axes([1, 0, 0.1, 1])
	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$')
	ax.set_zlabel('$z$')
	fig4.colorbar(surfacePlot, cax=cax, orientation='vertical')
	fig4.savefig(PLOTS +'Figure4.png')


	fig5, ax = plt.subplots(figsize=(8,8))
	ax.contourf(xx, yy, phi, cmap=cm.jet)
	ax.set_title('Updated Contour Plot of potential function $\phi$')
	ax.set_xlabel(r'X axis $\longrightarrow$')
	ax.set_ylabel(r'Y axis $\longrightarrow$')
	ax.grid(True)
	fig5.savefig(PLOTS + "Figure5.png")

	Jx = np.zeros_like(phi)
	Jy = np.zeros_like(phi)

	Jx[1:-1, 1:-1] = (phi[1:-1, 0:-2] - phi[1:-1, 2:])*0.5
	Jy[1:-1, 1:-1] = (phi[2:, 1:-1] - phi[0:-2, 1:-1])*0.5

	fig6, ax = plt.subplots(figsize=(8,8))
	#ax.scatter(xunit[ii[0]], yunit[ii[1]], color='r', s=10, label='$V = 1V$ region')
	ax.quiver(xx, yy, Jx, Jy)
	ax.set_title('Vector Plot of current')
	ax.set_xlabel(r'X axis $\longrightarrow$')
	ax.set_ylabel(r'Y axis $\longrightarrow$')
	ax.legend()
	ax.grid(True)
	fig6.savefig(PLOTS + "Figure6.png")
	
	plt.show()
main()
