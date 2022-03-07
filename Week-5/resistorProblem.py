"""
Title			: resistorProblem.py
Author			: Harish R EE20B044
Last Modified 	: Mar 5 2022
Purpose 		: Numerical solution to current flow distribution in a resistor
Inputs			: 2-D dimensions of copper plate, radius of central lead, 
					Number of iterations
"""

"""
Some comments on the plate : 
1 cm X 1 cm in size
X goes from -0.5 to 0.5
Y goes from -0.5 to 0.5
"""

import argparse
from pylab import *
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.pyplot as plt

def debug(name, value):
	global args
	if args.debug:
		print("The value of " + name)
		print(value)
		print('*'*5)

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
	xx, yy = meshgrid(np.linspace(-0.5, 0.5, args.length), np.flip(np.linspace(-0.5, 0.5, args.width)))
	ii = np.where(xx*xx + yy*yy <= 1.05*(args.radius**2/(args.length*args.width)))
	phi[ii] = 1.0
	
	debug("xx", xx)
	debug("yy", yy)
	debug("ii", ii)
	debug("phi", phi)

	fig, ax = plt.subplots(figsize =(8,8))
	cont = ax.contour(xx, yy, phi, 20)
	ax.scatter(ii[0], yy[ii[1]][0], 'ro')
	ax.set_title('Contour Plot of potential')
	ax.set_xlabel('x axis')
	ax.set_ylabel('y axis')

	plt.show()
	


main()
