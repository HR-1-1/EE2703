"""
Title			: EE20B044.py
Author			: Harish R EE20B044
Created			: May 12 2022
Last Modified	: May 12 2022
Purpose			: Find antenna currents in a half-wave dipole antenna
Inputs			: Kindly look at python3 <file-name.py> --help 
"""


import argparse
from pylab import *
import matplotlib.pyplot as plt

PLOTS = './'

parameters = {'axes.labelsize': 10, 
              'axes.titlesize': 12, 
              'legend.fontsize': 8,  
             'mathtext.fontset':'cm',
             'figure.figsize': (10,15)} 

def debug(var):
    global args
    if args.debug:
        for name, value in var.items():
            print("The value of " + name)
            print(value.round(5))
            print('*'*5)

"""
	Comments on the Antenna :
		Antenna is an array with indices 0 to 2N (2N+1 elems)
		Point N is the feed of the array
"""

l=0.5 #Quarter Wavelength in meters
c=2.9979e8 #Speed of light in m/s
mu0=4e-7*np.pi #Permeability of free space

def ampMat(n,r):
	"""
	Ampere Law matrix H = M*J
	"""
	return np.identity(n,dtype=float)/(2*np.pi*r)

def mkDist(z,r):
	"""
	Distance Vectors from observer and source
	"""
	Z1, Z2 = meshgrid(z,z)
	return np.sqrt(np.square(Z1-Z2)+r**2)

def main():
	parser = argparse.ArgumentParser(description="Program to calculate antenna currents in a half-wave dipole antenna", epilog="Feel free to checkout the source code")
	parser.add_argument('-n', '--num', nargs='?', type=int, help="Number of section in each half section of antenna", default=4)
	parser.add_argument('-i', '--current', nargs='?', type=float, help="Current injected into the antenna", default=1.0)
	parser.add_argument('-r', '--radius', nargs='?',  type=float, help="Radius of the wire", default=0.01)
	parser.add_argument('-d', '--debug', help="Better debugging in code", action="store_true")
	parser.add_argument('-p', '--plot', help="Display plots", action="store_true")
	
	global args 
	args = parser.parse_args()
	
	N = args.num
	Im = args.current
	a = args.radius

	lamda = l*4.0 #Wavelength in meters
	f = c/lamda #Frequency in Hz
	dz = l/N #Spacing of current samples
	k = 2*pi/lamda #Wave number
	z = np.linspace(-N*dz,N*dz,int(2*N+1)) #Points in antenna
	u = np.delete(z,N)[1:-1] #Locations of Unknown currents
	I = np.zeros(2*N+1) #Current vector corresponding to z

	#Assigning the boundary currents
	I[N] = Im
	I[np.where(abs(z)==N*dz)] = 0
	
	J = np.zeros(2*N-2) #Current vector corresponding to u
	
	M = ampMat(len(J),a) #Ampere law matrix for H computed at r=a

	debug({"z vector":z, "I":I, "J":J, "M":M})

	Rz = mkDist(z,a) #Distance from locations of all samle points
	Ru = mkDist(u,a) #Distance from locations of unknown currents	
	Rn = np.delete(Rz[N],N)[1:-1] #Distance from Feed-point
	
	debug({"Rz":Rz, "Ru":Ru, "Rm":Rn})

	P = np.divide(mu0*exp(-1j*k*Ru)*dz/(4*np.pi),Ru)
	PB = np.divide(mu0*exp(-1j*k*Rn)*dz/(4*np.pi),Rn)
	
	debug({"P*1e8":(P*1e8).round(2), "PB":PB})

	Q = -1*P*(a/mu0)*(np.divide(-1j*k,Ru)+np.divide(-1,np.square(Ru)))
	QB = -1*PB*(a/mu0)*(np.divide(-1j*k,Rn)+np.divide(-1,np.square(Rn)))
	
	debug({"Q":Q, "QB":QB})
	
	J = np.squeeze(np.matmul(inv(M-Q),QB[:,np.newaxis])*Im)
	
	#Add boundary currents to J and store it in I
	I[1:(N)] = real(J[0:(N-1)]) #Lower half dipole
	I[(N+1):(2*N)] = real(J[N-1:2*N-2]) #Upper half dipole

	#Calculate the assumed sinusoidal variation
	I_true = Im*sin(k*(l-abs(z)))
	
	#Prints
	debug({"P*1e8":(P*1e8).round(2),"True current": I_true,"Estimated current":I})
	
	fig1, ax1 = plt.subplots()
	ax1.plot(z, I, linewidth=1.5, label=r'Estimated current')
	ax1.plot(z, I_true, linewidth=1.5, label=r'Assumed sinusoidal current')
	ax1.set_title(r'Current Distribution in a thin-wire Dipole antenna')
	ax1.set_xlabel(r'z (meters)$\longrightarrow$')
	ax1.set_ylabel(r'Current (amps)$\longrightarrow$')
	ax1.legend(loc='best', shadow=True, framealpha=1)
	#ax1.set_xlim(0, 1e-3)
	ax1.grid()
	fig1.savefig(PLOTS+'Figure1.png')
	
	if(args.plot):
		plt.show()
	
	print("The RMSE Error in the estimation is {}".format(np.square(I-I_true).mean(0).round(2)))

	if(not args.plot and not args.debug):
		print("For user interaction, Kindly run : python3 EE20B044.py --help")
main()
