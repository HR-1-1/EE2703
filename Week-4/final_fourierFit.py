"""
Title   : fourierFit.py
Author  : Harish R EE20B044
Date    : Feb 23 2022
Purpose : 1. Linear fit a curve to its fourier series
		   2. Study the difference in the best fit and integration approach
Inputs  : Sit back and relax
"""

import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.pyplot import clabel
from scipy import special as sp
from scipy.linalg import lstsq
from scipy.integrate import quad
PATH = "./plots/"

"""

# Plotting the Function computed using Fourier Coefficients
semilogy(x, exp_fn_fourier, 'ro', label="Function using Fourier Coefficients")
ylim([pow(10, -1), pow(10, 4)])
legend()

grid()
show()


plot(x, coscos_fn_fourier, 'ro', label="Function using Fourier Coefficients")
legend(loc='upper right')

plt.axis([-5, 5, -0.5, 2])
grid()
show()
"""

def exp(x):
    return np.exp(x)

def ccos(x):
    return np.cos(np.cos(x))

def an(x, k, f):
    return f(x)*np.cos(k*x)


def bn(x, k, f):
    return f(x)*np.sin(k*x)

def get_cf(f):

    cf = []
    cf.append((quad(f, 0, 2*np.pi)[0])/(2*np.pi))
    for i in range(1, 26):
        cf.append((quad(an, 0, 2*np.pi, args=(i, f))[0])/np.pi)
        cf.append((quad(bn, 0, 2*np.pi, args=(i, f))[0])/np.pi)

    return cf

def get_cf_mat(Npts, Ncfs , x):
    M = np.zeros((Npts, Ncfs))  # allocate space for A
    M[:, 0] = 1  # col 1 is all ones
    for k in range(1, int((Ncfs+1)/2)):
        M[:, 2*k-1] = np.cos(k*x)  # cos(kx) column
        M[:, 2*k] = np.sin(k*x)  # sin(kx) column
    return M

def get_fourier_fn(c):
    A = matrix_create(400, 51, x)
    f = A.dot(c)
    return f

def lstsq_fit(fun, N):
	x = np.linspace(0, 2*np.pi, N)[:-1]
	b = fun(x)
	A = get_cf_mat(N-1, 51, x)
	c = lstsq(A, b)[0]
	est = A.dot(c)
	return est, c

def main():

	x = np.linspace(-2*np.pi, 4*np.pi, 400)
	x_ = np.linspace(0, 2*np.pi, 400)
	exp_cf = get_cf(exp)
	ccos_cf = get_cf(ccos)

	exp_fit, exp_cf_fit = lstsq_fit(exp, 401)
	ccos_fit, ccos_cf_fit = lstsq_fit(ccos, 401)

	exp_dev = np.abs(exp_cf_fit - exp_cf)
	ccos_dev = np.abs(ccos_cf_fit - ccos_cf)
	
	fig1, ax = plt.subplots(figsize=(8,8))
	ax.semilogy(x, exp(x), 'k', label="Original Function")
	ax.semilogy(x, exp(x % (2*np.pi)), 'b--', label="Expected Function from fourier series")
	ax.semilogy(x_, exp_fit, 'go', label="Function from Least Square fit Co-efficients")
	ax.legend()
	ax.set_title(r'Figure 1 : Plot of $e^{x}$')
	ax.set_xlabel(r'x$\longrightarrow$')
	ax.set_ylabel(r'$e^{x}\longrightarrow$')
	ax.grid(True)
	fig1.savefig(PATH + "Figure1.png")

	fig2, ax = plt.subplots(figsize=(8,8))
	ax.plot(x, ccos(x), 'k', linewidth=2, label="Original Function")
	ax.plot(x, ccos(x % (2*np.pi)), 'b--', label="Expected Function from fourier series")
	ax.plot(x_, ccos_fit, 'go', label="Functions from Least Square fit Co-efficients")
	ax.legend(loc='upper right')
	ax.set_title(r"Figure 2 : Plot of $\cos(\cos(x))$")
	ax.set_xlabel(r'x$\longrightarrow$')
	ax.set_ylabel(r'$\cos(\cos(x)\longrightarrow$')
	ax.grid(True)
	fig2.savefig(PATH + "Figure2.png")
	
	# Computing function using fourier coeff
	#exp_fn_fourier = compute_fn(exp_cf)
	#coscos_fn_fourier = compute_fn(ccos_cf)

	fig3, ax = plt.subplots(figsize=(8,8))
	ax.semilogy(np.abs(exp_cf[1::2]), 'yo', label=r"$a_{n}$ using Integration")
	ax.semilogy(np.abs(exp_cf[2::2]), 'ro', label=r"$b_{n}$ using Integration")
	ax.semilogy(np.abs(exp_cf_fit[1::2]), 'go', label=r"$a_{n}$ using Least Squares")
	ax.semilogy(np.abs(exp_cf_fit[2::2]), 'bo', label=r"$b_{n}$ using Least Squares")
	ax.legend()
	ax.set_title("Figure 3 : Fourier coefficients of $e^{x}$ (semi-log)")
	ax.set_xlabel(r'n$\longrightarrow$')
	ax.set_ylabel(r'Magnitude of coeffients$\longrightarrow$')
	ax.grid(True)
	fig3.savefig(PATH + 'Figure3.png')
	
	fig4, ax = plt.subplots(figsize=(8,8))
	ax.loglog(np.abs(exp_cf[1::2]), 'yo', label=r"$a_{n}$ using Integration")
	ax.loglog(np.abs(exp_cf[2::2]), 'ro', label=r"$b_{n}$ using Integration")
	ax.loglog(np.abs(exp_cf_fit[1::2]), 'go', label=r"$a_{n}$ using Least Squares")
	ax.loglog(np.abs(exp_cf_fit[2::2]), 'bo', label=r"$b_{n}$ using Least Squares")
	ax.legend(loc='upper right')
	ax.set_title("Figure 4 : Fourier coefficients of $e^{x}$ (Log-Log)")
	ax.set_xlabel(r'n$\longrightarrow$')
	ax.set_ylabel(r'Magnitude of coeffients$\longrightarrow$')
	ax.grid(True)
	fig4.savefig(PATH + 'Figure4.png')
	

	fig5, ax = plt.subplots(figsize=(8,8))
	ax.semilogy(np.abs(ccos_cf[1::2]), 'yo', label=r"$a_{n}$ using Integration")
	ax.semilogy(np.abs(ccos_cf[2::2]), 'ro', label=r"$b_{n}$ using Integration")
	ax.semilogy(np.abs(ccos_cf_fit[1::2]), 'go', label=r"$a_{n}$ using Least Squares")
	ax.semilogy(np.abs(ccos_cf_fit[2::2]), 'bo', label=r"$b_{n}$ using Least Squares")
	ax.legend()
	ax.set_title("Figure 5 : Fourier coefficients of $\cos(\cos{x})$ (semi-log)")
	ax.set_xlabel(r'n$\longrightarrow$')
	ax.set_ylabel(r'Magnitude of coeffients$\longrightarrow$')
	ax.grid(True)
	fig5.savefig(PATH + 'Figure5.png')
	

	fig6, ax = plt.subplots(figsize=(8,8))
	ax.loglog(np.abs(ccos_cf[1::2]), 'yo', label=r"$a_{n}$ using Integration")
	ax.loglog(np.abs(ccos_cf[2::2]), 'ro', label=r"$b_{n}$ using Integration")
	ax.loglog(np.abs(ccos_cf_fit[1::2]), 'go', label=r"$a_{n}$ using Least Squares")
	ax.loglog(np.abs(ccos_cf_fit[2::2]), 'bo', label=r"$b_{n}$ using Least Squares")
	ax.legend(loc='upper right')
	ax.set_title("Figure 4 : Fourier coefficients of $\cos(\cos{x})$ (Log-Log)")
	ax.set_xlabel(r'n$\longrightarrow$')
	ax.set_ylabel(r'Magnitude of coeffients$\longrightarrow$')
	ax.grid(True)
	fig6.savefig(PATH + 'Figure6.png')
	
	print("Maximum deviation in exp coefficients : ", np.max(exp_dev))
	print("Maximum deviation in cos_cos coefficients : ", np.max(ccos_dev))
	
	plt.show()

main()
