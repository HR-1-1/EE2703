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

	print("Maximum deviation in exp coefficients : ", np.max(exp_dev))
	print("Maximum deviation in cos_cos coefficients : ", np.max(ccos_dev))	
