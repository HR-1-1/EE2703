########################################################
#Title   : xtraSpice.py
#Author  : Harish R EE20B044
#Date    : Feb 2 2022
#Purpose : Circuit solver
#Inputs  : A dot Netlist file
########################################################

import sys
import re
import numpy as np
import warnings
from math import pi

CIRCUIT = '.circuit'
END = '.end'
AC = '.ac'

def P2R(A,phi):
	return A * ( np.cos(phi) + np.sin(phi)*1j )   

def helper_mna(x,m,n1,n2,v1,v2):

	x[n1][v1] += m
	x[n1][v2] += -m
	x[n2][v1] += -m
	x[n2][v2] += m

def get_obj(name):
	for elem in solver.elements:
		if elem.name == name:
			return elem

class passive:

	def __new__(cls, name, n1, n2, val):
		
		if float(val)<=0:
			warnings.warn(" Value of impedence {} is zero".format(name))
		
		if n1==n2:
			warnings.warn("Resistor {} is connected across the identitical nodes".format(name))
			return None

		instance = super().__new__(cls)
		instance.name = name
		instance.n1 = n1
		instance.n2 = n2
		instance.val = float(val)	
		
		return instance
	
	def fill_mna(self, mna_G, mna_ind, unk, freq=None):
		
		if self.name[0] == 'R':
			val = 1/self.val
		elif self.name[0] == 'C':
			val = 2j*pi*freq*self.val
		elif self.name[0] == 'L':
			val =  -1j/(2*pi*freq*self.val)

		helper_mna(mna_G, val, unk[self.n1][1], 
					unk[self.n2][1], unk[self.n1][1], unk[self.n2][1])

class ind_current_src:

	def __init__(self, name, n1, n2,mode,val,phase=None):
		self.name = name
		self.n1 = n1
		self.n2 = n2	
		self.mode = mode
		if mode == 'dc':
			self.val = float(val)
		elif mode == 'ac':
			self.val = P2R(float(val)/2, float(phase))
		else:
			print("Enter a proper mode")

	def fill_mna(self, mna_G, mna_ind, unk, freq=None):
		mna_ind[unk[self.n1][1]] += -self.val
		mna_ind[unk[self.n2][1]] += self.val

class ind_voltage_src:

	def __init__(self, name, n1, n2, mode, val, phase=None):
		self.name = name
		self.n1 = n1
		self.n2 = n2
		self.mode = mode
		if mode == 'dc':
			self.val = float(val)
		elif mode == 'ac':
			self.val = P2R(float(val)/2, float(phase))
		
	def fill_mna(self, mna_G, mna_ind, unk, freq=None):
		
		mna_G[unk[self.n1][1]][unk[self.name][1]] += 1
		mna_G[unk[self.n2][1]][unk[self.name][1]] += -1
		mna_G[unk[self.name][1]][unk[self.n1][1]] += 1
		mna_G[unk[self.name][1]][unk[self.n2][1]] += -1

		mna_ind[unk[self.name][1]] += self.val

class vcvs:

	def __init__(self, name, n1, n2, n3, n4, val):
		self.name = name
		self.n1 = n1
		self.n2 = n2
		self.n3 = n3
		self.n4 = n4
		self.val = float(val)
	
	def fill_mna(self, mna_G, mna_ind, unk, freq=None):
		mna_G[unk[self.n1][1]][unk[self.name][1]] += 1
		mna_G[unk[self.n2][1]][unk[self.name][1]] += -1
		mna_G[unk[self.name][1]][unk[self.n1][1]] += 1
		mna_G[unk[self.name][1]][unk[self.n2][1]] += -1
		mna_G[unk[self.name][1]][unk[self.n3][1]] += -self.val
		mna_G[unk[self.name][1]][unk[self.n4][1]] += self.val
	
class vccs:
	
	def __init__(self, name, n1, n2, n3, n4, val):
		self.name = name
		self.n1 = n1
		self.n2 = n2
		self.n3 = n3
		self.n4 = n4
		self.val = float(val)
	
	def fill_mna(self, mna_G, mna_ind, unk, freq=None):
		helper_mna(mna_G, self.val, unk[self.n1][1], 
					unk[self.n2][1], unk[self.n3][1], unk[self.n4][1])

class ccvs:
	
	def __init__(self, name, n1, n2, v0, val):
		self.name = name
		self.n1 = n1
		self.n2 = n2
		self.v0 = v0
		self.val = float(val)

	def fill_mna(self, mna_G, mna_ind, unk, freq=None):
		
		v = get_obj(self.v0)
		mna_G[unk[self.n1][1]][unk[self.name][1]] += 1
		mna_G[unk[self.n2][1]][unk[self.name][1]] += -1

		mna_G[unk[self.name][1]][unk[self.n1][1]] += 1
		mna_G[unk[self.name][1]][unk[self.n2][1]]  += -1

		mna_G[unk[self.name][1]][unk[v.name][1]] += -self.val

class cccs:
	
	def __init__(self, name, n1, n2, v0, val):
		self.name = name
		self.n1 = n1
		self.n2 = n2
		self.v0 = v0
		self.val = float(val)

	def fill_mna(self, mna_G, mna_ind, unk, freq=None):
		
		v = get_obj(self.v0)
		mna_G[unk[self.n1][1]][unk[v.name][1]] += self.val
		mna_G[unk[self.n2][1]][unk[v.name][1]] += -self.val

class xtraSpice: # netlist solver class
	
	def __init__(self, netlist):
		self.netlist = netlist
		self.mode = None
		self.freq = None
		self.file_parser()
	
	def get_mode(self, opdir):
		 
		if opdir.split('#')[0].strip().split()[0] == AC:
			self.mode = 'ac'
			self.freq = float(opdir.split('#')[0].strip().split()[-1])
		else:
			pass
	
	def file_parser(self): # parses the netlist file and checks its validity
	
		try:
			with open(self.netlist) as f:
				lines = f.readlines()
				start = -1; end = -2

				for line in lines:	# extracting circuit definition start and end lines

					if CIRCUIT == line.split('#')[0].strip():
						start = lines.index(line)
					elif END == line.split('#')[0].strip():	
						end = lines.index(line)
						break		
				
				if start>=end or start*end<0:	# validating circuit block
					print("Invalid circuit definition")
					sys.exit(0)
				
				for line in lines[end:]:
					if AC == line.split()[0]:
						self.get_mode(line)
				
			self.ckt_def = lines[start+1:end]
			
		except IOError : # catch incorrect filename error
			print('Invalid file')
			sys.exit
	
	def element_extracter(self): # creates objects for elements and stores them in a list

		self.elements = []
		
		# elems is a dictionary of elements and classes
		elems = {"R": passive,
			"L": passive,
			"C": passive,
			"V": ind_voltage_src,
			"I": ind_current_src,
			"E": vcvs, 
			"G": vccs,
			"H": ccvs, 
			"F": cccs}
	
		for line in self.ckt_def:
			comp = line.split('#')[0].split()
			 
			try:
				self.elements.append(elems[comp[0][0]](*comp))
			except Exception as e:
				if e is KeyError:
					print("Kindly follow the naming convention for Elements")
				else:
					pass
		
	def get_nodes(self): # get the nodes and form the unknown vector
		
		node_list = []

		for element in self.elements:	
			node_list.extend([element.n1, element.n2])
		
		node_list = list(set(node_list))

		if "GND" not in node_list:
			print("Ground node not defined")
			sys.exit(0)
		else:
			node_list.remove("GND")

		self.unk = {y : ['NV',x] for x,y in enumerate(node_list)}
			
		for elem in self.elements:
			if isinstance(elem, (ind_voltage_src, vcvs, ccvs)):
				self.unk[elem.name] = ['VSCA',len(self.unk)]
			
		self.unk['GND'] = ['GND',len(self.unk)]

	def mna_solver(self): # solve the mna matrix 
		
		self.mna_G = np.zeros((len(self.unk),len(self.unk))).astype(complex)
		self.mna_ind = np.zeros(len(self.unk)).astype(complex)

		for elem in self.elements:
			elem.fill_mna(self.mna_G, self.mna_ind, self.unk, self.freq)
		
		try:
			self.mna_res = np.linalg.solve(self.mna_G[:-1,:-1], self.mna_ind[:-1])
		except np.linalg.LinAlgError as e:
			print("Sorry! This ciruit raises the following error : {}".format(e))
	
	def solve(self):
		self.element_extracter()
		self.get_nodes()
		self.mna_solver()
	
	def display(self, debug =False): # display the result
	
		abbrev = {'NV':"Node Voltage",
				'VSCA' : "Voltage Source current"}
		
		if debug==True:
			print(self.mna_G)
			print(self.mna_ind)
			print(self.unk)

		for un,idx in self.unk.items():
			try:
				if un != 'GND':	
					print("Value of {}, {} is {:.2E} {}".format(abbrev[idx[0]],un,self.mna_res[idx[1]],idx[0][-1]))
			except Exception as e:
				print("Error!", e)

			
if len(sys.argv) != 2:
	print('\nUsage: python3 %s <inputfile>' % sys.argv[0])
	sys.exit()

solver = xtraSpice(sys.argv[1])
solver.solve()
solver.display()
