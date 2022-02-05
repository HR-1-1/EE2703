########################################################
#Title   : spice_part_2.py
#Author  : Harish R EE20B044
#Date    : Feb 2 2022
#Purpose : A circuit solver
#Inputs  : A Netlist file
########################################################

import sys
import re
import numpy as np
import warnings

from spice_part_1 import file_parser

def helper_mna(x,m,n1,n2,v1,v2):

	x[n1][v1] = m
	x[n1][v2] = -m
	x[n2][v1] = -m
	x[n2][v2] = m

class resistor:
	def __new__(cls, name, n1, n2, val):
		
		if float(val)<=0:
			warnings.warn(" Value of resistor {} is zero".format(name))
			return None
		
		if n1==n2:
			warnings.warn("Resistor {} is connected across the identitical nodes".format(name))
			return None

		instance = super().__new__(cls)
		instance.name = name
		instance.n1 = n1
		instance.n2 = n2
		instance.val = float(val)
	
		return instance
	
	def fill_mna(self, mna_G, mna_ind, unk):
		helper_mna(mna_G, 1/self.val, unk[self.n1][1], 
					unk[self.n2][1], unk[self.n1][1], unk[self.n2][1])
	
class ind_current_src:
	def __init__(self, name, n1, n2, val):
		self.name = name
		self.n1 = n1
		self.n2 = n2
		self.val = float(val)

	def fill_mna(self, mna_G, mna_ind, unk):
		mna_ind[unk[self.n1][1]] = self.val
		mna_ind[unk[self.n2][1]] = self.val

class ind_voltage_src:
	def __init__(self, name, n1, n2, val):
		self.name = name
		self.n1 = n1
		self.n2 = n2
		self.val = float(val)
	
	def fill_mna(self, mna_G, mna_ind, unk):

		mna_G[unk[self.n1][1]][unk[self.name][1]] = 1
		mna_G[unk[self.n2][1]][unk[self.name][1]] = -1
		mna_G[unk[self.name][1]][unk[self.n1][1]] = 1
		mna_G[unk[self.name][1]][unk[self.n2][1]] = -1

		mna_ind[unk[self.name][1]] = self.val

class xtraSpice:
	
	def __init__(self, ckt_def):
		self.ckt_def = ckt_def

	def element_extracter(self):

		self.elements = []

		for line in self.ckt_def:
			comp = line.split('#')[0].split()
			if comp[0][0]=='R':
				self.elements.append(resistor(*comp))
			elif comp[0][0]=='I':
				self.elements.append(ind_current_src(*comp))
			elif comp[0][0]=='V':
				self.elements.append(ind_voltage_src(*comp))
			else :
				print("Kindly follow the naming convention for Elements")
	
	def get_nodes(self):
		
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
			if isinstance(elem, ind_voltage_src):
				self.unk[elem.name] = ['VSCA',len(self.unk)]
		
		self.unk["GND"] = ['GND',len(self.unk)]

	def mna_solver(self):
		
		self.mna_G = np.zeros((len(self.unk),len(self.unk)))
		self.mna_ind = np.zeros(len(self.unk))

		for elem in self.elements:
			elem.fill_mna(self.mna_G, self.mna_ind, self.unk)
		
		#print(mna_G, mna_ind)

		self.mna_res = np.linalg.solve(self.mna_G[:-1,:-1], self.mna_ind[:-1])
	
	def solve(self):
		self.element_extracter()
		self.get_nodes()
		self.mna_solver()
	
	def display(self):
	
		abbrev = {'NV':"Node Voltage",
				'VSCA' : "Voltage Source current",
				'CCA' : "Controlling current"}
		
		for un,idx in self.unk.items():
			try:
				print("Value of {}, {} is {} {}".format(abbrev[idx[0]],un,self.mna_res[idx[1]],idx[0][-1]))
			except Exception as e:
				pass

def main():
	ckt_def = file_parser()
	solver = xtraSpice(ckt_def)
	solver.solve()
	solver.display()

main()
