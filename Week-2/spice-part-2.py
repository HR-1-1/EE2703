#######################################################
#Title   : spice-part-2.py
#Author  : Harish R EE20B044
#Date    : Feb 2 2022
#Purpose :
#Inputs  :
########################################################

import sys
import re
import numpy as np
from spice_part_1 import file_parser

class resistor:
	def __init__(self, name, node1, node2, value):
		self.name = name
		self.node1 = node1
		self.node2 = node2
		self.value = value
	
class ind_current_src:
	def __init__(self, name, node1, node2, value):
		self.name = name
		self.node1 = node1
		self.node2 = node2
		self.value = value
	
class ind_voltage_src:
	def __init__(self, name, node1, node2, value):
		self.name = name
		self.node1 = node1
		self.node2 = node2
		self.value = value

def nodes(elements):
	
	node_list = []
	i=1
	for element in elements:	
		node_list.extend([element.node1, element.node2])
	
	node_list = list(set(node_list))

	if "GND" not in node_list:
		print("Ground node not defined")
		sys.exit(0)
	else:
		node_list.remove("GND")

	node_dict = {y : x for x,y in enumerate(node_list)}
	return node_dict

def main():
	ckt_def = file_parser()
	elements = []
	for line in ckt_def:
		comp = line.split('#')[0].split()
		if comp[0][0]=='R':
			elements.append(resistor(*comp))
		elif comp[0][0]=='I':
			elements.append(ind_current_src(*comp))
		elif comp[0][0]=='V':
			elements.append(ind_voltage_src(*comp))
		else :
			print("sorry")
	
	node_dict=nodes(elements)
	print(node_dict)	
main()		
			

