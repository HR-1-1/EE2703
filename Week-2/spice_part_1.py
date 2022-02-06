##############################################################
#Title   : spice-part-1.py
#Author  : Harish R EE20B044
#Date    : Jan 26 2022
#Purpose : 1. Analyze a Netlist file
#		   2. Traverse the circuit definition from last element \
#		      to first and print each line with words reversed
#Inputs  : A valid .netlist file [Given as a CommandLine input]
################################################################

import sys

CIRCUIT = '.circuit'
END = '.end'
AC = '.ac'
def netlist_analyzer(ckt_def):
	
	'''
	The Rules
	1. Neglect All lines till we meet a dot command
	2. .circuit starts the circuit definition and .end ends it
	3. Anythong that follows the # charecter in a line is comment
	4. A valid netlist file has onyy one circuit definiton
	5. Allowed Elements : R L C V I E[VCVS] G[VCCS] H[CCVS] F[CCCS] 
	'''
	#print("Analysis")

	# elems is a dictionary of elements and their symbols
	elems = {"R": "Resistor",
			"L": "Inductor",
			"C": "Capacitor",
			"V": "Independent Voltage Source",
			"I": "Independent Current Source",
			"E": "VCVS", 
			"G": "VCCS",
			"H": "CCVS", 
			"F": "CCCS"}
	
	for line in ckt_def:
		words = line.split('#')[0].split() # words is a list containing each entry of the MNA matrix
		source = words[0][0] # source contains the letter code indicating type of element
		
		print("Name of the Element: ", elems[source])
		print("\n".join(["Node {}: {}".format(i+1,word) for i, word in enumerate(words[1:len(words)-1])]))
		print("Value: {}".format(words[-1]))
		print('--'*5)
	return None

def reverser(ckt_def):
	
	print("\nReversal\n")
	for line in reversed([' '.join(reversed(line.split('#')[0].split())) for line in ckt_def]):
		print(line) 
		
	print('*'*10+'\n')

def file_parser():
	
	if len(sys.argv) != 2:
		print('\nUsage: python3 %s <inputfile>' % sys.argv[0])
		sys.exit()
	
	try:
		with open(sys.argv[1]) as f:
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
			
			opdir=[]
			
			for line in lines[end:]:
				if AC == line.split('#')[0].strip():
					opdir.append(line)
			
		return lines[start+1:end], opdir

	except IOError : # catch incorrect filename error
		print('Invalid file')
		sys.exit()

#ckt_def = file_parser()
#reverser(ckt_def)
#netlist_analyzer(ckt_def)
