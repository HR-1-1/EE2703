# Spice Netlist Part 1
## Context
A spice netlist consists of lines of the following forms:\
	`....
	.circuit
	name n1 n2 value # comment
	name n1 n2 n3 n4 value # comment
	name n1 n2 vname value # comment
	.end
	....`
The circuit definition consists of lines, each of which defines a branch of the circuit.
## Assignment
- Accept the name of netlist file as commandline.
- Analyze the tokens and determine the from node, the to node, the type of element and
the value.
- Traverse the circuit definition from last element to first and print out each line with
words in reverse order.
- Check for all corner cases.
## Usage Instructions
Clone the repo and execute the following commands in order. [This assumes you're in a environment with access to bash]
- `cd EE2703/WEEK-1`
- `chmod +x test_suite.sh`
- `./test_suite.sh spice-part-1.py`

