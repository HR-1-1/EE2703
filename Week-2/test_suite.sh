#!/bin/bash

if [ -z $1 ]; then
        echo "Filename is empty. Correct Usage : ${0##*/} <python-file>"
        exit 0
fi

for filename in $PWD/test_suite/*.netlist; do
	echo "Running test case $filename"
	python3 $1 $filename
	echo "-------------------"
done
