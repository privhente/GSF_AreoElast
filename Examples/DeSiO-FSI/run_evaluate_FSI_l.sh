#!/bin/bash
# run_examples.sh
# Batch file to run all DeSiO-Aero Examples. on Luis Cluster
clear
mypath="$(pwd)"
rm "output_l.log"
filecheck="check.log"
echo "Start running DeSiO-FSI examples in " $mypath >> $mypath/output_l.log
# Loop over directories in DeSiO_Example
for f1 in $mypath/*; do
	# if file is directory then
    	if [ -d "$f1" ]; then
		for f2 in $f1/*; do
			if [ -d "$f2" ]; then
				f3=$f2
				cd "$f3"
				name=$(cat "$filecheck") #'cat $file' is assigned to the $name variable
				echo 'results from ' $f3 ': ' $name >> $mypath/output_l.log
			fi
		done
	fi
done
