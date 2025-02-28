#!/bin/bash
# run_examples.sh
# Batch file to run all DeSiO-Examples. on Luis Cluster
clear
mypath="$(pwd)"
rm "output_l.log" 
echo "Start running DeSiO - examples in " $mypath >> $mypath/output_l.log
# Loop over directories in DeSiO_Example
for f1 in $mypath/*; do
	# if file is directory then
    	if [ -d "$f1" ]; then
		for f2 in $f1/*; do
			if [ -d "$f2" ]; then
				f3=$f2
				cd "$f3"
				cp $mypath/run_cluster_l.sh $f3
				# send job to cluster, defined in run_cluster.sh
				sbatch run_cluster_l.sh
			fi
		done
	fi
done
