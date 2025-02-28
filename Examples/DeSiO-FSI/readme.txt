============================================================================
How to run all FSI examples on Luis-Cluster.
============================================================================
	Please follow the steps:

	1. open command shell terminal
	2. change working directory to current directory: ~$ cd $path
	3. run all examples on Cluster by entering: ~$ sbatch run_examples_FSI_l.sh
	4. if all examples are ready, evaluate results: ~$ sh run_evaluate_FSI_l.sh
	5. compare results written in output_l.log with output_ref_l.log

============================================================================
How to run all FSI examples on windows.
============================================================================
	Please follow the steps:

	1. copy DeSiO.bat in current directory
	2. start running examples: run_examples_FSI_w.bat
	3. results will be written in output_w.log 
	4. compare results with output_ref_w.log