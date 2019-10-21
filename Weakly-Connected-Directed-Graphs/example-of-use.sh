#!/bin/bash

if [ ! -e img ]; then mkdir img; fi;
if [ ! -e img/data_dbol ]; then mkdir img/data_dbol; fi;
if [ ! -e img/data_duol ]; then mkdir img/data_duol; fi;

#************************* GENERATION OF PH********************
echo "Generation of  Random Graphs with DUO the Algorithm";
	if [ -e data_duol/results_nmi.txt ]; then rm data_duol/results_nmi.txt; fi; 
	for i in 100 101 102 103 104 105; do #different seed
		./duol -c 30 -n 0 -s $i > /dev/null;  
	done;
  
echo "Generation of  Random Graphs with DBO the Algorithm";
	if [ -e data_dbol/results_nmi.txt ]; then rm data_dbol/results_nmi.txt; fi; 
	for i in 100 101 102 103 104 105; do #different seeds
		./dbol -c 30 -n 0 -s $i > /dev/null;  
	done;
  
