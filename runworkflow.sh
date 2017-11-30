#! /bin/bash

if [ ! -e part_I.py ]; then
	echo "Fatal Error: Workflow file missing (part_Ia.py)"
	exit
fi

if [ ! -e part_II.py ]; then
    echo "Fatal Error: Workflow file missing (part_IIa.py)"
    exit
fi

if [ ! -e config.json ]; then
	echo "Fatal Error: Workflow file missing (config.json)"
	exit
fi

if [ -e part_I.py -a -e part_II.py -a -e config.json ]; then
	echo "Status: All files present, beginning workflow"
	
	python3 part_I.py
	if (( $? == 0 )); then
	    echo "Fatal Error: Unable to execute part_I.py"
	    exit
	fi

	python3 part_II.py
	if (( $? == 0 )); then
	    echo "Fata Error: Unable to execute part_II.py"
	    exit
	fi
fi