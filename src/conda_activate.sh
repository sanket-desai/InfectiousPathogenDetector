#!/bin/bash
medaka_env=$(conda info --envs | grep "medaka")
if [[ -z "$medaka_env" ]]; then
	echo "Error: Please provide a valid virtual environment. For a list of valid virtual environment, please see 'conda env list' "
	exit   
else 
	conda activate medaka
	medaka_variant -i $1 -f $2 -o $3 -t $4 -p -d
fi;
