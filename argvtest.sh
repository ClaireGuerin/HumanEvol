#!/bin/bash
while IFS=',' read -r line || [[ -n $line ]]; do
    #echo $line
	python testargv.py $line
done < $1