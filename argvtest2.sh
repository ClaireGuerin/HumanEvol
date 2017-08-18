#!/bin/bash
while IFS=',' read -r line || [[ -n $line ]]; do
    #echo $line
	python testargv2.py $line
done < $1