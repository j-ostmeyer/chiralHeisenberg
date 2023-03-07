#!/bin/bash

for syst in mh rhmc; do
	cd $syst
	rm -f sqrt_phi2.csv

	for dir in */*/; do
		l=$(echo "${dir%%/*}" | cut -d '_' -f 2)
		b=$(echo "${dir%/}" | cut -d '_' -f 3)

		if [ -e "${dir}fort.14" ]; then
			dat="${dir}fort.14"
		else
			dat="${dir}condmod"
		fi

		res=$(awk '{print $1}' "${dat}" | comp-avg -n | cut -f 2,4,5,7)

		echo "$l	$b	$res" >> sqrt_phi2.csv

	done

	cd ..
done
