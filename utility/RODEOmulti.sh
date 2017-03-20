#!/bin/bash

echo "RODEO with multiple inputs"

for a in input/*.inp; do

	b=${a#*/}
	c=${b%.inp}

#	echo $a
#	echo $b
#	echo $c
	echo "Executing: ./rodeo22.pl -l $a -o output/$c.html -p hmm -c config/thiopeptide.conf -csva output/$c.csv -x"
	./rodeo22.pl -l $a -o output/$c.html -p hmm -c config/thiopeptide.conf -csva output/$c.csv -x

done
