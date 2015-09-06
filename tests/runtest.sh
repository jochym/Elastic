#!/bin/bash

mkdir tmp
cd tmp
ln -s ../parcalc ../elastic ../test_*.py .

if [ ".$1" = "." ] ; then
	for t in test_*.py ; do 
		python $t
	done
else
	python test_$1.py
fi

cd ..
rm -rf tmp

