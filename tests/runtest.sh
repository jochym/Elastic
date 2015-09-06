#!/bin/bash

mkdir tmp
cd tmp
ln -s ../parcalc ../elastic ../test.py .
python test.py
cd ..
rm -rf tmp

