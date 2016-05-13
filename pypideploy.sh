#!/bin/bash

set -x

echo 'PyPI deploynment script'

srv=$1
PW=$2
shift
shift

cat  >~/.pypirc <<-EOF 
[distutils]
index-servers =
    pypi
    pypitest

[pypi]
repository = https://pypi.python.org/pypi
username:jochym
password:$PW

[pypitest]
repository = https://testpypi.python.org/pypi
username:jochym
password:$PW

EOF

cat ~/.pypirc

echo $PW

python setup.py register -r $srv
python setup.py $* upload -r $srv
