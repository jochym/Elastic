#!/bin/bash

set -x

echo 'PyPI deploynment script'

srv=$1
shift

cat  >~/.pypirc <<-EOF 
[distutils]
index-servers =
    pypi
    pypitest

[pypi]
repository = https://pypi.python.org/pypi
username:jochym
password:$PASSWORD

[pypitest]
repository = https://testpypi.python.org/pypi
username:jochym
password:$TESTPASSWORD

EOF

cat ~/.pypirc

python setup.py register -r $srv
python setup.py $* upload -r $srv
