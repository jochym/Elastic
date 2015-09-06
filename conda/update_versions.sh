#!/bin/bash

VER=`(cd .. ; python -c "from setuptools_scm import get_version ; print(get_version())")`

echo "Building version: $VER"

sed "s/version: XXXX/version: \"$VER\"/g" <meta.yaml.tmpl >meta.yaml

sed "s/version = XXXX/version = \"$VER\"/g" <../doc/source/conf.py.tmpl  >../doc/source/conf.py
