#!/bin/bash

VER=`(cd .. ; python -c "from setuptools_scm import get_version ; print(get_version())")`

echo "Building version: $VER"

sed "s/version: XXXX/version: \"$VER\"/g" <meta.yaml.tmpl >meta.yaml

sed "s/version = XXXX/version = \"$VER\"/g" <../docs/source/conf.py.tmpl  >../docs/source/conf.py

echo $VER >VERSION

