#!/bin/bash

VER=`(cd .. ; python -c "from setuptools_scm import get_version ; print(get_version())")`

echo $VER

sed "s/version: XXXX/version: $VER/g" <meta.yaml.tmpl >meta.yaml

