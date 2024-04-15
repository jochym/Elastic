# -*- coding: utf-8 -*-

from __future__ import annotations
from setuptools import setup
from setuptools_scm import ScmVersion
from setuptools_scm.version import no_local_version


def clean_version(version: ScmVersion) -> str:
    from setuptools_scm.version import guess_next_version
    from setuptools_scm.version import only_version
    
    return version.format_next_version(#guess_next_version,
                                        only_version, 
                                        "{guessed}.{distance}")


setup(use_scm_version={"version_scheme": clean_version, 
                        "local_scheme": no_local_version})

