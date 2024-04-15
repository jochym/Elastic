# -*- coding: utf-8 -*-

from __future__ import annotations
from setuptools import setup
from setuptools_scm import ScmVersion
from setuptools_scm.version import get_no_local_node

def current_version(version: ScmVersion) -> str:
    from setuptools_scm.version import guess_next_version
    
    return version.format_next_version(guess_next_version,
                                        "{tag}.{distance}")


setup(use_scm_version={"version_scheme": current_version, 
                        "local_scheme": get_no_local_node})

