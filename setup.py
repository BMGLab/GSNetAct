"""Shim for tooling that still invokes setup.py directly.

All packaging metadata lives in ``pyproject.toml``.  The previous version read
``README_PACKAGE.md``, which is git-ignored and absent from the repository, so
``pip install .`` failed on a clean clone with ``FileNotFoundError``.
"""

from setuptools import setup

setup()
