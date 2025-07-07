from setuptools import setup, find_packages

with open("README_PACKAGE.md", "r") as f:
    desc = f.read()

setup(name="gsnetact",

      version="0.0.6",

      description="The GSNetAct Python Package.",

      long_description=desc,

      long_description_content_type="text/markdown",

      packages=find_packages(),

      entry_points={

          'console_scripts': [

              'makeGeneSets=gsnetact.Utils.makeJsonFile:makeJson_console'

              ]

          },

      )
