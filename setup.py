from setuptools import setup, find_packages

with open("README_PACKAGE.md", "r") as f:
    desc = f.read()

setup(name="gsnetact", 
      
      version="0.0.3",

      description="The GSNetAct Python Package.",
        
      long_description=desc,
    
      long_description_content_type="text/markdown",

      packages=find_packages(),
            
      )
      
