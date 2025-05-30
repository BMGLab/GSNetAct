from setuptools import setup, find_packages

with open("README.md","r") as f:
    desc = f.read()

setup(name="gsnetact", 
      
      version="0.0.1",

      description="The pathway scoring algorithm as a python package.",
        
      long_description=desc,
    
      long_description_content_type="text/markdown",

      packages=find_packages(),
            
      )
      
