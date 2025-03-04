from setuptools import setup, find_packages

setup(name = "pathWayScoring", 
      
      version = "0.0.1",

      description = "The pathway scoring algorithm as a python package.",
    
      package_dir = {"" : "pathwayScoring"},

      packages= find_packages(where="pathwayScoring"),
            
      )
      
