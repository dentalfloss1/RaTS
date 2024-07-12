import glob

from setuptools import setup, Extension, Distribution, find_packages

module1 = Extension('RaTS/sbpl',['RaTS/sbplmodule.c'],
                     extra_link_args=["-lm","-lgsl","-lgslcblas"])

setup(name="RaTS",
      version="1.0",
      description="Radio Transient Simulations",
      author="Sarah Chastain",
      author_email="sarahichastain@gmail.com",
      packages=["RaTS"],
      scripts = glob.glob('scripts/*.py'),
      python_requires ='>3.6',
      install_requires = ["astropy","bitarray","colorama","cycler","fonttools","kiwisolver","matplotlib","numpy","packaging","Pillow","pyerfa","pyparsing","python-dateutil","PyYAML","scipy","six","tqdm"],
      test_suite="tests",
      ext_modules=[module1])


