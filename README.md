# transient-simulations

## Requirements

Either chromium or firefox with geckodriver (https://github.com/mozilla/geckodriver)

Python 3.6 or greater with the following libraries:
* numpy
* scipy
* bokeh
* tqdm 
* selenium

Should be platform independent

## Installing

1. Clone the repository.
2. Ensure that you have either chromium or firefox with geckdriver installed
3. Run ```python3 -m pip install pipenv --user``` to do a user install of pipenv (replacing python3 with whatever the proper alias is). 
4. cd into the simulations subdirectory and run ```pipenv install```

Optionally, you can try installing the above packages individually and create a virtual environment, but troubleshooting will be more challenging.


## Running the Simulation

1. Edit the config.ini file to your liking
2. Optionally specify a observation file (or fill out the observation parameters in the config.ini file)
3. Run simulation using:
``` pipenv run python3 simulate.py ```


## Adding lightcurves

Lightcurves are imported dynamically by calling whatever is in the lightcurve type field in the config.ini file. 
For example, if there is a lightcurve class file called "example.py"  then all one has to do is specify "example"
in the lightcurvetype variable. 

The structure of these files should be easy to copy by taking the existing lightcurves as examples. The procedure 
generally should be as follows:

1. Specify whether the lightcurve has definite edges
2. Define the earliest and latest critical times that the lightcurve can be simulated for. 
*Note: Lightcurves with a definite beginning must have a critical time at
the beginning of the lightcurve*
3. Specify function for the integrated flux
4. Specify functions for the lines of the expected probability of 1 

