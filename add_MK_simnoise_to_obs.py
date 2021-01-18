import numpy as np
from argparse import ArgumentParser
from datetime import datetime

starttime = datetime.now().strftime('%y%m%d%Hh%Mm%Ss')
parser = ArgumentParser()
parser.add_argument('obsfile', help='supply a list of observations')
args = parser.parse_args()
observations = np.loadtxt(args.obsfile,dtype={'names': ('dateobs', 'duration', 'field'), 'formats': ('U32','f8','U32')})
scalarfactor = 4e-4*np.sqrt(8)
MKnoises = np.array([np.random.normal(scalarfactor/np.sqrt(o), scale=0.08*scalarfactor/np.sqrt(o)) for o in observations['duration']]) #  VERY ROUGH approximation for noise from meerkat images 
with open('obsnoise'+starttime+'.txt', 'a+') as f:
        for date,duration,noise,field in zip(observations['dateobs'],observations['duration'],MKnoises,observations['field']):
                f.write("{},{},{},{}\n".format(date, duration,noise,field))       
        f.flush()
print('wrote obsnoise'+starttime+'.txt')