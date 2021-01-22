import numpy as np
from argparse import ArgumentParser
from datetime import datetime, timedelta

FOV = 1.4 # degrees 
minint = 8 # minimum imaging time for meerkat is 8 seconds
starttime = datetime.now().strftime('%y%m%d%Hh%Mm%Ss')
parser = ArgumentParser()
parser.add_argument('obsfile', help='supply a list of observations')
parser.add_argument('--multiconfig',action='store_true')
args = parser.parse_args()
multiconfig = args.multiconfig
observations = np.loadtxt(args.obsfile, delimiter = ',', dtype={'names': ('dateobs', 'duration', 'ra', 'dec'), 'formats': ('U32','f8','f8', 'f8')})
scalarfactor = 4e-4*np.sqrt(8)
MKnoises = np.array([np.random.normal(scalarfactor/np.sqrt(o), scale=0.08*scalarfactor/np.sqrt(o)) for o in observations['duration']]) #  VERY ROUGH approximation for noise from meerkat images 
with open('obsnoise'+starttime+'.txt', 'a+') as f:
        for date,duration,noise,r,d in zip(observations['dateobs'],observations['duration'],MKnoises,observations['ra'],observations['dec']):
                f.write("{},{},{},{},{},{}\n".format(date, duration,noise, r, d, FOV))       
        f.flush()
if multiconfig:
    num_minims = int(np.sum([d//minint for d in observations['duration']]))
    # print(num_minims, num_minims2)
    min_obs = np.zeros(num_minims, dtype={'names': ('dateobs', 'duration', 'ra', 'dec'), 'formats': ('U32','f8','f8', 'f8')})
    min_obs['duration'] = minint
    for i in range(len(observations)):
        o = observations[i]
        for j in range(int(o['duration']//minint)):
            currentindex = int(np.sum([d//minint for d in observations['duration'][0:i]]) + j)
            min_obs['dateobs'][currentindex] = (datetime.fromisoformat(o['dateobs']) + timedelta(seconds=minint*j)).strftime("%Y-%m-%dT%H:%M:%S.%f+00:00")
            min_obs['ra'][currentindex] = o['ra']
            min_obs['dec'][currentindex] = o['dec']
    minMKnoises = np.array([np.random.normal(scalarfactor/np.sqrt(o), scale=0.08*scalarfactor/np.sqrt(o)) for o in min_obs['duration']])
    with open('obsnoise'+starttime+'.txt', 'a+') as f:
        for date,duration,noise,r,d in zip(min_obs['dateobs'],min_obs['duration'],minMKnoises,min_obs['ra'],min_obs['dec']):
                f.write("{},{},{},{},{},{}\n".format(date, duration,noise, r, d,FOV))       
        f.flush()
    datetimes = np.array([datetime.fromisoformat(date) for date in observations['dateobs']])
    datetimes.sort()
    start = 0
    # print(datetimes[start])
    # print(datetimes[start+1:] - datetimes[start])
    for i in range(8):
        try:
            sameday = np.where((datetimes[(start):] - datetimes[start])<timedelta(days=1))[0] + start
            start = max(sameday)+1
            daynoise = np.random.normal(scalarfactor/np.sqrt(np.sum(observations['duration'][sameday])), scale=0.08*scalarfactor/np.sqrt(np.sum(observations['duration'][sameday])))
            with open('obsnoise'+starttime+'.txt', 'a+') as f:
                for date,duration,r,d in zip(observations['dateobs'][sameday],observations['duration'][sameday],observations['ra'][sameday],observations['dec'][sameday]):
                    f.write("{},{},{},{},{},{}\n".format(date, duration,daynoise, r, d,FOV))       
                f.flush()
        except IndexError:
            break

    
    # for i in range(len(datetimes)):
        
    # print(datetimes[(([datetimes[1:] - datetimes[0]).total_seconds()/3600/24) < 24])
    # iterate, parse with datetime, add times using durations, bin larger, write file etc all that stuff

print('wrote obsnoise'+starttime+'.txt')