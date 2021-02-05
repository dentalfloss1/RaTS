import configparser
import random
import time
import datetime
import numpy as np
import os
import glob
import argparse
import math
import warnings
from bokeh.plotting import figure, show, output_file
from bokeh.models import LinearColorMapper, SingleIntervalTicker, ColorBar, Title
from bokeh.io import export_png
import scipy.interpolate as interpolate
from tqdm import tqdm
from astropy import units as u 
from astropy.coordinates import SkyCoord
from scipy.special import binom
from concurrent import futures

def observing_strategy(obs_setup, det_threshold, nobs, obssens, obssig, obsinterval, obsdurations):
    if obs_setup is not None:
        tstart, tdur, sens, ra, dec  = np.loadtxt(obs_setup, unpack=True, delimiter = ',',
            dtype={'names': ('start', 'duration','sens', 'ra', 'dec'), 'formats': ('U32','f8','f8','f8','f8')})
        
        tstart = np.array([time.mktime(datetime.datetime.strptime(t, "%Y-%m-%dT%H:%M:%S.%f+00:00").timetuple())/(3600*24) for t in tstart])
        sortkey = np.argsort(tstart)
        obs = np.zeros((len(tstart),3))
        obs[:,0]+=tstart[sortkey]
        obs[:,1]+=tdur[sortkey]/3600/24 # convert into days
        obs[:,2]+=sens[sortkey]
        pointing = np.array([np.array([r,d]) for r,d in zip(ra[sortkey],dec[sortkey])])
    else: # Enter 'trial mode' according to specified cadence
        observations = []
        for i in range(nobs):
            tstart = (time.mktime(datetime.datetime.strptime("2019-08-08T12:50:05.0", "%Y-%m-%dT%H:%M:%S.%f").timetuple()) + 60*60*24*obsinterval*i)/(3600*24)    #in days at an interval of 7 days
            tdur = obsdurations
            sens = random.gauss(obssens, obssig) * det_threshold # in Jy. Choice of Gaussian and its characteristics were arbitrary.
            observations.append([tstart,tdur,sens])
        observations = np.array(observations,dtype=np.float64)
        # pointing = np.array([SkyCoord(ra=275.0913169*u.degree, dec=7.185135679*u.degree, frame='icrs') for l in observations])
        pointing = np.array([np.array([275.0913169,7.185135679]) for l in observations])
        # pointing[0:4]-=1
        pointing = pointing[observations[:,0].argsort()]
        obs = observations[observations[:,0].argsort()] # sorts observations by date
        FOV = np.array([1.5 for l in observations]) # make FOV for all observations whatever specified here, 1.5 degrees for example
        pointing[0:4,0]+=2 # Make an offset between some pointings
        pointFOV = np.zeros((len(observations),3))
        pointFOV[:,0:2] += pointing
        pointFOV[:,2] += FOV
        uniquepoint = np.unique(pointFOV,axis=0)
        uniquesky = SkyCoord(ra=uniquepoint[:,0],dec=uniquepoint[:,1], unit='deg', frame='fk5')
        print(uniquesky)
        ### Next two lines set up variables for defining region properties. Region array is of maximum theoretical length assuming no more than double overlapping fovs (triple or more never get computed ##
        numrgns = len(uniquepoint) 
        regions = np.zeros(np.uint32(numrgns + binom(numrgns,2)), dtype={'names': ('ra', 'dec','identity', 'area'), 'formats': ('f8','f8','U32','f8')})
        for i in range(numrgns): # Label the individual pointings. These regions are for example region 1 NOT 2 and 2 NOT 1 
            regions['identity'][i] = str(i)
            regions['ra'][i] = uniquepoint[i,0]
            regions['dec'][i] = uniquepoint[i,1]
            regions['area'][i] = (4*np.pi*np.sin(uniquepoint[i,2]*np.pi/180/2)**2)
            leftoff = i + 1
        for i in range(len(uniquesky)): # Label intersections: For example: 1 AND 2
            for j in range(i+1,len(uniquesky)):
                if uniquesky[i].separation(uniquesky[j]).deg <= (uniquepoint[i,2] + uniquepoint[j,2]):
                    d = uniquesky[i].separation(uniquesky[j]).rad
                    r1 = uniquepoint[i,2]*np.pi/180
                    r2 = uniquepoint[j,2]*np.pi/180
                    gamma = np.arctan((np.cos(r2)/np.cos(r1)/np.sin(d)) - (1/np.tan(d)))
                    pa = uniquesky[i].position_angle(uniquesky[j])
                    # https://arxiv.org/ftp/arxiv/papers/1205/1205.1396.pdf
                    # and https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
                    fullcone1 = 4*np.pi*np.sin(r1/2)**2 
                    cutchord1 = 2*(np.arccos(np.sin(gamma)/np.sin(r1)) - np.cos(r1)*np.arccos(np.tan(gamma)/np.tan(r1))) 
                    fullcone2 = 4*np.pi*np.sin(r2/2)**2
                    cutchord2 = 2*(np.arccos(np.sin(gamma)/np.sin(r2)) - np.cos(r2)*np.arccos(np.tan(gamma)/np.tan(r2))) 
                    
                    centerreg = uniquesky[i].directional_offset_by(pa, gamma*u.radian)
                    
                    regions['identity'][leftoff] = str(i)+'&'+str(j)
                    regions['ra'][leftoff] = centerreg.ra.deg
                    regions['dec'][leftoff] = centerreg.dec.deg
                    regions['area'][leftoff] = cutchord1 + cutchord2
                    leftoff+=1
    return obs, pointFOV, regions[regions['identity'] != '']
    

def generate_pointings(n_sources, pointFOV, regions):
    uniquepointFOV = np.unique(pointFOV,axis=0)
    maxFOV = np.max(pointFOV[:,2])
    rng = np.random.default_rng() 
    minra = min(uniquepointFOV[:,0])
    mindec = min(uniquepointFOV[:,1])
    maxra = max(uniquepointFOV[:,0])
    maxdec = max(uniquepointFOV[:,1])
    uniqueskycoord = SkyCoord(ra=uniquepointFOV[:,0]*u.deg, dec=uniquepointFOV[:,1]*u.deg, frame='fk5')
    ra_random = lambda n: (rng.random(int(n))*(min(maxra + maxFOV,360) - max(minra - maxFOV,0)) + max(minra - maxFOV,0)) * u.degree
    dec_random = lambda n: (rng.random(int(n))*(min(maxdec + maxFOV,90)  - max(mindec - maxFOV, -90)) + max(mindec - maxFOV, -90)) * u.degree
    rollforskycoord = lambda n: SkyCoord(ra=ra_random(n), dec=dec_random(n), frame='fk5')
    ## This algorithm makes a somewhat large boundary that should encompass all pointings and generate sources within the large boundary. ##
    ## Then it rejects the ones not in a FOV of an observations and "re-rolls the dice" until they are all in FOV ##
    c = rollforskycoord(n_sources)
    sourceinfo = np.zeros(n_sources, dtype={'names': ('ra', 'dec','identity', 'area'), 'formats': ('f8','f8','U32','f8')})
    accept = np.zeros(c.shape,dtype=bool)
    reject = np.logical_not(accept) 
    while True:
        for i in range(len(uniqueskycoord)):
            cond = np.zeros(len(c),dtype=bool)
            cond[reject] += c[reject].separation(uniqueskycoord[i]).deg <= (uniquepointFOV[i,2])
            sourceinfo['identity'][reject & cond] = str(i)+'&'
        sourceinfo['identity'][reject] = (np.char.rstrip(sourceinfo['identity'][reject], U'&'))
        for i in range(len(regions)):
            key = sourceinfo['identity'] == regions['identity'][i]
            sourceinfo['ra'][reject & key]+=regions['ra'][i]
            sourceinfo['dec'][reject & key]+=regions['dec'][i]
            sourceinfo['area'][reject & key]+=regions['area'][i]
            # print(key)
        accept = sourceinfo['identity'] != ''
        reject = np.logical_not(accept)
        print("Rejected sources: ",len(c[reject]))
        if len(c[reject])==0:
            break
        c.data.lon[reject] = ra_random(len(c[reject]))
        c.data.lat[reject] = dec_random(len(c[reject]))
        c.cache.clear()
    
    return sourceinfo

def generate_sources(n_sources, start_survey, end_survey, fl_min, fl_max, dmin, dmax, lightcurve):
    bursts = np.zeros((n_sources, 3), dtype=np.float64) # initialise the numpy array
    bursts[:,1] = np.absolute(np.power(10, np.random.uniform(np.log10(dmin), np.log10(dmax), n_sources))) # random number for duration
    bursts[:,2] = np.absolute(np.power(10, np.random.uniform(np.log10(fl_min), np.log10(fl_max), n_sources))) # random number for flux
    return bursts
    
def generate_start(bursts, potential_start, potential_end, n_sources):
    bursts[:,0] = np.random.uniform(potential_start, potential_end, n_sources)
    bursts = bursts[bursts[:,0].argsort()] # Order on critical times

    return bursts
    

    

def detect_bursts(obs, flux_err, det_threshold, extra_threshold, sources, gaussiancutoff, edges, fluxint, file, dump_intermediate, write_source, sourceinfo, pointFOV):
    
    simpointings = SkyCoord(ra=sourceinfo['ra'], dec=sourceinfo['dec'],unit='deg',frame='fk5')
    sourceinfo['identity']
    print("Detecting simulated transients in the observations")
    ## pointingcond is a function that detects whether simulated transients are in our point of view or not ##
    pointingcond = lambda f: simpointings.separation(SkyCoord(ra=p[0], dec=p[1], unit='deg', frame='fk5')).deg < f 
    if edges[0] == 1 and edges[1] == 1: # TWO edges
        edgecondmask = lambda start_obs, end_obs: (sources[:,0] + sources[:,1] > start_obs) & (sources[:,0] < end_obs) #
    elif edges[0] == 0 and edges[1] == 1: # Only ending edge, In this case, critical time is the end time
        edgecondmask = lambda start_obs, end_obs: (sources[:,0] > start_obs)
    elif edges[0] == 1 and edges[1] == 0: # Only starting edge 
        edgecondmask = lambda start_obs, end_obs: (sources[:,0] < end_obs)
    elif edges[0] == 0 and edges[1] == 0:
        edgecondmask = lambda start_obs, end_obs: (sources[:,0] == sources[:,0])
    else:
        print("Invalid edges, edges=",edges)
        exit()
    condmask = lambda start_obs, end_obs, p, f: pointingcond(f) & edgecondmask(start_obs,end_obs)
    def detloop(o,p):
       
        f = p[2]
        single_candidate = np.zeros(len(sources),dtype=bool)
        flux_int = np.zeros(len(sources),dtype='f8')
        start_obs, duration, sensitivity = o
        extra_sensitivity = sensitivity * (det_threshold + extra_threshold) / det_threshold
        end_obs = start_obs + duration
        single_candidate = condmask(start_obs, end_obs, p, f)       
        F0_o = sources[single_candidate][:,2]
        error = np.sqrt((abs(random.gauss(F0_o * flux_err, 0.05 * F0_o * flux_err)))**2 + (sensitivity/det_threshold)**2) 
        F0 =random.gauss(F0_o, error) # Simulate some variation in source flux of each source.
        F0[(F0<0)] = F0[(F0<0)]*0
        tau = sources[single_candidate][:,1] # characteristic durations
        tcrit = sources[single_candidate][:,0] # critical times
         # How much time passed in which the transient was on, but not being observed in the survey.
        flux_int[single_candidate] += fluxint(F0, tcrit, tau, end_obs, start_obs)
        candidates = (flux_int > sensitivity) & (single_candidate)  
        extra_candidates = np.array(single_candidate & (flux_int > extra_sensitivity),dtype=bool)
        return candidates, extra_candidates
    file_cleanup = glob.glob("*.npy")
    for f in file_cleanup:
        os.remove(f)
    index =0
    for o,p in tqdm(zip(obs,pointFOV),total=len(obs)):
        candidates, extra_candidates = detloop(o,p)
        np.save("candidates"+str(index),candidates)
        np.save("extra_candidates"+str(index),candidates)
        index+=1
    single_detection = np.zeros((len(obs),len(sources)),dtype=bool)
    extra_detection = np.zeros((len(obs),len(sources)),dtype=bool)
    for i in range(index):
        single_detection[i,:] = np.load('candidates'+str(i)+'.npy')
        extra_detection[i,:] = np.load('extra_candidates'+str(i)+'.npy')
    file_cleanup = glob.glob("*.npy")
    for f in file_cleanup:
        os.remove(f)
    detections = []
    numdet = np.zeros(len(sources),dtype=int)
    numdetextra = np.zeros(len(sources),dtype=int)
    for p in zip(np.unique(pointFOV, axis=0)):
        threetruthpoint = np.sum(pointFOV==p,axis=1)==3
        nobs = np.sum(threetruthpoint)
        numdet += np.sum(single_detection[threetruthpoint,:],axis=0)
        numdetextra += np.sum(extra_detection[threetruthpoint,:],axis=0)
        numdet[numdet==nobs]*=0
        numdetextra[numdetextra==nobs]*=0
    detections = (numdet > 0) & (numdetextra > 0)
    if dump_intermediate:
        with open(file + '_DetTrans', 'w') as f:
            f.write('# Tstart\tDuration\tFlux\n')        ## INITIALISE LIST OF DETECTED TRANSIENTS
            write_source(file + '_DetTrans', sources[detections])
            print("Written Detected Sources")
    return sources[detections], detections

def statistics(fl_min, fl_max, dmin, dmax, det, all_simulated):

    flux_bins = np.geomspace(fl_min, fl_max, num=int(round((np.log10(fl_max)-np.log10(fl_min))/0.05)), endpoint=True)
    dur_ints = np.geomspace(dmin, dmax, num=int(round((np.log10(dmax)-np.log10(dmin))/0.05)), endpoint=True)

    fluxes = np.array([],dtype=np.float32)
    durations = np.array([],dtype=np.float32)

    stats = np.zeros(((len(flux_bins)-1)*(len(dur_ints)-1), 5), dtype=np.float32)
    
    alls = np.array([],dtype=np.uint32)
    dets = np.array([],dtype=np.uint32)
    
    allhistarr, _, _ = np.histogram2d(all_simulated[:,1], all_simulated[:,2],bins = [dur_ints,flux_bins])
    dethistarr, _, _ = np.histogram2d(det[:,1], det[:,2],bins = [dur_ints,flux_bins])
    probabilities = dethistarr/allhistarr
    
    durations = dur_ints[:-1] + dur_ints[1:]/2
    fluxes = flux_bins[:-1] + flux_bins[1:]/2   

    stats[:,0] = np.repeat(durations, len(fluxes))
    stats[:,1] = np.tile(fluxes, len(durations))
    stats[:,2] = probabilities.flatten()
    stats[:,3] = dethistarr.flatten()
    stats[:,4] = allhistarr.flatten()

    return stats


def plots(obs, file, extra_threshold, det_threshold, flux_err, toplot, gaussiancutoff, lclines):
    lightcurve = ''
    toplot[:,0] = np.log10(toplot[:,0])
    toplot[:,1] = np.log10(toplot[:,1])

    gaps = np.array([],dtype=np.float32)
    for i in range(len(obs)-1):
        gaps = np.append(gaps, obs[i+1,0] - obs[i,0] + obs[i,1])
        # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
    min_sens = min(obs[:,2])
    max_sens = max(obs[:,2])
    extra_thresh = max_sens / det_threshold * (extra_threshold + det_threshold)
    sens_last = obs[-1,2]
    sens_maxgap = obs[np.where((gaps[:] == max(gaps)))[0]+1 ,2][0]

    durmax = obs[-1,0] + obs[-1,1] - obs[0,0]
    day1_obs = obs[0,1]
    max_distance = max(gaps)

    dmin=min(toplot[:,0])
    dmax=max(toplot[:,0])
    flmin=min(toplot[:,1])
    flmax=max(toplot[:,1])

    xs = np.arange(dmin, dmax, 1e-3)
    ys = np.arange(flmin, flmax, 1e-3)

    day1_obs_x = np.empty(len(ys))
    day1_obs_x.fill(day1_obs)
    
    sensmin_y = np.empty(len(xs))
    sensmin_y.fill(min_sens)
    
    sensmax_y = np.empty(len(xs))
    sensmax_y.fill(max_sens)

    extra_y = np.empty(len(xs))
    extra_y.fill(extra_thresh)

    X = np.linspace(dmin, dmax, num = 1000)
    Y = np.linspace(flmin, flmax, num = 1000)

    X, Y = np.meshgrid(X, Y)

    Z = interpolate.griddata(toplot[:,0:2], toplot[:,2], (X, Y), method='linear')
    p = figure(title="Probability Contour Plot",tooltips = [("X", "$X"), ("Y", "$Y"), ("value", "@image")], x_axis_type = "log", y_axis_type = "log")
    p.x_range.range_padding = p.y_range.range_padding = 0
    color_mapper = LinearColorMapper(palette="Viridis256",low = 0.0, high = 1.0)
    color_bar = ColorBar(color_mapper=color_mapper, ticker=SingleIntervalTicker(interval = 0.1), label_standoff=12, border_line_color=None, location=(0,0))
    p.image(image=[Z], x=np.amin(10**xs), y=np.amin(10**ys), dw=(np.amax(10**xs)-np.amin(10**xs)), dh=(np.amax(10**ys)-np.amin(10**ys)),palette="Viridis256")
    durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices = lclines(xs, ys, durmax, max_distance, flux_err, obs)   
    p.line(10**xs[durmax_y_indices], durmax_y[durmax_y_indices],  line_width=2, line_color = "red")
    p.line(10**xs[maxdist_y_indices], maxdist_y[maxdist_y_indices],  line_width=2, line_color = "red")
    if durmax_x[0]!=' ':
        p.line(10**durmax_x, 10**ys,   line_width=2, line_color = "red")
    if maxdist_x[0]!=' ':    
        p.line(10**maxdist_x, 10**ys,  line_width=2, line_color = "red")
    if (np.amin(day1_obs_x) > np.amin(10**ys)): p.line(day1_obs_x, 10**ys,  line_width=2, line_color = "red")
    if (sensmin_y[0] > np.amin(10**ys)): p.line(10**xs, sensmin_y,  line_width=2, line_color = "red")
    if (sensmax_y[0] > np.amin(10**ys)): p.line(10**xs, sensmax_y,  line_width=2, line_color = "red")
    if (extra_y[0] > np.amin(10**ys)): p.line(10**xs, extra_y,  line_width=2, line_color = "red")
    p.add_layout(color_bar, 'right')
    p.add_layout(Title(text="Duration (days)", align="center"), "below")
    p.add_layout(Title(text="Transient Flux Density (Jy)", align="center"), "left")
    p.toolbar.logo = None
    p.toolbar_location = None
    p.toolbar.active_drag = None
    p.toolbar.active_scroll = None
    p.toolbar.active_tap = None

    output_file(file + "_ProbContour.html", title = "Probability Contour")
    export_png(p, filename=file + "_ProbContour.png")
    show(p)
    #f = open( 'durmax_y.log', 'w' )
    #for element in durmax_y:
    #    f.write(str(element)+'\n')
    #f.close()
    #f = open( 'maxdist_y.log', 'w' )
    #for element in maxdist_y:
    #    f.write(str(element)+'\n')
    #f.close()
    #f = open( 'xs.log', 'w' )
    #for element in xs:
    #    f.write(str(element)+',')
    #f.close()


def gausscdf(x, t):
    # in the future eliminate the hard-coded 10 here. Replace with calculation based on contents of observations
    return ((x*np.sqrt(np.pi))/(5.0*np.sqrt(2)))*norm.cdf(t, loc = x/2.0, scale = (x/(2*5.0)))

def halfgauss1cdf(x, t):
    # same above with the 5 hint: (2tau/ 10)
    return ((x*np.sqrt(np.pi)*np.sqrt(2))/(5.0))*norm.cdf(t, loc = 0, scale = x/5.0)

def choppedgausscdf(x, t, gaussiancutoff):
    return ((x*np.sqrt(np.pi))/(gaussiancutoff*np.sqrt(2)))*norm.cdf(t, loc = x/2.0, scale = (x/(2*gaussiancutoff)))

