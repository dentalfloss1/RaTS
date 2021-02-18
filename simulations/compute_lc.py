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
from bitarray import bitarray

def observing_strategy(obs_setup, det_threshold, nobs, obssens, obssig, obsinterval, obsdurations):
    """Parse observation file or set up trial mode. Return array of observation info and a regions observed"""
    
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
        
        observations = np.array(observations,dtype=np.float32)
        # pointing = np.array([SkyCoord(ra=275.0913169*u.degree, dec=7.185135679*u.degree, frame='icrs') for l in observations])
        pointing = np.array([np.array([275.0913169,7.185135679]) for l in observations])
        # pointing[0:4]-=1
        pointing = pointing[observations[:,0].argsort()]
        obs = observations[observations[:,0].argsort()] # sorts observations by date
        FOV = np.array([1.5 for l in observations]) # make FOV for all observations whatever specified here, 1.5 degrees for example
        pointing[0:6,0]+=2 # Make an offset between some pointings
        pointFOV = np.zeros((len(observations),3))
        pointFOV[:,0:2] += pointing
        pointFOV[:,2] += FOV
        uniquepoint = np.unique(pointFOV,axis=0)
        uniquesky = SkyCoord(ra=uniquepoint[:,0],dec=uniquepoint[:,1], unit='deg', frame='fk5')
        ### Next two lines set up variables for defining region properties. Region array is of maximum theoretical length assuming no more than double overlapping fovs (triple or more never get computed ##
        numrgns = len(uniquepoint) 
        regions = np.zeros(np.uint32(numrgns + binom(numrgns,2)), dtype={'names': ('ra', 'dec','identity', 'area', 'timespan'), 'formats': ('f8','f8','U32','f8', 'f8')})
        regions['timespan'][0] = observations[5,0] + observations[5,1] - observations[0,0]
        regions['timespan'][1] = observations[-1,0] + observations[-1,1] - observations[5,0] - observations[5,1]
        regions['timespan'][2] = observations[-1,0] + observations[-1,1] - observations[0,0]
        for i in range(numrgns): # Label the individual pointings. These regions are for example region 1 NOT 2 and 2 NOT 1 
            regions['identity'][i] = str(i)
            regions['ra'][i] = uniquepoint[i,0]
            regions['dec'][i] = uniquepoint[i,1]
            regions['area'][i] = (4*np.pi*np.sin(uniquepoint[i,2]*np.pi/180/2)**2)*(180/np.pi)**2 # Assumes single circular regions, for multiple pointings or other shapes this needs altering
            leftoff = i + 1
        for i in range(len(uniquesky)-1): # Label intersections: For example: 1 AND 2
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
                    regions['area'][leftoff] = (cutchord1 + cutchord2)*(180/np.pi)**2
                    leftoff+=1
    return obs, pointFOV, regions[regions['identity'] != '']
    

def generate_pointings(n_sources, pointFOV, regions):
    """Simulate pointings for each simulated source. Use a monte-carlo like method to roll the dice to determine position. Return a numpy array and bitarray"""
    
    uniquepointFOV = np.unique(pointFOV, axis=0)
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
    btarr = bitarray(n_sources*len(uniquepointFOV))
    btarr.setall(False)
    for i in range(len(uniqueskycoord)):
        seps = c.separation(uniqueskycoord[i]).deg <= (uniquepointFOV[i,2])
        print(seps)
    ra_randfield = lambda n, ptfov: (rng.random(int(n))*(min(ptfov[0] + ptfov[2],360) - max(ptfov[0] - ptfov[2],0)) + max(ptfov[0] - ptfov[2],0)) * u.degree
    dec_randfield = lambda n, ptfov: (rng.random(int(n))*(min(ptfov[1] + ptfov[2],90)  - max(ptfov[1] - ptfov[2], -90)) + max(ptfov[1] - ptfov[2], -90)) * u.degree    
    rollforskycoordfield = lambda n, ptfov: SkyCoord(ra=ra_randfield(n, ptfov), dec=dec_randfield(n, ptfov), frame='fk5')
    totalFOV = np.sum(regions['area'][0:len(uniquepointFOV)]) - np.sum(regions['area'][len(uniquepointFOV):len(regions)])
    print(totalFOV)
    print(regions['area'])
    for i in range(len(uniqueskycoord)):
        targetnum = int(n_sources*regions['timespan'][i]*regions['area'][i]/regions['timespan'][2]/regions['area'][i])# - int(overlappingnum)
        print(targetnum)
        sc = rollforskycoordfield(targetnum, uniquepointFOV[i])
        accept = np.zeros(sc.shape,dtype=bool)
        reject = np.logical_not(accept)
        cond = np.zeros(len(sc),dtype=bool)
        while True:
            cond[reject] +=(sc[reject].separation(uniqueskycoord[i]).deg < (uniquepointFOV[i,2])) 
            if np.sum(~cond)==0:
                break
            reject = ~cond
            accept = ~reject
            sc.data.lon[reject] = ra_randfield(len(sc[reject]),uniquepointFOV[i])
            sc.data.lat[reject] = dec_randfield(len(sc[reject]),uniquepointFOV[i])
            sc.cache.clear()
        print(i*n_sources + btarr[0:i*n_sources].count(bitarray('1')), (i*n_sources + btarr[0:i*n_sources].count(bitarray('1')) + int(targetnum)))
        btarr[(i*n_sources + btarr[0:i*n_sources].count(bitarray('1'))):(i*n_sources + btarr[0:i*n_sources].count(bitarray('1')) + int(targetnum))] = True
    print(btarr[0:n_sources].count(bitarray('1')))
    print(btarr[n_sources:2*n_sources].count(bitarray('1')))
    return btarr

def generate_sources(n_sources, start_survey, end_survey, fl_min, fl_max, dmin, dmax, lightcurve):
    """Generate characteristic fluxes and characteristic durations for simulated sources. Return as numpy array"""
    
    bursts = np.zeros((n_sources, 3), dtype=np.float32) # initialise the numpy array
    bursts[:,1] = np.absolute(np.power(10, np.random.uniform(np.log10(dmin), np.log10(dmax), n_sources))) # random number for duration
    bursts[:,2] = np.absolute(np.power(10, np.random.uniform(np.log10(fl_min), np.log10(fl_max), n_sources))) # random number for flux
    return bursts
    
def generate_start(bursts, potential_start, potential_end, n_sources):
    """Generate characteristic times to go with durations and fluxes. Return modified numpy array"""
    
    bursts[:,0] = np.random.uniform(potential_start, potential_end, n_sources)
    bursts = bursts[bursts[:,0].argsort()] # Order on critical times

    return bursts
    

def detect_bursts(obs, flux_err,  det_threshold, extra_threshold, sources, gaussiancutoff, edges, fluxint, file, dump_intermediate, write_source, pointFOV, ptbtarr):
    """Detect simulated sources by using a series of conditionals along with the integrated flux calculation. Returns detected sources and the boolean array to get source indices"""
    
    start = datetime.datetime.now()
    uniquepointFOV = np.unique(pointFOV,axis=0)
    # This bitarray is defined here because it tracks overall detections through all observations
    detbtarr = bitarray(len(sources))
    detbtarr.setall(False)
    for p in range(len(uniquepointFOV)):
        
        if edges[0] == 1 and edges[1] == 1: # TWO edges
            edgecondmask = lambda start_obs, end_obs: (sources[:,0] + sources[:,1] > start_obs) & (sources[:,0] < end_obs) #
        elif edges[0] == 0 and edges[1] == 1: # Only ending edge, In this case, critical time is the end time
            edgecondmask = lambda start_obs, end_obs: (sources[:,0] > start_obs)
        elif edges[0] == 1 and edges[1] == 0: # Only starting edge 
            edgecondmask = lambda start_obs, end_obs: (sources[:,0] < end_obs)
        elif edges[0] == 0 and edges[1] == 0:
            edgecondmask = lambda start_obs, end_obs: (sources[:,0] == sources[:,0])
        candbitarr = bitarray(len(sources)*len(obs))
        for i in tqdm(range(len(obs))): # bitarray stores whether or not sources fall within observations
            candbitarr[i*len(sources):(i+1)*len(sources)] = bitarray(list(edgecondmask(obs[i,0],obs[i,0] + obs[i,1]))) & ptbtarr[p*len(sources):(p+1)*len(sources)]
     
        end = datetime.datetime.now()
        
        print('conditionals elapsed: ',end-start)
        start = datetime.datetime.now()
        candidates = bitarray(len(sources)) # True if source meets candidate criteria including det_threshold
        candidates.setall(False)
        extra_candidates = bitarray(len(sources)) # True if source meets candidate criteria plus extra threshold
        extra_candidates.setall(False)
        detallbtarr = bitarray(len(sources))
        detallbtarr.setall(True) # Bitarray that determines if source is detected in every observation. If it is, we set it to "not detected" since it is constant.
        for i in tqdm(range(len(obs))):
            flux_int = np.zeros((len(sources)),dtype=np.float32)
            candind = np.array(candbitarr[i*len(sources):(i+1)*len(sources)].search(bitarray([True]))) # Turn the candbitarr into indices. Clunky, but it's the best way to do it I think.
            if candind.size == 0: # No candidates!
                detallbtarr.setall(False) # Otherwise will reject all previous candidates
                break
            F0_o = sources[candind][:,2]
            tcrit = sources[candind][:,0]
            tau = sources[candind][:,1]
            end_obs = obs[i,0]+obs[i,1]
            start_obs = obs[i,0]
            error = np.sqrt((abs(random.gauss(F0_o * flux_err, 0.05 * F0_o * flux_err)))**2 + (obs[i,2]/det_threshold)**2) 
            F0 =random.gauss(F0_o, error)
            F0[(F0<0)] = F0[(F0<0)]*0
            flux_int[candind] += fluxint(F0, tcrit, tau, end_obs, start_obs) # uses whatever class of lightcurve supplied: tophat, ered, etc      
############################# CHECK SENSITIVTY CALC ###############################
            sensitivity = obs[i,2]
            candidates |= bitarray(list(flux_int > obs[i,2])) # Do a bitwise or to determine if candidates meet flux criteria. Sources rejected by edge criteria above are at zero flux anyway
            extra_sensitivity =  obs[i,2] * (det_threshold + extra_threshold) / det_threshold
            extra_candidates |= bitarray(list(flux_int > extra_sensitivity))
            candidates &= extra_candidates # Weird, I know. We really just use the extra detection criteria, but this is a holdover
            detallbtarr &= candidates # End result is true if all are true. In other words, detected in every obs.
        candidates &= ~detallbtarr # We do a bitwise not and the result is that we only keep candidates that were not detected in every obs
        detbtarr |= candidates # Have some persistence between pointings in the detections
        end = datetime.datetime.now()
        print('indexing elapsed: ',end-start)
    detections = np.zeros(len(sources), dtype=bool)
    detections[detbtarr.search(bitarray([True]))] = True # Again, this is how we turn a bitarray into indices
    return sources[detections], detections

def statistics(fl_min, fl_max, dmin, dmax, det, all_simulated):
    """Calculate probabilities based on detections vs simulated, return a numpy array"""
    
    flux_bins = np.geomspace(fl_min, fl_max, num=int(round((np.log10(fl_max)-np.log10(fl_min))/0.05)), endpoint=True)
    dur_ints = np.geomspace(dmin, dmax, num=int(round((np.log10(dmax)-np.log10(dmin))/0.05)), endpoint=True)

    fluxes = np.array([],dtype=np.float32)
    durations = np.array([],dtype=np.float32)

    stats = np.zeros(((len(flux_bins)-1)*(len(dur_ints)-1), 5), dtype=np.float32)
    
    alls = np.array([],dtype=np.uint32)
    dets = np.array([],dtype=np.uint32)
    
    allhistarr, _, _ = np.histogram2d(all_simulated[:,1], all_simulated[:,2],bins = [dur_ints,flux_bins])
    dethistarr, _, _ = np.histogram2d(det[:,1], det[:,2],bins = [dur_ints,flux_bins])
    
    print(np.sort(dethistarr))
    try: # If there were not many sources in the field we may get a 0/0, therefore we set the 0/0 to equal 0. This just appears as a "hole" in the surface density map
        probabilities = dethistarr/allhistarr
    except RuntimeWarning:
        probabilities = np.zeros(dethistarr.shape)
        findzero = np.argwhere(allhistarr==0)
        allhistarr[findzero]+=1
        probabilities = dethistarr/allhistarr
        
    durations = dur_ints[:-1] + dur_ints[1:]/2
    fluxes = flux_bins[:-1] + flux_bins[1:]/2   

    # output to static HTML file
    output_file("line.html")

    poo = figure(plot_width=400, plot_height=400, x_axis_type = "log", y_axis_type = "log")

    # add a circle renderer with a size, color, and alpha
    poo.circle(dur_ints[:-1], dethistarr[:,50], size=10, color="navy", alpha=0.5)

    # show the results
    show(poo)

    stats[:,0] = np.repeat(durations, len(fluxes))
    stats[:,1] = np.tile(fluxes, len(durations))
    stats[:,2] = probabilities.flatten()
    stats[:,3] = dethistarr.flatten()
    stats[:,4] = allhistarr.flatten()

    return stats


def plots(obs, file, extra_threshold, det_threshold, flux_err, toplot, gaussiancutoff, lclines):
    """Using stats, observations, and what not, generate a plot using bokeh"""
    
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

