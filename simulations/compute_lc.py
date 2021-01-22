import configparser
import random
import time
import datetime
import numpy as np
import os
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
    else:
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
        obs = observations[observations[:,0].argsort()] # sorts observations by day
        FOV = np.array([1.5 for l in observations])

    return obs, pointing, FOV 

def generate_pointings(n_sources, uniquepoint, FOV):
#Ov erlapping fields
# esp stuff like VLITE
# non-overlaping doesn'at matter so much 
    maxFOV = 1.5
    rng = np.random.default_rng()
    minra = min(uniquepoint[:,0])
    mindec = min(uniquepoint[:,1])
    maxra = max(uniquepoint[:,0])
    maxdec = max(uniquepoint[:,1])
    uniquepointap = SkyCoord(ra=uniquepoint[:,0], dec=uniquepoint[:,1], unit='deg',frame='fk5') 
    # print(uniquepointap)
    # exit()
    matchind = []
    # while len(matchind)<n_sources:
    ra_random = lambda n: (rng.random(int(n))*(min(maxra + maxFOV,360) - max(minra - maxFOV,0)) + max(minra - maxFOV,0)) * u.degree
    dec_random = lambda n: (rng.random(int(n))*(min(maxdec + maxFOV,90)  - max(mindec - maxFOV, -90)) + max(mindec - maxFOV, -90)) * u.degree
    rollforskycoord = lambda n: SkyCoord(ra=ra_random(n), dec=dec_random(n), frame='fk5')
    c = rollforskycoord(n_sources)
    while True:
        accept = np.where(np.array([unique.separation(c).deg for unique in uniquepointap])<1.5)[1]
        mask = np.zeros(c.shape, dtype='bool')
        mask[accept] = True
        imask = np.logical_not(mask)
        print("Rejected sources: ",len(c[imask]))
        if len(c[imask])==0:
            break
        c.data.lon[imask] = ra_random(len(c[imask]))
        c.data.lat[imask] = dec_random(len(c[imask]))
        c.cache.clear()
    
    return c

def generate_sources(n_sources, start_survey, end_survey, fl_min, fl_max, dmin, dmax, lightcurve):
    bursts = np.zeros((n_sources, 3), dtype=np.float64) # initialise the numpy array
    bursts[:,1] = np.absolute(np.power(10, np.random.uniform(np.log10(dmin), np.log10(dmax), n_sources))) # random number for duration
    bursts[:,2] = np.absolute(np.power(10, np.random.uniform(np.log10(fl_min), np.log10(fl_max), n_sources))) # random number for flux
    return bursts
    
def generate_start(bursts, potential_start, potential_end, n_sources):
    bursts[:,0] = np.random.uniform(potential_start, potential_end, n_sources)
    bursts = bursts[bursts[:,0].argsort()] # Order on critical times

    return bursts
    

    

def detect_bursts(obs, flux_err, det_threshold, extra_threshold, sources, gaussiancutoff, edges, fluxint, file, dump_intermediate, write_source, pointing, simpointings):
    lightcurve = ''
    single_detection = np.array([],dtype=np.uint32)
    extra_detection = np.array([],dtype=np.uint32)
    print("Detecting simulated transients in the observations")
    for i in tqdm(range(len(obs))):
        o = obs[i]
        p = pointing[i]
        FOV = 1.5
        start_obs, duration, sensitivity = o
        extra_sensitivity = sensitivity * (det_threshold + extra_threshold) / det_threshold
        end_obs = start_obs + duration
        
        pointingcond = np.where(simpointings.separation(SkyCoord(ra=p[0],dec=p[1], unit='deg',frame='fk5')).deg < FOV)[0]
        pointmask = np.zeros(sources.shape, dtype='bool')
        pointmask[pointingcond] = True

        # The following handles edge cases
        if edges[0] == 1 and edges[1] == 1: # TWO edges
            edgecond = np.where((sources[:,0] + sources[:,1] > start_obs) & (sources[:,0] < end_obs))[0]
        elif edges[0] == 0 and edges[1] == 1: # Only ending edge, In this case, critical time is the end time
            edgecond = np.where(sources[:,0] > start_obs)[0]
        elif edges[0] == 1 and edges[1] == 0: # Only starting edge 
            edgecond = np.where(sources[:,0] < end_obs)[0]
        elif edges[0] == 0 and edges[1] == 0:
            edgecond = np.where(sources[:,0] == sources[:,0])[0]
        else:
            print("Invalid edges, edges=",edges)
            exit()
        
        edgemask = np.zeros(sources.shape, dtype='bool')
        edgemask[edgecond] = True
        
        single_candidate = np.where(pointmask & edgemask)[0]
        
        # filter on integrated flux
        F0_o = sources[single_candidate][:,2]
        error = np.sqrt((abs(random.gauss(F0_o * flux_err, 0.05 * F0_o * flux_err)))**2 + (sensitivity/det_threshold)**2) 
        # F0 =random.gauss(F0_o, error) # Simulate some variation in source flux of each source.
        F0 = F0_o
        F0[(F0<0)] = F0[(F0<0)]*0
        for i in range(len(F0)):
            if F0[i]<0:
                print(F0[i],F0_o[i], error[i])
        # F0 = F0_o
        tau = sources[single_candidate][:,1] # characteristic durations
        tcrit = sources[single_candidate][:,0] # critical times
         # How much time passed in which the transient was on, but not being observed in the survey.
        
        flux_int = fluxint(F0, tcrit, tau, end_obs, start_obs)
        # print(len(single_candidate[flux_int > sensitivity]))
        # print(np.sort(single_candidate))
        candidates = single_candidate[(flux_int > sensitivity)]
        # print(np.sort(flux_int[(flux_int > sensitivity)]/sensitivity))
        # for c in candidates:
            # if sources[c,2] < sensitivity:
                # print(sources[c,2], c)
        extra_candidates = np.array(single_candidate[flux_int > extra_sensitivity])

        single_detection = np.append(single_detection, candidates)
        extra_detection = np.append(extra_detection, extra_candidates)
        
    dets = unique_count(single_detection)
    detections = dets[0][np.where(dets[1] < len(obs))[0]]
    detections = detections[(np.in1d(detections, extra_detection))] # in1d is deprecated consider upgrading to isin
    # for s in sources[detections]:
       # if s[2] < sensitivity:
           # print(s[2], sensitivity, extra_sensitivity)
    
    if dump_intermediate:
        with open(file + '_DetTrans', 'w') as f:
            f.write('# Tstart\tDuration\tFlux\n')        ## INITIALISE LIST OF DETECTED TRANSIENTS
            write_source(file + '_DetTrans', sources[detections])
            print("Written Detected Sources")
    return sources[detections]

def unique_count(a):
    unique, inverse = np.unique(a, return_inverse=True)
    count = np.zeros(len(unique), np.int)
    np.add.at(count, inverse, 1)
    return [unique, count]


def statistics(fl_min, fl_max, dmin, dmax, det, all_simulated):

    flux_ints = np.geomspace(fl_min, fl_max, num=int(round((np.log10(fl_max)-np.log10(fl_min))/0.05)), endpoint=True)
    dur_ints = np.geomspace(dmin, dmax, num=int(round((np.log10(dmax)-np.log10(dmin))/0.05)), endpoint=True)

    fluxes = np.array([],dtype=np.float32)
    durations = np.array([],dtype=np.float32)
#    probabilities = np.array([],dtype=np.float32)

#    stats = np.zeros(((len(flux_ints)-1)*(len(dur_ints)-1), 3), dtype=np.float32)
    stats = np.zeros(((len(flux_ints)-1)*(len(dur_ints)-1), 5), dtype=np.float32)
    
    alls = np.array([],dtype=np.uint32)
    dets = np.array([],dtype=np.uint32)
    
    doit_f = True
    for i in range(len(dur_ints) - 1):

        all_dur = np.array([],dtype=np.uint32)
        det_dur = np.array([],dtype=np.uint32)

        all_dur = np.append(all_dur, np.where((all_simulated[:,1] >= dur_ints[i]) & (all_simulated[:,1] < dur_ints[i+1]))[0])
        det_dur = np.append(det_dur, np.where((det[:,1] >= dur_ints[i]) & (det[:,1] < dur_ints[i+1]))[0])
        durations = np.append(durations, dur_ints[i] + dur_ints[i+1]/2.)

        for m in range(len(flux_ints) - 1):
            dets = np.append(dets, float(len(np.where((det[det_dur,2] >= flux_ints[m]) & (det[det_dur,2]< flux_ints[m+1]))[0])))
            alls = np.append(alls, float(len(np.where((all_simulated[all_dur,2] >= flux_ints[m]) & (all_simulated[all_dur,2] < flux_ints[m+1]))[0])))

            if doit_f == True:
                fluxes = np.append(fluxes,  (flux_ints[m] + flux_ints[m+1])/2.)

        doit_f = False

    probabilities = np.zeros(len(alls),dtype=np.float32)
    toprint_alls = alls[np.where((alls[:] != 0.))]
    toprint_dets = dets[np.where((alls[:] != 0.))]

    probabilities[np.where((alls[:] != 0.))] = dets[np.where((alls[:] != 0.))] / alls[np.where((alls[:] != 0.))]
    
    stats[:,0] = np.repeat(durations, len(fluxes))
    stats[:,1] = np.tile(fluxes, len(durations))
    stats[:,2] = probabilities
    stats[:,3] = dets
    stats[:,4] = alls



#    write_source(file + '_Stat', stats)
    

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

