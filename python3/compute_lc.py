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
from scipy.stats import norm
from scipy.special import erf


def observing_strategy(obs_setup, det_threshold, nobs, obssens, obssig, obsinterval, obsdurations):
    observations = []
    try:
        with open(obs_setup, 'r') as f:
            for line in f:
                if len(line) == 0 or line[0] == '#':
                    continue
                cols = line.split('\t')
                tstart = time.mktime(datetime.datetime.strptime(cols[0], "%Y-%m-%dT%H:%M:%S.%f").timetuple())/(3600*24)    #in days
                tdur = float(cols[1])    # in days
                sens = float(cols[2]) * det_threshold    # in Jy
                observations.append([tstart, tdur, sens])
    except TypeError:
        for i in range(nobs):
            tstart = (time.mktime(datetime.datetime.strptime("2019-08-08T12:50:05.0", "%Y-%m-%dT%H:%M:%S.%f").timetuple()) + 60*60*24*obsinterval*i)/(3600*24)    #in days at an interval of 7 days
            tdur = obsdurations
            sens = random.gauss(obssens, obssig) * det_threshold # in Jy. Choice of Gaussian and its characteristics were arbitrary.
            observations.append([tstart,tdur,sens])
            
    observations = np.array(observations,dtype=np.float64)
    obs = observations[observations[:,0].argsort()] # sorts observations by day

    return obs

def generate_sources(n_sources, start_survey, end_survey, fl_min, fl_max, dmin, dmax, lightcurve):
    bursts = np.zeros((n_sources, 3), dtype=np.float64) # initialise the numpy array
    bursts[:,1] = np.absolute(np.power(10, np.random.uniform(np.log10(dmin), np.log10(dmax), n_sources))) # random number for duration
    bursts[:,2] = np.absolute(np.power(10, np.random.uniform(np.log10(fl_min), np.log10(fl_max), n_sources))) # random number for flux
    return bursts
    
def generate_start(bursts, potential_start, potential_end, n_sources):
    bursts[:,0] = np.random.uniform(potential_start, potential_end, n_sources)
    bursts = bursts[bursts[:,0].argsort()] # Order on critical times
    # if(lightcurve == "gaussian") or (lightcurve == "choppedgaussian"):  
        # bursts[:,0] = np.random.uniform(start_survey - bursts[:,1], end_survey + bursts[:,1], n_sources)
    # elif((lightcurve == "tophat") or (lightcurve == "fred") or (lightcurve == "halfgaussian1") or (lightcurve == "parabolic")):
        # bursts[:,0] = np.random.uniform(start_survey - bursts[:,1], end_survey, n_sources) # randomize the start times based partly on the durations above
    # elif(lightcurve == "wilma"):
        # bursts[:,0] = np.random.uniform(start_survey - bursts[:,1], end_survey + bursts[:,1], n_sources) # randomize the start times based partly on the durations above
    # elif(lightcurve == "ered"): 
        # bursts[:,0] = np.random.uniform(start_survey - bursts[:,1], end_survey + bursts[:,1], n_sources) # randomize the start times based partly on the durations above
    return bursts
    

    

def detect_bursts(obs, flux_err, det_threshold, extra_threshold, sources, gaussiancutoff, edges, fluxint):
    lightcurve = ''
    single_detection = np.array([],dtype=np.uint32)
    extra_detection = np.array([],dtype=np.uint32)

    for o in obs:
        start_obs, duration, sensitivity = o
        extra_sensitivity = sensitivity * (det_threshold + extra_threshold) / det_threshold
        end_obs = start_obs + duration
        
        # The following handles edge cases
        if edges[0] == 1 and edges[1] == 1: # TWO edges
            single_candidate = np.where((sources[:,0] + sources[:,1] > start_obs) & (sources[:,0] < end_obs))[0]
        elif edges[0] == 0 and edges[1] == 1: # Only ending edge, In this case, critical time is the end time
            single_candidate = np.where(sources[:,0] > start_obs)[0]
        elif edges[0] == 1 and edges[1] == 0: # Only starting edge 
            single_candidate = np.where(sources[:,0] < end_obs)[0]
        elif edges[0] == 0 and edges[1] == 0:
            single_candidate = np.where(sources[:,0] == sources[:,0])[0]
        else:
            print("Invalid edges, edges=",edges)
        
        # if ((lightcurve == 'fred') or (lightcurve == 'halfgaussian1')):
            # Transients start before observation ends
            # single_candidate = np.where(sources[:,0] < end_obs)[0]
        # elif lightcurve == 'wilma':
            # Transients do not end before the observation starts 
            # single_candidate = np.where(sources[:,0] + sources[:,1] > start_obs)[0]
        # elif lightcurve == 'tophat'  or lightcurve == 'parabolic':
            # Transients end after observation starts and start before observation ends
            # single_candidate = np.where((sources[:,0] + sources[:,1] > start_obs) & (sources[:,0] < end_obs))[0] # Take all of the locations that qualify as a candidate. The zero index is a wierd python workaround
        # elif lightcurve == 'gaussian' or lightcurve == 'choppedgaussian':
            # single_candidate = np.where(sources[:,0] == sources[:,0])[0]
        # elif lightcurve == 'ered':
            # single_candidate = np.where((sources[:,0] + sources[:,1] > start_obs) & (sources[:,0] < end_obs) )[0]


        # filter on integrated flux
        F0_o = sources[single_candidate][:,2]
        error = np.sqrt((abs(random.gauss(F0_o * flux_err, 0.05 * F0_o * flux_err)))**2 + (sensitivity/det_threshold)**2) 
        F0 =random.gauss(F0_o, error) # Simulate some variation in source flux of each source.
        # F0 = F0_o
        tau = sources[single_candidate][:,1] # characteristic durations
        tcrit = sources[single_candidate][:,0] # critical times
         # How much time passed in which the transient was on, but not being observed in the survey.
        
        flux_int = fluxint(F0, tcrit, tau, end_obs, start_obs)
        
        if lightcurve == 'fred':
            tend = end_obs - t_burst # Burst is never really "Off" so, how long is it on? 
            flux_int = np.multiply(F0, np.multiply(tau, np.divide(np.exp(-np.divide(tstart,tau)) - np.exp(-np.divide(tend,tau)), (end_obs-start_obs))))
        
        elif lightcurve == 'parabolic':
            tend = np.minimum(t_burst + tau, end_obs) - t_burst
            deltat = tend-tstart
            f1 = F0*deltat
            
            tpk = tau/2.0
            deltat2 = tend - tpk
            deltat1 = tstart - tpk
            deltat2c = np.power(deltat2,3)
            deltat1c = np.power(deltat1,3)
            deltadelta = deltat2c - deltat1c
            fdelt = F0*deltadelta
            fdelt1 = fdelt/3.0
            durhalf = tau/2.0
            durc = np.power(durhalf,3)
            fdelt2 = fdelt1/durc
            
            flux_int_pre = f1 - fdelt2
            flux_int = flux_int_pre/(end_obs-start_obs)
            
            

        elif lightcurve == 'wilma':
            t1 = np.minimum(end_obs,t_burst + tau) - (t_burst + tau)
            t2 = start_obs - (t_burst + tau)
            t3 = np.divide(t1,tau)
            e1 = np.exp(t3)
            try: 
                e2 = np.exp(np.divide(t2,tau))
            except: 
                print(np.sort(np.divide(t2,tau), axis=-1, kind='quicksort', order=None))
                print(np.sort(t2, axis=-1, kind='quicksort', order=None))
                print(np.sort(tau, axis=-1, kind='quicksort', order=None))
                input("Press enter to continue")
                e2 = e1
            e3 = np.divide(e1,(end_obs-start_obs))
            e4 = np.divide(e2,(end_obs-start_obs))
            etot = e3-e4
            f1 = np.multiply(F0,etot)
            f2 = np.multiply(f1,tau)
            flux_int = f2
            
            
            # try: 
                # flux_int = np.multiply(F0, np.multiply(tau, np.divide(np.exp(-np.divide(np.minimum(end_obs_arr,t_burst)-t_burst,tau)-np.exp(-np.divide(start_obs - t_burst,tau))), (end_obs-start_obs))))
            # except: 
                # flux_int = 0
                # print(np.multiply(F0, np.multiply(tau, np.divide(np.exp(-np.divide(np.minimum(end_obs_arr,t_burst)-t_burst,tau)-np.exp(-np.divide(start_obs - t_burst,tau))), (end_obs-start_obs)))))
        elif lightcurve == 'ered':
            tend = end_obs - t_burst # Burst is never really "Off" so, how long is it on? 
            flux_int_fred = np.multiply(F0, np.multiply(tau, np.divide(np.exp(-np.divide(tstart,tau)) - np.exp(-np.divide(tend,tau)), (end_obs-start_obs))))
            # flux_int_fred = 0
            
            t1 = np.minimum(end_obs,t_burst + tau) - (t_burst + tau)
            t2 = start_obs - (t_burst + tau)
            t3 = np.divide(t1,tau)
            e1 = np.exp(t3)
            try: 
                e2 = np.exp(np.divide(t2,tau))
            except: 
                print(np.sort(np.divide(t2,tau), axis=-1, kind='quicksort', order=None))
                print(np.sort(t2, axis=-1, kind='quicksort', order=None))
                print(np.sort(tau, axis=-1, kind='quicksort', order=None))
                input("Press enter to continue")
                e2 = e1
            e3 = np.divide(e1,(end_obs-start_obs))
            e4 = np.divide(e2,(end_obs-start_obs))
            etot = e3-e4
            f1 = np.multiply(F0,etot)
            f2 = np.multiply(f1,tau)
            flux_int_wilma = f2
            flux_int = flux_int_fred + flux_int_wilma
        
        elif lightcurve == 'choppedgaussian':
            sensrat = sensitivity/F0
            natlog = np.log(sensrat)
            doublenatlog = -2.0*natlog
            doublenatlog[doublenatlog<=0] = 0.5**2
            sqnatlog = np.sqrt(doublenatlog)
            xpre = sqnatlog
            xpre[xpre > gaussiancutoff] = gaussiancutoff
            x = xpre
            tend = np.minimum(t_burst + tau, end_obs) - t_burst
            # doublepi = 2.0*np.pi
            sqpi = np.sqrt(np.pi)
            scaletau = tau/(x*np.sqrt(2))
            dimensionless = scaletau/(end_obs - start_obs)
            fpre1 = F0 * dimensionless
            fpre2 = fpre1 * sqpi
            cdf2 = norm.cdf(end_obs , loc = t_burst + (tau/2.0), scale = (0.5*tau)/x)
            cdf1 = norm.cdf(start_obs , loc = t_burst + (tau/2.0), scale = (0.5*tau)/x)
            diffcdf = cdf2 - cdf1 
            flux_int = fpre2*diffcdf
        elif lightcurve == 'gaussian':
            sensrat = sensitivity/F0
            natlog = np.log(sensrat)
            doublenatlog = -2.0*natlog
            doublenatlog[doublenatlog<=0] = 0.5**2
            sqnatlog = np.sqrt(doublenatlog)
            x = sqnatlog
            # print(x)
            tend = np.minimum(t_burst + tau, end_obs) - t_burst
            # doublepi = 2.0*np.pi
            sqpi = np.sqrt(np.pi)
            scaletau = tau/(x*np.sqrt(2))
            dimensionless = scaletau/(end_obs - start_obs)
            fpre1 = F0 * dimensionless
            fpre2 = fpre1 * sqpi
            cdf2 = norm.cdf(end_obs , loc = t_burst + (tau/2.0), scale = (0.5*tau)/x)
            cdf1 = norm.cdf(start_obs , loc = t_burst + (tau/2.0), scale = (0.5*tau)/x)
            diffcdf = cdf2 - cdf1 
            flux_int = fpre2*diffcdf
            # flux_int = np.sqrt(2.0*np.pi)*(tau/(gaussiancutoff*(end_obs-start_obs)))*np.multiply(F0, norm.cdf(end_obs , loc = t_burst + (tau/2.0), scale = tau/gaussiancutoff)-norm.cdf(start_obs , loc = t_burst + (tau/2.0), scale = tau/gaussiancutoff))
            # flux_int = np.multiply(F0,(-1.0/2.0)*erf((3.0*(-2.0*end_obs+2.0*t_burst+tau))/(np.sqrt(2)*tau))+(1.0/2.0)*erf((3.0*(-2.0*start_obs+2.0*t_burst+tau))/(np.sqrt(2)*tau)))
        elif lightcurve == 'halfgaussian1':
            sensrat = sensitivity/F0
            natlog = np.log(sensrat)
            doublenatlog = -2.0*natlog
            doublenatlog[doublenatlog<=0] = 0.5**2
            sqnatlog = np.sqrt(doublenatlog)
            x = sqnatlog
            # print(x)
            tend = np.minimum(t_burst + tau, end_obs) - t_burst
            # doublepi = 2.0*np.pi
            sqpi = np.sqrt(np.pi)
            scaletau = (tau*np.sqrt(2))/(x)
            dimensionless = scaletau/(end_obs - start_obs)
            fpre1 = F0 * dimensionless
            fpre2 = fpre1 * sqpi
            cdf2 = norm.cdf(end_obs , loc = t_burst , scale = tau/x)
            cdf1 = norm.cdf(np.maximum(t_burst, start_obs) , loc = t_burst , scale = tau/x)
            diffcdf = cdf2 - cdf1 
            flux_int = fpre2*diffcdf
            
            # flux_int = np.sqrt(2.0*np.pi)*(tau/(5.0*(end_obs-start_obs)))*np.multiply(F0, norm.cdf(end_obs , loc = t_burst, scale = tau/5.0)-norm.cdf(np.maximum(t_burst, start_obs) , loc = t_burst, scale = tau/5.0))
            # flux_int = np.multiply(F0,(-1.0/2.0)*erf((3.0*(-2.0*end_obs+2.0*t_burst+tau))/(np.sqrt(2)*tau))+(1.0/2.0)*erf((3.0*(-2.0*start_obs+2.0*t_burst+tau))/(np.sqrt(2)*tau)))
        candidates = single_candidate[(flux_int > sensitivity)]
        extra_candidates = np.array(single_candidate[flux_int > extra_sensitivity])

        single_detection = np.append(single_detection, candidates)
        extra_detection = np.append(extra_detection, extra_candidates)
        
    dets = unique_count(single_detection)
    detections = dets[0][np.where(dets[1] < len(obs))[0]]
    detections = detections[(np.in1d(detections, extra_detection))] # in1d is deprecated consider upgrading to isin
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
    if (lightcurve == 'fred'):
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs / np.power(10,x) / (np.exp(-(durmax - day1_obs + np.power(10,x)) /  np.power(10,x)) - np.exp(-((durmax + np.power(10,x)) / np.power(10,x)))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, (((1. + flux_err) * sens_maxgap * day1_obs) /  np.power(10,x))   / (np.exp(-(max_distance / np.power(10,x))) - np.exp(-(max_distance + day1_obs) / np.power(10,x))))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)    

    # elif (lightcurve == 'tophat'):
        # durmax_x = np.empty(len(ys))
        # durmax_x.fill(np.log10(durmax))
        # maxdist_x = np.empty(len(ys))
        # maxdist_x.fill(np.log10(max_distance))
        
    elif (lightcurve == 'parabolic'):
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        maxdist_x = np.empty(len(ys))
        maxdist_x.fill(np.log10(max_distance))
        
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs /  (day1_obs - (1.0/(3*np.power(np.power(10,x)/2.0,2.0)))*(np.power(durmax + np.power(10,x), 3.0) - np.power(durmax - day1_obs + np.power(10.0, x), 3.0))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y = np.append(maxdist_y, (1. + flux_err) * sens_maxgap * day1_obs /  (day1_obs - (1.0/(3*np.power(np.power(10,x)/2.0,2.0)))*(np.power(max_distance + day1_obs, 3.0) - np.power(max_distance, 3.0))))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)    
    elif (lightcurve == 'gaussian'):
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, ((1. + flux_err) * sens_last * day1_obs  ) / (gausscdf(np.power(10,x),durmax + np.power(10,x)) - gausscdf(np.power(10,x), durmax - day1_obs + np.power(10,x))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, ((1. + flux_err) * sens_last * day1_obs ) / (gausscdf(np.power(10,x),max_distance + day1_obs ) - gausscdf(np.power(10,x), max_distance)))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)
    elif (lightcurve == 'choppedgaussian'):
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, ((1. + flux_err) * sens_last * day1_obs  ) / (choppedgausscdf(np.power(10,x),durmax + np.power(10,x), gaussiancutoff) - choppedgausscdf(np.power(10,x), durmax - day1_obs + np.power(10,x), gaussiancutoff)))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, ((1. + flux_err) * sens_last * day1_obs ) / (choppedgausscdf(np.power(10,x),max_distance + day1_obs, gaussiancutoff) - choppedgausscdf(np.power(10,x), max_distance,gaussiancutoff)))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)
#        maxdist_y = np.array([],dtype=np.float64)
#        for x in xs:
#            try:
#                durmax_y = np.append(durmax_y, ((1. + flux_err) * sens_last * day1_obs  ) / (choppedgausscdf(np.power(10,x),durmax + np.power(10,x), gaussiancutoff) - choppedgausscdf(np.power(10,x), durmax - day1_obs + np.power(10,x), gaussiancutoff)))
#            except:
#                durmax_y = np.append(durmax_y, np.inf)
#            try:
#                maxdist_y =  np.append(maxdist_y, ((1. + flux_err) * sens_last * day1_obs ) / (choppedgausscdf(np.power(10,x),max_distance + day1_obs, gaussiancutoff) - choppedgausscdf(np.power(10,x), max_distance, gaussiancutoff)))
#            except:
#                maxdist_y = np.append(maxdist_y, np.inf)
#    
    elif (lightcurve == 'halfgaussian1'):
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, ((1. + flux_err) * sens_last * day1_obs  ) / (halfgauss1cdf(np.power(10,x),durmax + np.power(10,x)) - halfgauss1cdf(np.power(10,x), durmax - day1_obs + np.power(10,x))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, ((1. + flux_err) * sens_last * day1_obs ) / (gausscdf(np.power(10,x),max_distance + day1_obs) - gausscdf(np.power(10,x), max_distance)))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)
    
    elif (lightcurve == 'wilma'):  
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs / np.power(10,x) / (np.exp(-(durmax - day1_obs + np.power(10,x)) /  np.power(10,x)) - np.exp(-((durmax + np.power(10,x)) / np.power(10,x)))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, (((1. + flux_err) * sens_maxgap * day1_obs) /  np.power(10,x))   / (np.exp(-(max_distance / np.power(10,x))) - np.exp(-(max_distance + day1_obs) / np.power(10,x))))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)
                
    elif (lightcurve == 'ered'):  
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs / np.power(10,x) / (2.0 -  np.exp(-(durmax - day1_obs + np.power(10,x)) /  np.power(10,x)) - np.exp(-((durmax + np.power(10,x)) / np.power(10,x)))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, (((1. + flux_err) * sens_maxgap * day1_obs) /  np.power(10,x))   / (2.0 - np.exp(-(max_distance / np.power(10,x))) - np.exp(-(max_distance + day1_obs) / np.power(10,x))))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)
                
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
    if True: 
        if lightcurve == 'fred':
            durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) & (durmax_y > np.amin(10**ys)))[0]
            maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
            p.line(10**xs[durmax_y_indices], durmax_y[durmax_y_indices],  line_width=2, line_color = "red")
            p.line(10**xs[maxdist_y_indices], maxdist_y[maxdist_y_indices], line_width=2, line_color = "red")
        elif lightcurve == 'wilma':
            durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) & (durmax_y > np.amin(10**ys)))[0]
            maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
            p.line(10**xs[durmax_y_indices], durmax_y[durmax_y_indices],  line_width=2, line_color = "red")
            p.line(10**xs[maxdist_y_indices], maxdist_y[maxdist_y_indices], line_width=2, line_color = "red")
        elif lightcurve == 'ered':
            durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) & (durmax_y > np.amin(10**ys)))[0]
            maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
            p.line(10**xs[durmax_y_indices], durmax_y[durmax_y_indices],  line_width=2, line_color = "red")
            p.line(10**xs[maxdist_y_indices], maxdist_y[maxdist_y_indices], line_width=2, line_color = "red")
            durmax_x = np.empty(len(ys))
            durmax_x.fill(np.log10(durmax))
            maxdist_x = np.empty(len(ys))
            maxdist_x.fill(np.log10(max_distance))
            p.line(10**durmax_x, 10**ys,   line_width=2, line_color = "red")
            p.line(10**maxdist_x, 10**ys,  line_width=2, line_color = "red")
        # elif lightcurve == 'tophat':
            # p.line(10**durmax_x, 10**ys,   line_width=2, line_color = "red")
            # p.line(10**maxdist_x, 10**ys,  line_width=2, line_color = "red")
        elif (lightcurve == 'gaussian') or (lightcurve == 'halfgaussian1') or (lightcurve == 'choppedgaussian'):
            durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) &  (durmax_y > np.amin(10**ys)))[0]
            maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
            # print(durmax_y)
            # print(maxdist_y)
            p.line(10**xs[durmax_y_indices], durmax_y[durmax_y_indices],  line_width=2, line_color = "red")
            p.line(10**xs[maxdist_y_indices], maxdist_y[maxdist_y_indices],  line_width=2, line_color = "red")
            p.line(10**durmax_x, 10**ys,   line_width=2, line_color = "red")
            # p.line(10**xs, durmax_y,  line_width=2, line_color = "red")
            # p.line(10**xs, maxdist_y,   line_width=2, line_color = "red")
        elif lightcurve == 'parabolic':
            durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) &  (durmax_y > np.amin(10**ys)))[0]
            maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
            p.line(10**xs[durmax_y_indices], durmax_y[durmax_y_indices],  line_width=2, line_color = "red")
            p.line(10**xs[maxdist_y_indices], maxdist_y[maxdist_y_indices],  line_width=2, line_color = "red")
            p.line(10**durmax_x, 10**ys,   line_width=2, line_color = "red")
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

