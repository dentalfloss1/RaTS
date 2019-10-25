#!/usr/local/anaconda2/bin/python

import parset
import random
import time
import datetime
import numpy as np
import os
import argparse

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib import rc
import pylab
rc('text', usetex=True)


def get_configuration():
    """Rerturns a populated configuration"""
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--parset", help="Parset filename")
    argparser.add_argument("--observations", help="Observation filename")    # format date YYYY-MM-DDTHH:MM:SS.mmmmmm \t dur (days) \t sens (Jy)
    argparser.add_argument("--transient_type", help="Transient lightcurve")    # format: 'fred' or 'tophat'
    argparser.add_argument("--keep", action='store_true', help="Keep previous bursts")
    argparser.add_argument("--stat_plot", action='store_true', help="Performs statistics and plots results")
    argparser.add_argument("--just_plot", action='store_true', help="Just plots results")


    return argparser.parse_args()

def write_source(filename, bursts):
    with open(filename, 'a') as f:
        for burst in bursts:
            f.write("{0}\t{1}\t{2}\n".format(burst[0], burst[1], burst[2]))            ### for python 2.6 (krieger)
#            f.write("{}\t{}\t{}\n".format(burst[0], burst[1], burst[2]))            ### for python 2.7 (struis)
        f.flush()

def initialise():
    config = get_configuration()
    obs_setup = config.observations
    lightcurve = config.transient_type

    if lightcurve != 'tophat' and lightcurve != 'fred':
        print 'Type of transient not recognised.\nUse tophat or fred.'
        exit()

    # Transients parameters
    parse = parset.Parset(config.parset)
#    parse = lofar.parameterset.parameterset(config.parset)
    n_sources = parse.getInt("n_sources")
    fl_min = parse.getFloat("fl_min")                    #Flux in same units as the flux in the observations file
    fl_max = parse.getFloat("fl_max")
    flux_err = parse.getFloat("flux_err")                #fractional
    dmin = parse.getFloat("dmin")                        #Duration in same units as the duration in the observations file
    dmax = parse.getFloat("dmax")
    det_threshold = parse.getInt("det_threshold")
    extra_threshold = parse.getInt("extra_threshold")
    file = parse.getString("file")

    file = file + '_' + lightcurve
    obs = observing_strategy(obs_setup, det_threshold)

    if config.just_plot:
        plots(obs, file, extra_threshold, det_threshold, flux_err, lightcurve)
        print 'done'
        exit()        

    if not config.keep:
        with open(file + '_SimTrans', 'w') as f:
            f.write('# Tstart\tDuration\tFlux\n')        ## INITIALISE LIST OF SIMILATED TRANSIENTS
    
        with open(file + '_DetTrans', 'w') as f:
            f.write('# Tstart\tDuration\tFlux\n')        ## INITIALISE LIST OF DETECTED TRANSIENTS
    

    start_time = obs[0][0]
    end_time = obs[-1][0] + obs[-1][1]

    simulated = generate_sources(n_sources, file, start_time, end_time, fl_min, fl_max, dmin, dmax)
    detect_bursts(obs, file, flux_err, det_threshold, extra_threshold, simulated, lightcurve)
    
    if config.stat_plot:
        with open(file + '_Stat', 'w') as f:
            f.write('# Duration\tFlux\tProbability\n')        ## INITIALISE LIST OF STATISTICS
        statistics(file, fl_min, fl_max, dmin, dmax)
        plots(obs, file, extra_threshold, det_threshold, flux_err, lightcurve)

def observing_strategy(obs_setup, det_threshold):
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
        for i in range(46):
            tstart = (time.mktime(datetime.datetime.strptime("2019-08-08T12:50:05.0", "%Y-%m-%dT%H:%M:%S.%f").timetuple()) + 60*60*24*7*i)/(3600*24)    #in days at an interval of 7 days
            tdur = 0.009
            sens = random.gauss(0.0000317, 0.0000046) * det_threshold # in Jy. Choice of Gaussian and its characteristics were arbitrary.
            observations.append([tstart,tdur,sens])
    observations = np.array(observations,dtype=np.float64)
    obs = observations[observations[:,0].argsort()] # sorts observations by day

    return obs

def generate_sources(n_sources, file, start_time, end_time, fl_min, fl_max, dmin, dmax):
    bursts = np.zeros((n_sources, 3), dtype=np.float64)
    bursts[:,1] = np.absolute(np.power(10, np.random.uniform(np.log10(dmin), np.log10(dmax), n_sources)))
    bursts[:,2] = np.absolute(np.power(10, np.random.uniform(np.log10(fl_min), np.log10(fl_max), n_sources)))
    bursts[:,0] = np.random.uniform(start_time - bursts[:,1], end_time, n_sources)
    
    bursts = bursts[bursts[:,0].argsort()]
    write_source(file + '_SimTrans', bursts)
    
    print "Written Simulated Sources"
    return bursts

def detect_bursts(obs, file, flux_err, det_threshold, extra_threshold, sources, lightcurve):
    detection = np.array([],dtype=np.uint32)
    extra_detection = np.array([],dtype=np.uint32)

    for o in obs:
        start_obs, duration, sensitivity = o
        extra_sensitivity = sensitivity * (det_threshold + extra_threshold) / det_threshold
        end_obs = start_obs + duration
        
        if lightcurve == 'fred':
            ## Transients start before observation ends
            candidate = np.where(sources[:,0] < end_obs)[0]
        elif lightcurve == 'tophat':
            ## Transients end after observation starts and start before observation ends
            candidate = np.where((sources[:,0] + sources[:,1] > start_obs) & (sources[:,0] < end_obs))[0]

        # filter on integrated flux
        F0_o = sources[candidate][:,2]
        error = np.sqrt((abs(random.gauss(F0_o * flux_err, 0.05 * F0_o * flux_err)))**2 + (sensitivity/det_threshold)**2)#
        F0 = random.gauss(F0_o, error)

        tau = sources[candidate][:,1]
        t_burst = sources[candidate][:,0]
        tstart = np.maximum(t_burst, start_obs) - t_burst

        if lightcurve == 'fred':
            tend = end_obs - t_burst
            flux_int = np.multiply(F0, np.multiply(tau, np.divide(np.exp(-np.divide(tstart,tau)) - np.exp(-np.divide(tend,tau)), (end_obs-start_obs))))
        
        elif lightcurve == 'tophat':
            tend = np.minimum(t_burst + tau, end_obs) - t_burst
            flux_int = np.multiply(F0, np.divide((tend - tstart), (end_obs-start_obs)))
        
        candidates = candidate[(flux_int > sensitivity)]
        extra_candidates = np.array(candidate[flux_int > extra_sensitivity])

        detection = np.append(detection, candidates)
        extra_detection = np.append(extra_detection, extra_candidates)
        
        
    dets = unique_count(detection)
    detections = dets[0][np.where(dets[1] < len(obs))[0]]
    
    detections = detections[(np.in1d(detections, extra_detection))]
    write_source(file + '_DetTrans', sources[detections])
    
    print "Written Detected Sources"
    return sources[detections]

def unique_count(a):
    unique, inverse = np.unique(a, return_inverse=True)
    count = np.zeros(len(unique), np.int)
    np.add.at(count, inverse, 1)
    return [unique, count]


def statistics(file, fl_min, fl_max, dmin, dmax):
    all = np.loadtxt(file + '_SimTrans')
    det = np.loadtxt(file + '_DetTrans')

    flux_ints = np.geomspace(fl_min, fl_max, num=(np.log10(fl_max)-np.log10(fl_min))/0.05, endpoint=True)
    dur_ints = np.geomspace(dmin, dmax, num=(np.log10(dmax)-np.log10(dmin))/0.05, endpoint=True)

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

        all_dur = np.append(all_dur, np.where((all[:,1] >= dur_ints[i]) & (all[:,1] < dur_ints[i+1]))[0])
        det_dur = np.append(det_dur, np.where((det[:,1] >= dur_ints[i]) & (det[:,1] < dur_ints[i+1]))[0])
        durations = np.append(durations, dur_ints[i] + dur_ints[i+1]/2.)

        for m in range(len(flux_ints) - 1):
            dets = np.append(dets, float(len(np.where((det[det_dur,2] >= flux_ints[m]) & (det[det_dur,2]< flux_ints[m+1]))[0])))
            alls = np.append(alls, float(len(np.where((all[all_dur,2] >= flux_ints[m]) & (all[all_dur,2] < flux_ints[m+1]))[0])))

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
    write_stat(file + '_Stat', stats)

    print "Written Statistics"

    return stats

def write_stat(filename, bursts):
    with open(filename, 'a') as f:
        for burst in bursts:
            f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(burst[0], burst[1], burst[2], burst[3], burst[4]))        # for python 2.6 (krieger)
#            f.write("{}\t{}\t{}\t{}\t{}\n".format(burst[0], burst[1], burst[2], burst[3], burst[4]))        # for python 2.7 (struis)
        f.flush()


def plots(obs, file, extra_threshold, det_threshold, flux_err, lightcurve):
    toplot = np.loadtxt(file + '_Stat')

    toplot[:,0] = np.log10(toplot[:,0])
    toplot[:,1] = np.log10(toplot[:,1])

    gaps = np.array([],dtype=np.float32)
    for i in range(len(obs)-1):
        gaps = np.append(gaps, obs[i+1,0] - obs[i,0] + obs[i,1])

    min_sens = min(obs[:,2])
    max_sens = max(obs[:,2])
    extra_thresh = max_sens / det_threshold * (extra_threshold + det_threshold)
    sens_last = obs[-1,2]
    sens_maxgap = obs[np.where((gaps[:] == max(gaps)))[0]+1 ,2][0]

    durmax = obs[-1,0] + obs[-1,1] - obs[0,0]
    day1_obs = obs[0,1]
    max_distance = max(gaps)

    xlabel = 'Log Transient Duration [days]'
    ylabel = 'Log Transient Flux Density [Jy]'
    plotname = 'probability_contour'

    fig = plt.figure()
    pylab.xlabel(r'{Log Transient Duration [days]', {'color':'k', 'fontsize':16})
    pylab.ylabel(r'{Log Transient Flux Density [Jy]', {'color':'k', 'fontsize':16})
    pylab.xticks(fontsize=18)
    pylab.yticks(fontsize=18)

    dmin=min(toplot[:,0])
    dmax=max(toplot[:,0])
    flmin=min(toplot[:,1])
    flmax=max(toplot[:,1])

    xs = np.arange(dmin, dmax, 0.01)
    ys = np.arange(flmin, flmax, 0.01)


    if lightcurve == 'fred':
        durmax_y = np.log10((1. + flux_err) * sens_last * day1_obs / np.power(10, xs) / (np.exp(-(durmax - day1_obs + np.power(10, xs)) / np.power(10, xs)) - np.exp(-((durmax + np.power(10, xs)) / np.power(10, xs)))))
        maxdist_y = np.log10((1. + flux_err) * sens_maxgap * day1_obs / np.power(10, xs) / (np.exp(-(max_distance / np.power(10, xs))) - np.exp(-(max_distance + day1_obs) / np.power(10, xs))))

        tmp1 = np.zeros((0,),dtype = np.float64)
        tmp2 = np.zeros((0,), dtype = np.float64)
        for i in range(len(xs)):
            tmp1 = np.append(tmp1,(1. + flux_err) * sens_last * day1_obs / xs[i] / (np.exp(-(durmax - day1_obs + xs[i]) /  xs[i]) - np.exp(-((durmax + xs[i]) / xs[i]))))
            tmp2 =  np.append(tmp2,(((1. + flux_err) * sens_maxgap * day1_obs) /  xs[i])   / (np.exp(-(max_distance / xs[i])) - np.exp(-(max_distance + day1_obs) / xs[i])))
        
        nottmp1 = (1. + flux_err) * sens_last * day1_obs / xs / (np.exp(-(durmax - day1_obs + xs) /  xs) - np.exp(-((durmax + xs) / xs)))
        nottmp2 = (((1. + flux_err) * sens_maxgap * day1_obs) /  xs)   / (np.exp(-(max_distance / xs)) - np.exp(-(max_distance + day1_obs) / xs))
        print("durmax_y, notdurmax_y, maxdist_y, notmaxdist_y\n")      
        for i in range(len(durmax_y)):
            print tmp1[i], nottmp1[i], tmp2[i], nottmp2[i]
    elif lightcurve == 'tophat':
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        maxdist_x = np.empty(len(ys))
        maxdist_x.fill(np.log10(max_distance))

    day1_obs_x = np.empty(len(ys))
    day1_obs_x.fill(np.log10(day1_obs))
    
    sensmin_y = np.empty(len(xs))
    sensmin_y.fill(np.log10(min_sens))
    
    sensmax_y = np.empty(len(xs))
    sensmax_y.fill(np.log10(max_sens))

    extra_y = np.empty(len(xs))
    extra_y.fill(np.log10(extra_thresh))

    ax = plt.gca()
    ax.set_xlim(dmin, dmax)
    ax.set_ylim(flmin, flmax)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    xi = np.linspace(dmin, dmax, 100)
    yi = np.linspace(flmin, flmax, 100)

    X, Y = np.meshgrid(xi, yi)
    Z = interpolate.griddata(toplot[:,0:2], toplot[:,2], (X,Y), method = 'linear')
    levels = np.linspace(0.000001, 1.01, 500)
#    levels = np.linspace(min(probabilities), max(probabilities), 100)

    surf = plt.contourf(X,Y,Z, levels=levels, cmap=cm.copper_r)
    cbar = plt.colorbar(surf,fraction=0.04, pad=0.01)
    cbar.set_ticklabels([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
#    cbar.ax.tick_params(labelsize=18, size=18)
    cbar.set_label(label='Probability',   weight='bold')
    
    if lightcurve == 'fred':
        ax.plot(xs, durmax_y, 'r-', color='r', linewidth=2)
        ax.plot(xs, maxdist_y, 'r-', color='r', linewidth=2)

    elif lightcurve == 'tophat':
        ax.plot(durmax_x, ys, 'r-', color='r', linewidth=2)
        ax.plot(maxdist_x, ys, 'r-', color='r', linewidth=2)

    ax.plot(day1_obs_x, ys, 'r-', color='r', linewidth=2)
    ax.plot(xs, sensmin_y, 'r-', color='r', linewidth=2)
    ax.plot(xs, sensmax_y, 'r-', color='r', linewidth=2)
    ax.plot(xs, extra_y, 'r-', color='r', linewidth=2)

    ax.tick_params(axis='both', direction='out')
    ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
    ax.get_yaxis().tick_left()

    plt.savefig(file + '_ProbContour.pdf')
    plt.close()

if __name__ == "__main__":
    initialise()
    print 'done'
