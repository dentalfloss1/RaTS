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
rc('text', usetex=False)


def get_configuration():
    """Returns a populated configuration"""
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--observations", help="Observation filename")    # format date YYYY-MM-DDTHH:MM:SS.mmmmmm, dur (sec), sens (Jy)
    argparser.add_argument("--transient_type", help="Transient lightcurve")    # format: 'fred' or 'tophat'
    argparser.add_argument("--keep", action='store_true', help="Keep previous bursts")
    argparser.add_argument("--stat_plot", action='store_true', help="Performs statistics and plots results")
    argparser.add_argument("--just_plot", action='store_true', help="Just plots results")


    return argparser.parse_args()

def write_source(filename, bursts):
    with open(filename, 'a') as f:
        for burst in bursts:
            f.write("{}\t{}\t{}\n".format(burst[0], burst[1], burst[2]))            
        f.flush()

def initialise():
    
# Transients parameters
    n_sources = np.long(2e6)  # integer number of sources to be simulated
    fl_min = np.float(1e-5)     #Minimum simulated flux, in same units as the flux in the observations file
    fl_max = np.float(1)   #Maximum simulated flux ,in same units as the flux in the observations file
    flux_err = np.float(0.1)    #Fractional error in the simulated flux
    dmin = np.float(0.001)      #Minimum simulated duration, in same units as the duration in the observations file
    dmax = np.float(1e3)    #Maximum simulated duration, in same units as the duration in the observations file
    det_threshold = np.float(5) #detection threshold, to be multiplied by the noise in the images
    extra_threshold = np.float(3)  #integer, extra detection threshold, to be multiplied by the noise in the images -- used to be extra certain of the transient sources
    file = "output"   #Name to be used for the output files

    config = get_configuration() # parses arguments, see above
    if config.transient_type != 'tophat' and config.transient_type != 'fred':
        print('Type of transient not recognised.\nUse tophat or fred.')
        exit()
    file = file + '_' + config.transient_type

    obs = observing_strategy(config.observations, det_threshold)

    if config.just_plot: # Not sure why this code exists, replots a previously calculated file named file above
        plots(obs, file, extra_threshold, det_threshold, flux_err, config.transient_type)
        print('done')
        exit()        

    if not config.keep: # if you to save the burst we save it:
        with open(file + '_SimTrans', 'w') as f:
            f.write('# Tstart\tDuration\tFlux\n')        ## INITIALISE LIST OF SIMILATED TRANSIENTS
    
        with open(file + '_DetTrans', 'w') as f:
            f.write('# Tstart\tDuration\tFlux\n')        ## INITIALISE LIST OF DETECTED TRANSIENTS
    

    start_time = obs[0][0] #recall that obs has the form: start, duration, sensitivity. Therefore this is the start of the very first observation.
    end_time = obs[-1][0] + obs[-1][1] #The time that the last observation started + its duration. end_time-start_time = duration of survey

    simulated = generate_sources(n_sources, file, start_time, end_time, fl_min, fl_max, dmin, dmax) # two functions down
    detect_bursts(obs, file, flux_err, det_threshold, extra_threshold, simulated, config.transient_type)
    
    if config.stat_plot:
        with open(file + '_Stat', 'w') as f:
            f.write('# Duration\tFlux\tProbability\n')        ## INITIALISE LIST OF STATISTICS
        statistics(file, fl_min, fl_max, dmin, dmax)
        plots(obs, file, extra_threshold, det_threshold, flux_err, config.transient_type)

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
    #The following two functions pick a random number that is evenly spaced logarithmically
    bursts[:,1] = np.power(10, np.random.uniform(np.log10(dmin), np.log10(dmax), n_sources)) # random number for duration
    bursts[:,2] = np.power(10, np.random.uniform(np.log10(fl_min), np.log10(fl_max), n_sources)) # random number for flux
    bursts[:,0] = np.random.uniform(start_time - bursts[:,1], end_time, n_sources) # randomize the start times based partly on the durations above
    
    bursts = bursts[bursts[:,0].argsort()] # Order on start times
    write_source(file + '_SimTrans', bursts) #file with starttime\tduration\tflux
    
    print("Written Simulated Sources")
    return bursts

def detect_bursts(obs, file, flux_err, det_threshold, extra_threshold, sources, lightcurve):
    single_detection = np.array([],dtype=np.uint32)
    extra_detection = np.array([],dtype=np.uint32)

    for o in obs:
        start_obs, duration, sensitivity = o
        extra_sensitivity = sensitivity * (det_threshold + extra_threshold) / det_threshold
        end_obs = start_obs + duration
        
        # The following handles edge cases
        if lightcurve == 'fred':
            ## Transients start before observation ends
            single_candidate = np.where(sources[:,0] < end_obs)[0]
        elif lightcurve == 'tophat':
            ## Transients end after observation starts and start before observation ends
            single_candidate = np.where((sources[:,0] + sources[:,1] > start_obs) & (sources[:,0] < end_obs))[0] # Take all of the locations that qualify as a candidate. The zero index is a wierd python workaround

        # filter on integrated flux
        F0_o = sources[single_candidate][:,2]
        error = np.sqrt((abs(random.gauss(F0_o * flux_err, 0.05 * F0_o * flux_err)))**2 + (sensitivity/det_threshold)**2) 
        F0 = random.gauss(F0_o, error) # Simulate some variation in source flux of each source.

        tau = sources[single_candidate][:,1] # Durations
        t_burst = sources[single_candidate][:,0] # start times
        tstart = np.maximum(t_burst, start_obs) - t_burst # How much time passed in which the transient was on, but not being observed in the survey.

        if lightcurve == 'fred':
            tend = end_obs - t_burst # Burst is never really "Off" so, how long is it on? 
            flux_int = np.multiply(F0, np.multiply(tau, np.divide(np.exp(-np.divide(tstart,tau)) - np.exp(-np.divide(tend,tau)), (end_obs-start_obs))))
        
        elif lightcurve == 'tophat':
            tend = np.minimum(t_burst + tau, end_obs) - t_burst
            flux_int = np.multiply(F0, np.divide((tend - tstart), (end_obs-start_obs)))
        
        candidates = single_candidate[(flux_int > sensitivity)]
        extra_candidates = np.array(single_candidate[flux_int > extra_sensitivity])

        single_detection = np.append(single_detection, candidates)
        extra_detection = np.append(extra_detection, extra_candidates)
        
    dets = unique_count(single_detection)
    detections = dets[0][np.where(dets[1] < len(obs))[0]]
    detections = detections[(np.in1d(detections, extra_detection))]
    write_source(file + '_DetTrans', sources[detections])
    
    print("Written Detected Sources")
    return sources[detections]

def unique_count(a):
    unique, inverse = np.unique(a, return_inverse=True)
    count = np.zeros(len(unique), np.int)
    np.add.at(count, inverse, 1)
    return [unique, count]


def statistics(file, fl_min, fl_max, dmin, dmax):
    all = np.loadtxt(file + '_SimTrans')
    det = np.loadtxt(file + '_DetTrans')

    flux_ints = np.linspace(np.log10(fl_min), np.log10(fl_max), num=(np.log10(fl_max)-np.log10(fl_min))/0.05, endpoint=True)
    dur_ints = np.linspace(np.log10(dmin), np.log10(dmax), num=(np.log10(dmax)-np.log10(dmin))/0.05, endpoint=True)

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

        all_dur = np.append(all_dur, np.where((np.log10(all[:,1]) >= dur_ints[i]) & (np.log10(all[:,1]) < dur_ints[i+1]))[0])
        det_dur = np.append(det_dur, np.where((np.log10(det[:,1]) >= dur_ints[i]) & (np.log10(det[:,1]) < dur_ints[i+1]))[0])
        durations = np.append(durations, np.power(10, (dur_ints[i] + dur_ints[i+1])/2.))

        for m in range(len(flux_ints) - 1):
            dets = np.append(dets, float(len(np.where((np.log10(det[det_dur,2]) >= flux_ints[m]) & (np.log10(det[det_dur,2]) < flux_ints[m+1]))[0])))
            alls = np.append(alls, float(len(np.where((np.log10(all[all_dur,2]) >= flux_ints[m]) & (np.log10(all[all_dur,2]) < flux_ints[m+1]))[0])))

            if doit_f == True:
                fluxes = np.append(fluxes, np.power(10, (flux_ints[m] + flux_ints[m+1])/2.))

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

    print("Written Statistics")

    return stats

def write_stat(filename, bursts):
    with open(filename, 'a') as f:
        for burst in bursts:
#            f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(burst[0], burst[1], burst[2], burst[3], burst[4]))        # for python 2.6 (krieger)
            f.write("{}\t{}\t{}\t{}\t{}\n".format(burst[0], burst[1], burst[2], burst[3], burst[4]))        # for python 2.7 (struis)
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
    pylab.xlabel(r'{Log Transient Duration [days]', {'color':'k'})
    pylab.ylabel(r'{Log Transient Flux Density [Jy]', {'color':'k'})
    pylab.xticks(fontsize=18)
    pylab.yticks(fontsize=18)

    dmin=min(toplot[:,0])
    dmax=max(toplot[:,0])
    flmin=min(toplot[:,1])
    flmax=max(toplot[:,1])

    xs = np.arange(dmin, dmax, 0.01, dtype = np.float32)
    ys = np.arange(flmin, flmax, 0.01, dtype = np.float32)
    print(dmin)
    print(xs)
    if lightcurve == 'fred':
        durmax_y = np.log10((1. + flux_err) * sens_last * day1_obs / np.power(10, xs) / (np.exp(-(durmax - day1_obs + np.power(10, xs)) / np.power(10, xs)) - np.exp(-((durmax + np.power(10, xs)) / np.power(10, xs)))))
        maxdist_y = np.log10((1. + flux_err) * sens_maxgap * day1_obs / np.power(10, xs) / (np.exp(-(max_distance / np.power(10, xs))) - np.exp(-(max_distance + day1_obs) / np.power(10, xs))))
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

    X, Y = np.mgrid[dmin:dmax:100j, flmin:flmax:100j]
    Z = interpolate.griddata(toplot[:,0:2], toplot[:,2], (X, Y), method='linear')

    levels = np.linspace(0.000001, 1.01, 500)
#    levels = np.linspace(min(probabilities), max(probabilities), 100)

    surf = plt.contourf(X,Y,Z, levels = levels, cmap=cm.copper_r)
    for c in surf.collections:
        c.set_edgecolor("face")
    cbar = plt.colorbar(surf,fraction=0.04, pad=0.01)
    cbar.set_ticklabels([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    cbar.set_label(label='Probability', weight='bold')
    
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

    plt.savefig(file + '_ProbContour.pdf',dpi=400)
    plt.close()

if __name__ == "__main__":
    initialise()
    print('done')
