import configparser
import random
import time
import datetime
import numpy as np
import os
import argparse
import math
import warnings
import compute_lc
import importlib
from bitarray import bitarray
warnings.simplefilter("error", RuntimeWarning)

def get_configuration():
    """Reads in command line flags and returns a populated configuration"""
    
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--observations", help="Observation filename")    # format date YYYY-MM-DDTHH:MM:SS.mmmmmm, dur (day), sens (Jy)
    argparser.add_argument("--keep", action='store_true', help="Keep previous bursts")
    # argparser.add_argument("--stat_plot", action='store_true', help="Performs statistics and plots results")
    # argparser.add_argument("--just_plot", action='store_true', help="Just plots results")
    # argparser.add_argument("--just_stat", action='store_true', help="Just plots results")


    return argparser.parse_args()

def read_ini_file(filename = 'config.ini'):
    """Reads config.ini and returns a dictionary of its contents"""
    
    params = configparser.ConfigParser()
    params.read(filename)

    return params

def write_source(filename, bursts):
    """Writes simulated sources, returns nothing"""
    
    with open(filename, 'a') as f:
        for burst in bursts:
            f.write("{}\t{}\t{}\n".format(burst[0], burst[1], burst[2]))            
        f.flush()
        
def write_stat(filename, bursts):
    """Writes statistics on detections, returns nothing"""
    
    with open(filename, 'a') as f:
        for burst in bursts:
#            f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(burst[0], burst[1], burst[2], burst[3], burst[4]))        # for python 2.6 (krieger)
            f.write("{}\t{}\t{}\t{}\t{}\n".format(burst[0], burst[1], burst[2], burst[3], burst[4]))        # for python 2.7 (struis)
        f.flush()



#currently must be "tophat" or "fred" or "gaussian" or "wilma" or "ered" or 'halfgaussian1' 'parabolic' or 'choppedgaussian'
# Main execution starts here
params = read_ini_file("config.ini") # read config.ini into params
config = get_configuration() # read command line input 
lightcurvetype = params['INITIAL PARAMETERS']['lightcurvetype'] 
lightcurve_obj = getattr(importlib.import_module(lightcurvetype), lightcurvetype) # import the lightcurve class specified in the config file 
lightcurve = lightcurve_obj()
obs, pointFOV, regions = compute_lc.observing_strategy(config.observations, 
    np.float(params['INITIAL PARAMETERS']['det_threshold']), 
    int(params['SIM']['nobs']), 
    np.float(params['SIM']['obssens']), 
    np.float(params['SIM']['obssig']), 
    np.float(params['SIM']['obsinterval']), 
    np.float(params['SIM']['obsdurations']))

# pybtarr is a bitarray of size n_sources by n_pointings. Bitarrays are always 1D and I do all of one pointing first then all the other pointing.
ptbtarr = compute_lc.generate_pointings(np.uint(np.float((params['INITIAL PARAMETERS']['n_sources']))),
    pointFOV,
    regions)
# def generate_sources(n_sources, file, start_time, end_time, fl_min, fl_max, dmin, dmax, dump_intermediate, lightcurve,gaussiancutoff):
bursts = compute_lc.generate_sources(np.uint(np.float((params['INITIAL PARAMETERS']['n_sources']))), #n_sources
    obs[0,0], #start_survey
    obs[-1,0] + obs[-1,1], #end_survey
    np.float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
    np.float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
    np.float(params['INITIAL PARAMETERS']['dmin']), # duration min
    np.float(params['INITIAL PARAMETERS']['dmax']),  #duration max
    lightcurvetype) # 

bursts = compute_lc.generate_start(bursts, 
    lightcurve.earliest_crit_time(obs[0,0],bursts[:,1]), # earliest crit time
    lightcurve.latest_crit_time(obs[-1,0] + obs[-1,1],bursts[:,1]),  # latest crit time
    np.uint(np.float((params['INITIAL PARAMETERS']['n_sources']))))

if config.keep:
    with open(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans', 'w') as f:
            f.write('# Tcrit\tcharacteristic\tPkFlux\n')        ## INITIALISE LIST OF SIMILATED TRANSIENTS
            write_source(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans' , bursts) #file with starttime\tduration\tflux
            print("Written Simulated Sources")


# det are the sources themselves while detbool is a numpy boolean array indexing all sources
det, detbool = compute_lc.detect_bursts(obs, 
    np.float(params['INITIAL PARAMETERS']['flux_err']), 
    np.float(params['INITIAL PARAMETERS']['det_threshold']) , 
    np.float(params['INITIAL PARAMETERS']['extra_threshold']), 
    bursts, 
    2, # gaussiancutoff 
    lightcurve.edges,# edges present ?
    lightcurve.fluxint, 
    params['INITIAL PARAMETERS']['file'],
    config.keep,
    write_source,
    pointFOV,
    ptbtarr) 
    

    
# stat is a large numpy array of stats like probability and bins 
stat = compute_lc.statistics(np.float(params['INITIAL PARAMETERS']['fl_min']), 
    np.float(params['INITIAL PARAMETERS']['fl_max']), 
    np.float(params['INITIAL PARAMETERS']['dmin']), 
    np.float(params['INITIAL PARAMETERS']['dmax']), 
    det, 
    bursts)
    
if config.keep:
    with open(params['INITIAL PARAMETERS']['file'] + '_Stat', 'w') as f:
        f.write('# Duration\tFlux\tProbability\tDets\tAlls\n')        ## INITIALISE LIST OF STATISTICS
        write_stat(params['INITIAL PARAMETERS']['file'] + '_Stat', stat)
        print("Written Statistics")

#Plot a combined probability plot for all fields
compute_lc.plots(obs, 
    params['INITIAL PARAMETERS']['file'], 
    np.float(params['INITIAL PARAMETERS']['extra_threshold']), 
    np.float(params['INITIAL PARAMETERS']['det_threshold']), 
    np.float(params['INITIAL PARAMETERS']['flux_err']),  
    stat, 
    2,
    lightcurve.lines)
    
uniquepointFOV = np.unique(pointFOV, axis=0)
#iterate over the individual pointings making statistics and probability plots
for i in range(len(uniquepointFOV)):
    r = regions[i]
    detbtarr = bitarray(list(detbool)) 
    #totind = ptbtarr[i*len(bursts):(i+1)*len(bursts)].search(bitarray('1')) # This and the following gets the indices needed to index the sources that match the region
    detind = (detbtarr & ptbtarr[i*len(bursts):(i+1)*len(bursts)]).search(bitarray('1'))
    stat = compute_lc.statistics(np.float(params['INITIAL PARAMETERS']['fl_min']), 
    np.float(params['INITIAL PARAMETERS']['fl_max']), 
    np.float(params['INITIAL PARAMETERS']['dmin']), 
    np.float(params['INITIAL PARAMETERS']['dmax']), 
    bursts[detind], 
    bursts)

    if config.keep:
        with open(params['INITIAL PARAMETERS']['file'] + '_Stat', 'w') as f:
            f.write('# Duration\tFlux\tProbability\tDets\tAlls\n')        ## INITIALISE LIST OF STATISTICS
            write_stat(params['INITIAL PARAMETERS']['file'] +r['identity'] + '_Stat', stat)
            print("Written Statistics")

    compute_lc.plots(obs, 
        params['INITIAL PARAMETERS']['file']+r['identity'], 
        np.float(params['INITIAL PARAMETERS']['extra_threshold']), 
        np.float(params['INITIAL PARAMETERS']['det_threshold']), 
        np.float(params['INITIAL PARAMETERS']['flux_err']),  
        stat, 
        2,
        lightcurve.lines)
        
leftoff = len(uniquepointFOV) # Same as above, but now make plots for the overlapping 
for i in range(len(uniquepointFOV)):
    for j in range(i+1,len(uniquepointFOV)):
        r = regions[leftoff]
        detbtarr = bitarray(list(detbool))
        # totind = (ptbtarr[i*len(bursts):(i+1)*len(bursts)] & ptbtarr[j*len(bursts):(j+1)*len(bursts)]).search(bitarray('1'))
        detind = (detbtarr & (ptbtarr[i*len(bursts):(i+1)*len(bursts)] & ptbtarr[j*len(bursts):(j+1)*len(bursts)])).search(bitarray('1'))
        stat = compute_lc.statistics(np.float(params['INITIAL PARAMETERS']['fl_min']), 
        np.float(params['INITIAL PARAMETERS']['fl_max']), 
        np.float(params['INITIAL PARAMETERS']['dmin']), 
        np.float(params['INITIAL PARAMETERS']['dmax']), 
        bursts[detind], 
        bursts)

        if config.keep:
            with open(params['INITIAL PARAMETERS']['file'] + '_Stat', 'w') as f:
                f.write('# Duration\tFlux\tProbability\tDets\tAlls\n')        ## INITIALISE LIST OF STATISTICS
                write_stat(params['INITIAL PARAMETERS']['file'] +r['identity'] + '_Stat', stat)
                print("Written Statistics")

        compute_lc.plots(obs, 
            params['INITIAL PARAMETERS']['file']+r['identity'], 
            np.float(params['INITIAL PARAMETERS']['extra_threshold']), 
            np.float(params['INITIAL PARAMETERS']['det_threshold']), 
            np.float(params['INITIAL PARAMETERS']['flux_err']),  
            stat, 
            2,
            lightcurve.lines)
        leftoff += 1
