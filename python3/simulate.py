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
warnings.simplefilter("error", RuntimeWarning)

def get_configuration():
    """Returns a populated configuration"""
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--observations", help="Observation filename")    # format date YYYY-MM-DDTHH:MM:SS.mmmmmm, dur (day), sens (Jy)
    argparser.add_argument("--keep", action='store_true', help="Keep previous bursts")
    # argparser.add_argument("--stat_plot", action='store_true', help="Performs statistics and plots results")
    # argparser.add_argument("--just_plot", action='store_true', help="Just plots results")
    # argparser.add_argument("--just_stat", action='store_true', help="Just plots results")


    return argparser.parse_args()

def read_ini_file(filename = 'config.ini'):
    params = configparser.ConfigParser()
    params.read(filename)

    return params

def write_source(filename, bursts):
    with open(filename, 'a') as f:
        for burst in bursts:
            f.write("{}\t{}\t{}\n".format(burst[0], burst[1], burst[2]))            
        f.flush()
        
def write_stat(filename, bursts):
    with open(filename, 'a') as f:
        for burst in bursts:
#            f.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(burst[0], burst[1], burst[2], burst[3], burst[4]))        # for python 2.6 (krieger)
            f.write("{}\t{}\t{}\t{}\t{}\n".format(burst[0], burst[1], burst[2], burst[3], burst[4]))        # for python 2.7 (struis)
        f.flush()



#currently must be "tophat" or "fred" or "gaussian" or "wilma" or "ered" or 'halfgaussian1' 'parabolic' or 'choppedgaussian'

params = read_ini_file("config.ini")
config = get_configuration()
lightcurvetype = params['INITIAL PARAMETERS']['lightcurvetype']
lightcurve_obj = getattr(importlib.import_module(lightcurvetype), lightcurvetype)
lightcurve = lightcurve_obj()
obs = compute_lc.observing_strategy(config.observations,
    np.float(params['INITIAL PARAMETERS']['det_threshold']), 
    int(params['SIM']['nobs']), 
    np.float(params['SIM']['obssens']), 
    np.float(params['SIM']['obssig']), 
    int(params['SIM']['obsinterval']), 
    np.float(params['SIM']['obsdurations']))

# def generate_sources(n_sources, file, start_time, end_time, fl_min, fl_max, dmin, dmax, dump_intermediate, lightcurve,gaussiancutoff):
bursts = compute_lc.generate_sources(np.uint(np.float((params['INITIAL PARAMETERS']['n_sources']))), #n_sources
    obs[0][0], #start_survey
    obs[-1][0] + obs[-1][1], #end_survey
    np.float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
    np.float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
    np.float(params['INITIAL PARAMETERS']['dmin']), # duration min
    np.float(params['INITIAL PARAMETERS']['dmax']),  #duration max
    lightcurvetype) # 

bursts = compute_lc.generate_start(bursts, 
    lightcurve.earliest_crit_time(obs[0][0],bursts[:,1]), # earliest crit time
    lightcurve.latest_crit_time(obs[-1][0] + obs[-1][1],bursts[:,1]),  # latest crit time
    np.uint(np.float((params['INITIAL PARAMETERS']['n_sources']))))

if config.keep:
    with open(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans', 'w') as f:
            f.write('# Tcrit\tcharacteristic\tPkFlux\n')        ## INITIALISE LIST OF SIMILATED TRANSIENTS
            write_source(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans' , bursts) #file with starttime\tduration\tflux
            print("Written Simulated Sources")



det = compute_lc.detect_bursts(obs, 
    np.float(params['INITIAL PARAMETERS']['flux_err']), 
    np.float(params['INITIAL PARAMETERS']['det_threshold']) , 
    np.float(params['INITIAL PARAMETERS']['extra_threshold']), 
    bursts, 
    2, # gaussiancutoff 
    lightcurve.edges,# edges present ?
    lightcurve.fluxint ) 
    
if config.keep:
        with open(params['INITIAL PARAMETERS']['file'] + '_DetTrans', 'w') as f:
            f.write('# Tstart\tDuration\tFlux\n')        ## INITIALISE LIST OF DETECTED TRANSIENTS
            write_source(params['INITIAL PARAMETERS']['file'] + '_DetTrans', sources[detections])
            print("Written Detected Sources")
    
if  lightcurvetype == 'ered': #lightcurvetype == 'gaussian' or lightcurvetype == 'choppedgaussian' or
    correctivefactor = 0.5
else:
    correctivefactor = 1
stat = compute_lc.statistics(np.float(params['INITIAL PARAMETERS']['fl_min']), 
    np.float(params['INITIAL PARAMETERS']['fl_max']), 
    correctivefactor*np.float(params['INITIAL PARAMETERS']['dmin']), 
    correctivefactor*np.float(params['INITIAL PARAMETERS']['dmax']), 
    det, 
    bursts)
    
if config.keep:
    with open(params['INITIAL PARAMETERS']['file'] + '_Stat', 'w') as f:
        f.write('# Duration\tFlux\tProbability\n')        ## INITIALISE LIST OF STATISTICS
        write_stat(params['INITIAL PARAMETERS']['file'] + '_Stat', stats)
        print("Written Statistics")

compute_lc.plots(obs, 
    params['INITIAL PARAMETERS']['file'], 
    np.float(params['INITIAL PARAMETERS']['extra_threshold']), 
    np.float(params['INITIAL PARAMETERS']['det_threshold']), 
    np.float(params['INITIAL PARAMETERS']['flux_err']),  
    stat, 
    2,
    lightcurve.lines)
