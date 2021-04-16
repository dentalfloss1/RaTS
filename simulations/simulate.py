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
from bokeh.plotting import figure, show, output_file, ColumnDataSource
from bokeh.models import LinearColorMapper, SingleIntervalTicker, ColorBar, Title, LogColorMapper, LogTicker, Range1d, HoverTool, Plot, VBar, LinearAxis, Grid, Line
from bokeh.io import export_png, curdoc, show
warnings.simplefilter("error", RuntimeWarning)

def get_configuration():
    """Reads in command line flags and returns a populated configuration"""
    
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--observations", help="Observation filename")    # format date YYYY-MM-DDTHH:MM:SS.mmmmmm, dur (day), sens (Jy)
    argparser.add_argument("--burstlength", help="All simulated transients this length")    # format date YYYY-MM-DDTHH:MM:SS.mmmmmm, dur (day), sens (Jy)
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
burstlength = np.float32(config.burstlength)
obs, pointFOV = compute_lc.observing_strategy(config.observations, 
    np.float(params['INITIAL PARAMETERS']['det_threshold']), 
    int(params['SIM']['nobs']), 
    np.float(params['SIM']['obssens']), 
    np.float(params['SIM']['obssig']), 
    np.float(params['SIM']['obsinterval']), 
    np.float(params['SIM']['obsdurations']))
uniquepointFOV = np.unique(pointFOV, axis=0)
regions, obssubsection = compute_lc.calculate_regions(pointFOV, obs)
overlapnums = []
for i in range(len(uniquepointFOV)):
    tsurvey = obs['start'][-1] + obs['duration'][-1] - obs['start'][0]
    startepoch = regions['start'][i]
    stopepoch = regions['stop'][i]
    targetnum = int(float(params['INITIAL PARAMETERS']['n_sources'])*(regions['timespan'][i]/tsurvey)) # Inner parenthesis is important for type conversion
    # pybtarr is a bitarray of size n_sources by n_pointings. Bitarrays are always 1D and I do all of one pointing first then all the other pointing.
    overlapnums.append(compute_lc.generate_pointings(targetnum,
        pointFOV,
        i))

    # exit()
    # def generate_sources(n_sources, file, start_time, end_time, fl_min, fl_max, dmin, dmax, dump_intermediate, lightcurve,gaussiancutoff):
    bursts = compute_lc.generate_sources(targetnum, #n_sources
        startepoch, #start_survey
        stopepoch, #end_survey
        np.float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
        np.float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
        np.float(params['INITIAL PARAMETERS']['dmin']), # duration min
        np.float(params['INITIAL PARAMETERS']['dmax']),  #duration max
        lightcurvetype,
        burstlength) # 

    bursts = compute_lc.generate_start(bursts, 
        lightcurve.earliest_crit_time(startepoch,bursts['chardur']), # earliest crit time
        lightcurve.latest_crit_time(stopepoch,bursts['chardur']),  # latest crit time
        targetnum)

    if config.keep:
        with open(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans', 'w') as f:
                f.write('# Tcrit\tcharacteristic\tPkFlux\n')        ## INITIALISE LIST OF SIMILATED TRANSIENTS
                write_source(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans' , bursts) #file with starttime\tduration\tflux
                print("Written Simulated Sources")


    # det are the sources themselves while detbool is a numpy boolean array indexing all sources
    det, detbool = compute_lc.detect_bursts(obs[obssubsection[i][0]:(obssubsection[i][1]+1)], 
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
        pointFOV) 
        

        
    # stat is a large numpy array of stats like probability and bins 
    stat = compute_lc.statistics(np.float(params['INITIAL PARAMETERS']['fl_min']), 
        np.float(params['INITIAL PARAMETERS']['fl_max']), 
        np.float(params['INITIAL PARAMETERS']['dmin']), 
        np.float(params['INITIAL PARAMETERS']['dmax']), 
        det, 
        bursts)
        
    # We now repeat all of the steps from simulating the sources to gathering statistics, but this time we 
    # use a single large value for transient duration and a single point in time for observations. We do this
    # to help eliminate false detections due to variations in observation sensitivity. 

    fake_obs = np.copy(obs)
    fake_obs['start'] = np.full(fake_obs['start'].shape, fake_obs['start'][0])
    import tophat
    tophatlc = tophat.tophat()
    fdbursts = compute_lc.generate_sources(targetnum, #n_sources
        startepoch, #start_survey
        stopepoch, #end_survey
        np.float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
        np.float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
        np.float(params['INITIAL PARAMETERS']['dmin']), # duration min
        np.float(params['INITIAL PARAMETERS']['dmax']),  #duration max
        "tophat",
        tsurvey) # 
    fdbursts['chartime'] += fake_obs['start'][0]

    # det are the sources themselves while detbool is a numpy boolean array indexing all sources
    fddet, fddetbool = compute_lc.detect_bursts(fake_obs[obssubsection[i][0]:(obssubsection[i][1]+1)], 
        np.float(params['INITIAL PARAMETERS']['flux_err']), 
        np.float(params['INITIAL PARAMETERS']['det_threshold']) , 
        np.float(params['INITIAL PARAMETERS']['extra_threshold']), 
        fdbursts, 
        2, # gaussiancutoff 
        tophatlc.edges,# edges present ?
        tophatlc.fluxint, 
        params['INITIAL PARAMETERS']['file'],
        False,
        write_source,
        pointFOV) 
        

        
    # stat is a large numpy array of stats like probability and bins 
    fdstat = compute_lc.statistics(np.float(params['INITIAL PARAMETERS']['fl_min']), 
        np.float(params['INITIAL PARAMETERS']['fl_max']), 
        np.float(params['INITIAL PARAMETERS']['dmin']), 
        np.float(params['INITIAL PARAMETERS']['dmax']), 
        fddet, 
        fdbursts)

    fl_max = np.float(params['INITIAL PARAMETERS']['fl_max'])
    fl_min = np.float(params['INITIAL PARAMETERS']['fl_min'])
    dmin = np.float(params['INITIAL PARAMETERS']['dmin'])
    dmax = np.float(params['INITIAL PARAMETERS']['dmax'])
    
    det_threshold = np.float(params['INITIAL PARAMETERS']['det_threshold'])
    extra_threshold = np.float(params['INITIAL PARAMETERS']['extra_threshold'])
    current_obs = obs[obssubsection[i][0]:(obssubsection[i][1]+1)]

    uldetections = int(params['INITIAL PARAMETERS']['uldetections'])


        
    
    
    cdet = fddet['charflux']

    mplretval = compute_lc.make_mpl_plots(regions['identity'][i],
                                        fl_min,
                                        fl_max,
                                        dmin,
                                        dmax,
                                        det_threshold,
                                        extra_threshold,
                                        current_obs,
                                        cdet,
                                        params['INITIAL PARAMETERS']['file'],
                                        np.float(params['INITIAL PARAMETERS']['flux_err']),
                                        np.copy(stat),
                                        2,
                                        lightcurve.lines,
                                        regions['area'][i],
                                        tsurvey,
                                        uldetections)
                                          

    
# exit()
    
overlaparray = np.array(overlapnums)
for i in range(len(uniquepointFOV),len(regions)):
    targetnum = np.sum(overlaparray[overlaparray[:,0]==i,1])
    if targetnum!=0:
        oindices = []
        for o in obssubsection:
            if o[2] in regions['identity'][i]:
                oindices.append([o[0],o[1]])
        tsurvey = obs['start'][-1] + obs['duration'][-1] - obs['start'][0]
        startepoch = regions['start'][i]
        stopepoch = regions['stop'][i]
        # pybtarr is a bitarray of size n_sources by n_pointings. Bitarrays are always 1D and I do all of one pointing first then all the other pointing.

        # exit()
        # def generate_sources(n_sources, file, start_time, end_time, fl_min, fl_max, dmin, dmax, dump_intermediate, lightcurve,gaussiancutoff):
        bursts = compute_lc.generate_sources(targetnum, #n_sources
            startepoch, #start_survey
            stopepoch, #end_survey
            np.float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
            np.float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
            np.float(params['INITIAL PARAMETERS']['dmin']), # duration min
            np.float(params['INITIAL PARAMETERS']['dmax']),  #duration max
            lightcurvetype,
            burstlength) # 

        bursts = compute_lc.generate_start(bursts, 
            lightcurve.earliest_crit_time(startepoch,bursts['chardur']), # earliest crit time
            lightcurve.latest_crit_time(stopepoch,bursts['chardur']),  # latest crit time
            targetnum)

        if config.keep:
            with open(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans', 'w') as f:
                    f.write('# Tcrit\tcharacteristic\tPkFlux\n')        ## INITIALISE LIST OF SIMILATED TRANSIENTS
                    write_source(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans' , bursts) #file with starttime\tduration\tflux
                    print("Written Simulated Sources")


        # det are the sources themselves while detbool is a numpy boolean array indexing all sources
        det, detbool = compute_lc.detect_bursts(obs[min(min(oindices)):(max(max(oindices))+1)], 
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
            pointFOV) 
            

            
        # stat is a large numpy array of stats like probability and bins 
        stat = compute_lc.statistics(np.float(params['INITIAL PARAMETERS']['fl_min']), 
            np.float(params['INITIAL PARAMETERS']['fl_max']), 
            np.float(params['INITIAL PARAMETERS']['dmin']), 
            np.float(params['INITIAL PARAMETERS']['dmax']), 
            det, 
            bursts)
                
        # We now repeat all of the steps from simulating the sources to gathering statistics, but this time we 
        # use a single large value for transient duration and a single point in time for observations. We do this
        # to help eliminate false detections due to variations in observation sensitivity. 

        fake_obs = np.copy(obs)
        fake_obs['start'] = np.full(fake_obs['start'].shape, fake_obs['start'][0])
        

        fdbursts = compute_lc.generate_sources(targetnum, #n_sources
            startepoch, #start_survey
            stopepoch, #end_survey
            np.float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
            np.float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
            np.float(params['INITIAL PARAMETERS']['dmin']), # duration min
            np.float(params['INITIAL PARAMETERS']['dmax']),  #duration max
            lightcurvetype,
            1e6*tsurvey) # 

        fdbursts = compute_lc.generate_start(fdbursts, 
            lightcurve.earliest_crit_time(startepoch,fdbursts['chardur']), # earliest crit time
            lightcurve.latest_crit_time(stopepoch,fdbursts['chardur']),  # latest crit time
            targetnum)

       
        # det are the sources themselves while detbool is a numpy boolean array indexing all sources
        fddet, fddetbool = compute_lc.detect_bursts(fake_obs[min(min(oindices)):(max(max(oindices))+1)], 
            np.float(params['INITIAL PARAMETERS']['flux_err']), 
            np.float(params['INITIAL PARAMETERS']['det_threshold']) , 
            np.float(params['INITIAL PARAMETERS']['extra_threshold']), 
            fdbursts, 
            2, # gaussiancutoff 
            lightcurve.edges,# edges present ?
            lightcurve.fluxint, 
            params['INITIAL PARAMETERS']['file'],
            False,
            write_source,
            pointFOV) 
            

            
        # stat is a large numpy array of stats like probability and bins 
        fdstat = compute_lc.statistics(np.float(params['INITIAL PARAMETERS']['fl_min']), 
            np.float(params['INITIAL PARAMETERS']['fl_max']), 
            np.float(params['INITIAL PARAMETERS']['dmin']), 
            np.float(params['INITIAL PARAMETERS']['dmax']), 
            fddet, 
            fdbursts)

        fl_max = np.float(params['INITIAL PARAMETERS']['fl_max'])
        fl_min = np.float(params['INITIAL PARAMETERS']['fl_min'])
        dmin = np.float(params['INITIAL PARAMETERS']['dmin'])
        dmax = np.float(params['INITIAL PARAMETERS']['dmax'])
        
        det_threshold = np.float(params['INITIAL PARAMETERS']['det_threshold'])
        extra_threshold = np.float(params['INITIAL PARAMETERS']['extra_threshold'])
        current_obs = obs[min(min(oindices)):(max(max(oindices))+1)]

        uldetections = int(params['INITIAL PARAMETERS']['uldetections'])


            
        
        
        cdet = fddet['charflux']

        mplretval = compute_lc.make_mpl_plots(regions['identity'][i].replace('&', 'and'),
                                            fl_min,
                                            fl_max,
                                            dmin,
                                            dmax,
                                            det_threshold,
                                            extra_threshold,
                                            current_obs,
                                            cdet,
                                            params['INITIAL PARAMETERS']['file'],
                                            np.float(params['INITIAL PARAMETERS']['flux_err']),
                                            np.copy(stat),
                                            2,
                                            lightcurve.lines,
                                            regions['area'][i],
                                            tsurvey,
                                            uldetections)
       
