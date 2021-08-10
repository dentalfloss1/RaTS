import configparser
import datetime
import numpy as np
import os
import argparse
import warnings
import compute_lc
import importlib
from bitarray import bitarray
# warnings.simplefilter("error", RuntimeWarning)

start = datetime.datetime.now()
def get_configuration():
    """Reads in command line flags and returns a populated configuration"""
    
    argparser = argparse.ArgumentParser()
    argparser.add_argument("--observations", help="Observation filename")    # format date YYYY-MM-DDTHH:MM:SS.mmmmmm, dur (day), sens (Jy)
    argparser.add_argument("--burstlength", help="All simulated transients this length (days)")   
    argparser.add_argument("--burstflux", help="All simulated transients this flux (Jy)")    
    argparser.add_argument("--keep", action='store_true', help="Keep previous bursts")
    argparser.add_argument("--configfile", default='config.ini', help="Configuration file. Default is config.ini")
    # argparser.add_argument("--stat_plot", action='store_true', help="Performs statistics and plots results")
    # argparser.add_argument("--just_plot", action='store_true', help="Just plots results")
    # argparser.add_argument("--just_stat", action='store_true', help="Just plots results")


    return argparser.parse_args()

def read_ini_file(filename):
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



#currently must be "tophat" or "fred" or "gaussian" or "wilma" or "ered" or  'parabolic' or 'choppedgaussian'
# Main execution starts here
config = get_configuration() # read command line input 
params = read_ini_file(config.configfile) # read config.ini into params
lightcurvetype = params['INITIAL PARAMETERS']['lightcurvetype'] 
lightcurve_obj = getattr(importlib.import_module(lightcurvetype), lightcurvetype) # import the lightcurve class specified in the config file 
lightcurve = lightcurve_obj()
burstlength = np.float32(config.burstlength)
burstflux = np.float32(config.burstflux)
obs, pointFOV = compute_lc.observing_strategy(config.observations, 
    float(params['INITIAL PARAMETERS']['det_threshold']), 
    int(params['SIM']['nobs']), 
    float(params['SIM']['obssens']), 
    float(params['SIM']['obssig']), 
    float(params['SIM']['obsinterval']), 
    float(params['SIM']['obsdurations']))
uniquepointFOV = np.unique(pointFOV, axis=0)
regions, obssubsection = compute_lc.calculate_regions(pointFOV, obs)
leftoff = len(uniquepointFOV)
overlapnums = np.zeros(len(regions), dtype={'names': ('name', 'sources'), 'formats': ('<U6','f8')})
obsmask = np.zeros((len(obs),len(uniquepointFOV)),dtype=bool)
for i in range(len(uniquepointFOV)):
    obsmask[:,i] = [np.all(p) for p in pointFOV[obssubsection[i][0]]==pointFOV]
    tsurvey = obs['start'][-1] + obs['duration'][-1] - obs['start'][0]
    startepoch = regions['start'][i]
    stopepoch = regions['stop'][i]
    targetnum = int(float(params['INITIAL PARAMETERS']['n_sources'])) # Inner parenthesis is important for type conversion
    # pybtarr is a bitarray of size n_sources by n_pointings. Bitarrays are always 1D and I do all of one pointing first then all the other pointing.
    # overlapnums, leftoff = compute_lc.generate_pointings(targetnum,
    #     pointFOV,
    #     i,
    #     leftoff,
    #     overlapnums) 
    # exit()
    # def generate_sources(n_sources, file, start_time, end_time, fl_min, fl_max, dmin, dmax, dump_intermediate, lightcurve,gaussiancutoff):
    bursts = compute_lc.generate_sources(targetnum, #n_sources
        startepoch, #start_survey
        stopepoch, #end_survey
        float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
        float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
        float(params['INITIAL PARAMETERS']['dmin']), # duration min
        float(params['INITIAL PARAMETERS']['dmax']),  #duration max
        lightcurvetype,
        burstlength,
        burstflux) # 

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
    det, detbool = compute_lc.detect_bursts(obs[obsmask[:,i]], 
        float(params['INITIAL PARAMETERS']['flux_err']), 
        float(params['INITIAL PARAMETERS']['det_threshold']) , 
        float(params['INITIAL PARAMETERS']['extra_threshold']), 
        bursts, 
        2, # gaussiancutoff 
        lightcurve.edges,# edges present ?
        lightcurve.fluxint, 
        params['INITIAL PARAMETERS']['file'],
        config.keep,
        write_source,
        pointFOV) 
        

        
    # stat is a large numpy array of stats like probability and bins 
    stat = compute_lc.statistics(float(params['INITIAL PARAMETERS']['fl_min']), 
        float(params['INITIAL PARAMETERS']['fl_max']), 
        float(params['INITIAL PARAMETERS']['dmin']), 
        float(params['INITIAL PARAMETERS']['dmax']), 
        det, 
        bursts)
        

    if not (np.isnan(burstlength) and np.isnan(burstflux)):
        print("Percent detected:", np.sum(detbool)/targetnum)
        print(np.sum(detbool))
        gaptime = 0
        for j in range(len(obs[obsmask[:,i]])):
            if obs[obsmask[:,i]]['gaps'][j]=='False':
                if j!=(len(obs[obsmask[:,i]])-1):
                    gaptime += obs[obsmask[:,i]]['start'][j+1] - (obs[obsmask[:,i]]['start'][j] + obs[obsmask[:,i]]['duration'][j])
            else:
                subobs, _ = compute_lc.observing_strategy(obs[obsmask[:,i]]['gaps'][j], float(params['INITIAL PARAMETERS']['det_threshold']), 1, 1, 1, 1, 1) # We are giving the scansfile name, so the other variables are unimportant, we set them to 1 
                for k in range(len(subobs)-1):
                    gaptime += subobs['start'][k+1] - (subobs['start'][k] + subobs['duration'][k])
                    if (k==len(subobs)-2) and (j!=len(obs[obsmask[:,i]])-1):
                        gaptime += obs[obsmask[:,i]]['start'][j+1] - (subobs['start'][k+1] + subobs['duration'][k+1])
        print("Gap percentage:", gaptime/(obs[obsmask[:,i]]['start'][-1] + obs[obsmask[:,i]]['duration'][-1] - obs[obsmask[:,i]]['start'][0]))
    else:
        # We now repeat all of the steps from simulating the sources to gathering statistics, but this time we 
        # use a single large value for transient duration and a single point in time for observations. We do this
        # to help eliminate false detections due to variations in observation sensitivity. 

        fake_obs = np.copy(obs)
        fake_obs['start'] = np.full(fake_obs['start'].shape, fake_obs['start'][0])
        fake_obs['gaps'] = 'False'
        import tophat
        tophatlc = tophat.tophat()
        fdbursts = compute_lc.generate_sources(targetnum, #n_sources
            startepoch, #start_survey
            stopepoch, #end_survey
            float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
            float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
            float(params['INITIAL PARAMETERS']['dmin']), # duration min
            float(params['INITIAL PARAMETERS']['dmax']),  #duration max
            "tophat",
            tsurvey,
            burstflux) # 
        fdbursts['chartime'] += fake_obs['start'][0]

        # det are the sources themselves while detbool is a numpy boolean array indexing all sources
        fddet, fddetbool = compute_lc.detect_bursts(fake_obs[obsmask[:,i]], 
            float(params['INITIAL PARAMETERS']['flux_err']), 
            float(params['INITIAL PARAMETERS']['det_threshold']) , 
            float(params['INITIAL PARAMETERS']['extra_threshold']), 
            fdbursts, 
            2, # gaussiancutoff 
            tophatlc.edges,# edges present ?
            tophatlc.fluxint, 
            params['INITIAL PARAMETERS']['file'],
            False,
            write_source,
            pointFOV) 
            

            
        # stat is a large numpy array of stats like probability and bins 
        fdstat = compute_lc.statistics(float(params['INITIAL PARAMETERS']['fl_min']), 
            float(params['INITIAL PARAMETERS']['fl_max']), 
            float(params['INITIAL PARAMETERS']['dmin']), 
            float(params['INITIAL PARAMETERS']['dmax']), 
            fddet, 
            fdbursts)

        fl_max = float(params['INITIAL PARAMETERS']['fl_max'])
        fl_min = float(params['INITIAL PARAMETERS']['fl_min'])
        dmin = float(params['INITIAL PARAMETERS']['dmin'])
        dmax = float(params['INITIAL PARAMETERS']['dmax'])
        
        det_threshold = float(params['INITIAL PARAMETERS']['det_threshold'])
        extra_threshold = float(params['INITIAL PARAMETERS']['extra_threshold'])
        current_obs = obs[obsmask[:,i]]

        
        detections = int(params['INITIAL PARAMETERS']['detections'])
        confidence = float(params['INITIAL PARAMETERS']['confidence'])/100

        cdet = fddet['charflux']
        compute_lc.make_mpl_plots(regions['identity'][i].replace('&', 'and'),
        fl_min,
        fl_max,
        dmin,
        dmax,
        det_threshold,
        extra_threshold,
        current_obs,
        cdet,
        params['INITIAL PARAMETERS']['file'],
        float(params['INITIAL PARAMETERS']['flux_err']),
        np.copy(stat),
        2,
        lightcurve.lines,
        regions['area'][i],
        tsurvey,
        detections,
        confidence,
        params['INITIAL PARAMETERS']['file'])

leftoff=0
obsmaskmulti = np.zeros((len(obs),len(regions)-len(uniquepointFOV)), dtype=bool)                    
for i in range(len(uniquepointFOV)-1):
    for j in range(i+1,len(uniquepointFOV)):
        if str(i)+'&'+str(j) in regions['identity']:
            obsmaskmulti[:,leftoff] = obsmask[:,i] | obsmask[:,j]
            leftoff+=1

for i in range(len(uniquepointFOV)-2):
    for j in range(i+1,len(uniquepointFOV)-1):
        for k in range(j+1,len(uniquepointFOV)):
            if str(i)+'&'+str(j)+'&'+str(k) in regions['identity']:
                obsmaskmulti[:,leftoff] = obsmask[:,i] | obsmask[:,j] | obsmask[:,k]

# overlaparray = np.array(len(overlapnums), dtype={'names': ('name', 'sources'), 'formats': ('str','f8')})
for i in range(len(uniquepointFOV),len(regions)):
    # targetnum = int(np.around(float(params['INITIAL PARAMETERS']['n_sources'])*(regions['area'][i]/regions['area'][0]))*(regions['timespan'][i]/regions['timespan'][0]))
    tsurvey = obs['start'][-1] + obs['duration'][-1] - obs['start'][0]
    startepoch = regions['start'][i]
    stopepoch = regions['stop'][i]
    # pybtarr is a bitarray of size n_sources by n_pointings. Bitarrays are always 1D and I do all of one pointing first then all the other pointing.

    # exit()
    # def generate_sources(n_sources, file, start_time, end_time, fl_min, fl_max, dmin, dmax, dump_intermediate, lightcurve,gaussiancutoff):
    bursts = compute_lc.generate_sources(targetnum, #n_sources
        startepoch, #start_survey
        stopepoch, #end_survey
        float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
        float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
        float(params['INITIAL PARAMETERS']['dmin']), # duration min
        float(params['INITIAL PARAMETERS']['dmax']),  #duration max
        lightcurvetype,
        burstlength,
        burstflux) # 

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
    det, detbool = compute_lc.detect_bursts(obs[obsmaskmulti[:,i-len(uniquepointFOV)]], 
        float(params['INITIAL PARAMETERS']['flux_err']), 
        float(params['INITIAL PARAMETERS']['det_threshold']) , 
        float(params['INITIAL PARAMETERS']['extra_threshold']), 
        bursts, 
        2, # gaussiancutoff 
        lightcurve.edges,# edges present ?
        lightcurve.fluxint, 
        params['INITIAL PARAMETERS']['file'],
        config.keep,
        write_source,
        pointFOV) 
        

        
    # stat is a large numpy array of stats like probability and bins 
    stat = compute_lc.statistics(float(params['INITIAL PARAMETERS']['fl_min']), 
        float(params['INITIAL PARAMETERS']['fl_max']), 
        float(params['INITIAL PARAMETERS']['dmin']), 
        float(params['INITIAL PARAMETERS']['dmax']), 
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
        float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
        float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
        float(params['INITIAL PARAMETERS']['dmin']), # duration min
        float(params['INITIAL PARAMETERS']['dmax']),  #duration max
        "tophat",
        tsurvey,
        burstflux) # 
    fdbursts['chartime'] += fake_obs['start'][0]

    # det are the sources themselves while detbool is a numpy boolean array indexing all sources
    fddet, fddetbool = compute_lc.detect_bursts(fake_obs[obsmaskmulti[:,i-len(uniquepointFOV)]], 
        float(params['INITIAL PARAMETERS']['flux_err']), 
        float(params['INITIAL PARAMETERS']['det_threshold']) , 
        float(params['INITIAL PARAMETERS']['extra_threshold']), 
        fdbursts, 
        2, # gaussiancutoff 
        tophatlc.edges,# edges present ?
        tophatlc.fluxint, 
        params['INITIAL PARAMETERS']['file'],
        False,
        write_source,
        pointFOV) 
        

        
    # stat is a large numpy array of stats like probability and bins 
    fdstat = compute_lc.statistics(float(params['INITIAL PARAMETERS']['fl_min']), 
        float(params['INITIAL PARAMETERS']['fl_max']), 
        float(params['INITIAL PARAMETERS']['dmin']), 
        float(params['INITIAL PARAMETERS']['dmax']), 
        fddet, 
        fdbursts)

    fl_max = float(params['INITIAL PARAMETERS']['fl_max'])
    fl_min = float(params['INITIAL PARAMETERS']['fl_min'])
    dmin = float(params['INITIAL PARAMETERS']['dmin'])
    dmax = float(params['INITIAL PARAMETERS']['dmax'])
    
    det_threshold = float(params['INITIAL PARAMETERS']['det_threshold'])
    extra_threshold = float(params['INITIAL PARAMETERS']['extra_threshold'])
    current_obs = obs[obsmaskmulti[:,i-len(uniquepointFOV)]]

    detections = int(params['INITIAL PARAMETERS']['detections'])
    confidence = float(params['INITIAL PARAMETERS']['confidence'])/100

        
    
    
    cdet = fddet['charflux']

    compute_lc.make_mpl_plots(regions['identity'][i].replace('&', 'and'),
        fl_min,
        fl_max,
        dmin,
        dmax,
        det_threshold,
        extra_threshold,
        current_obs,
        cdet,
        params['INITIAL PARAMETERS']['file'],
        float(params['INITIAL PARAMETERS']['flux_err']),
        np.copy(stat),
        2,
        lightcurve.lines,
        regions['area'][i],
        tsurvey,
        detections,
        confidence,
        params['INITIAL PARAMETERS']['file'])

end = datetime.datetime.now()
print("total runtime: ", end - start)