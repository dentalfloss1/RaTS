#!python
import configparser
import datetime
import numpy as np
import os
import argparse
import warnings
from RaTS import compute_lc
from RaTS import tophat
import importlib
from bitarray import bitarray
from tqdm import tqdm
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


if __name__=='__main__':
    #currently must be "tophat" or "fred" or "gaussian" or "wilma" or "ered" or  'parabolic' or 'choppedgaussian'
    # Main execution starts here
    config = get_configuration() # read command line input 
    params = read_ini_file(config.configfile) # read config.ini into params
    try: 
        lightcurvetype = params['INITIAL PARAMETERS']['lightcurvetype'] 
    except KeyError:
        print("No config.ini file found. Creating one with default values.")
        from RaTS.make_templates import configfilestring
        with open("config.ini","w") as f:
            f.write(configfilestring)
        exit()
    lightcurve_obj = getattr(importlib.import_module(f"RaTS.{lightcurvetype}"), lightcurvetype) # import the lightcurve class specified in the config file 
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
    tsurveylist= []
    fl_min = float(params['INITIAL PARAMETERS']['fl_min'])
    fl_max = float(params['INITIAL PARAMETERS']['fl_max'])
    dmin = float(params['INITIAL PARAMETERS']['dmin'])
    dmax = float(params['INITIAL PARAMETERS']['dmax'])
    for i in range(len(uniquepointFOV)):
        obsmask[:,i] = [np.all(p) for p in pointFOV[obssubsection[i][0]]==pointFOV]
        tsurvey = obs[obsmask[:,i]][-1][0] + obs[obsmask[:,i]][-1][1] - obs[obsmask[:,i]][0][0] 
        tsurveylist.append(tsurvey)
        startepoch = regions['start'][i]
        stopepoch = regions['stop'][i]
        targetnum = int(float(params['INITIAL PARAMETERS']['srcperbin'])) # Inner parenthesis is important for type conversion
        # Create bins
        flux_bins = np.geomspace(fl_min, fl_max, num=int(round((np.log10(fl_max)-np.log10(fl_min))/0.05)), endpoint=True)
        dur_ints = np.geomspace(dmin, dmax, num=int(round((np.log10(dmax)-np.log10(dmin))/0.05)), endpoint=True)
        durations = (dur_ints[:-1] + dur_ints[1:])/2
        fluxes = (flux_bins[:-1] + flux_bins[1:])/2
        stats = np.zeros(((len(flux_bins)-1)*(len(dur_ints)-1), 5), dtype=np.float32)
        # iterate over bins
        statcounter = 0
        srcsimtime = 0
        dettime = 0
        stattime = 0
        for ldurbin,rdurbin in tqdm(zip(dur_ints[:-1],dur_ints[1:]),total=len(dur_ints[:-1])):
            thisdur = (ldurbin+rdurbin)/2
            for lfluxbin,rfluxbin in zip(flux_bins[:-1],flux_bins[1:]):
                t1 = datetime.datetime.now()
                thisflux = (lfluxbin + rfluxbin)/2
                bursts = compute_lc.generate_sources(targetnum, #n_sources
                    startepoch, #start_survey
                    stopepoch, #end_survey
                    lfluxbin, #Flux min
                    rfluxbin, #Flux max
                    ldurbin, # duration min
                    rdurbin,  #duration max
                    lightcurvetype,
                    burstlength,
                    burstflux) # 
    
                bursts = compute_lc.generate_start(bursts, 
                    lightcurve.earliest_crit_time(startepoch,bursts['chardur']), # earliest crit time
                    lightcurve.latest_crit_time(stopepoch,bursts['chardur']),  # latest crit time
                    targetnum)
                t2 = datetime.datetime.now() 
                srcsimtime+= (t2-t1).total_seconds()
                if config.keep:
                    with open(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans', 'w') as f:
                            f.write('# Tcrit\tcharacteristic\tPkFlux\n')        ## INITIALISE LIST OF SIMILATED TRANSIENTS
                            write_source(params['INITIAL PARAMETERS']['file'] + '_' + lightcurvetype + '_SimTrans' , bursts) #file with starttime\tduration\tflux
                            print("Written Simulated Sources")
    
                t1 = datetime.datetime.now() 
                # det are the sources themselves while detbool is a numpy boolean array indexing all sources
                det, detbool = compute_lc.detect_bursts(obs[obsmask[:,i]], 
                    float(params['INITIAL PARAMETERS']['flux_err']), 
                    float(params['INITIAL PARAMETERS']['det_threshold']) , 
                    bursts, 
                    lightcurve.fluxint)
                t2 = datetime.datetime.now()
                dettime += (t2-t1).total_seconds()
                t1 = datetime.datetime.now() 
                stats[statcounter,0] = thisdur 
                stats[statcounter,1] = thisflux 
                stats[statcounter,2] = np.nan_to_num(np.sum(detbool)/targetnum) # probability for this bin
                statcounter+=1
                t2= datetime.datetime.now()
                stattime += (t2-t1).total_seconds()
              
                
        print(srcsimtime,"seconds simulating sources")
        print(dettime, "seconds detecting sources")
        print(stattime,"seconds aggregating stats")
        totaltime = srcsimtime + dettime +stattime
        print(100*srcsimtime/totaltime,"% of the time simulating sources")
        print(100*dettime/totaltime, "% of the time detecting sources")
        print(100*stattime/totaltime,"% of the time aggregating stats")
        if True: # not (np.isnan(burstlength) and np.isnan(burstflux)):
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
            tophatlc = tophat.tophat()
            fdstats = np.zeros(((len(flux_bins)-1)*(len(dur_ints)-1), 5), dtype=np.float32)
            # iterate over bins
            statcounter = 0
            for ldurbin,rdurbin in tqdm(zip(dur_ints[:-1],dur_ints[1:]),total=len(dur_ints[:-1])):
                thisdur = (ldurbin+rdurbin)/2
                for lfluxbin,rfluxbin in zip(flux_bins[:-1],flux_bins[1:]):
                    thisflux = (lfluxbin + rfluxbin)/2
                    fdbursts = compute_lc.generate_sources(targetnum, #n_sources
                        startepoch, #start_survey
                        stopepoch, #end_survey
                        lfluxbin, #Flux min
                        rfluxbin, #Flux max
                        ldurbin, # duration min
                        rdurbin,  #duration max
                        "tophat",
                        burstlength,
                        burstflux) # 
    
                    fdbursts['chartime'] += fake_obs['start'][0]
            
    
    
                    # det are the sources themselves while detbool is a numpy boolean array indexing all sources
                    fddet, fddetbool = compute_lc.detect_bursts(fake_obs[obsmask[:,i]], 
                        float(params['INITIAL PARAMETERS']['flux_err']), 
                        float(params['INITIAL PARAMETERS']['det_threshold']) , 
                        fdbursts, 
                        tophatlc.fluxint)
                    fdstats[statcounter,0] = thisdur 
                    fdstats[statcounter,1] = thisflux 
                    fdstats[statcounter,2] = np.nan_to_num(np.sum(fddetbool)/targetnum) # probability for this bin
                    statcounter+=1
                    

            fdlimit = np.where(fdstats[:,2] > 0.99 )
            for fdl in fdlimit:
                print(fdstats[fdl,0],fdstats[fdl,1])
        det_threshold = float(params['INITIAL PARAMETERS']['det_threshold'])
        extra_threshold = float(params['INITIAL PARAMETERS']['extra_threshold'])
        current_obs = obs[obsmask[:,i]]
    
            
        detections = int(params['INITIAL PARAMETERS']['detections'])
        confidence = float(params['INITIAL PARAMETERS']['confidence'])/100
        cdet = False 
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
        np.copy(stats),
        2,
        lightcurve.lines,
        regions['area'][i],
        tsurvey,
        detections,
        confidence,
        params['INITIAL PARAMETERS']['file'])
    if len(uniquepointFOV) > 1:
        combinedprobs = np.average([s[:,2] for s in statlist], weights = [np.full(statlist[0][:,2].shape, w) for w in regions['area']*np.array(tsurveylist)], axis=0)
        combinedname = regions['identity'][0]
        for i in range(1,len(uniquepointFOV)):
            combinedname += 'and'+str(regions['identity'][i])
        combinedstat = np.copy(stats)
        combinedstat[:,2] = combinedprobs
        compute_lc.make_mpl_plots(combinedname,
            fl_min,
            fl_max,
            bursts['chardur'].min(),
            bursts['chardur'].max(),
            det_threshold,
            extra_threshold,
            obs,
            False,
            params['INITIAL PARAMETERS']['file'],
            float(params['INITIAL PARAMETERS']['flux_err']),
            combinedstat,
            2,
            lightcurve.lines,
            np.sum(regions['area']),
            obs['start'][-1] + obs['duration'][-1] - obs['start'][0],
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
        if hasattr(lightcurve, 'docustompop'):
            print('Getting custom transient population for lightcurve')
            bursts = lightcurve.custompop(targetnum, #n_sources
                startepoch, #start_survey
                stopepoch, #end_survey
                float(params['INITIAL PARAMETERS']['fl_min']), #Flux min
                float(params['INITIAL PARAMETERS']['fl_max']), #Flux max
                float(params['INITIAL PARAMETERS']['dmin']), # duration min
                float(params['INITIAL PARAMETERS']['dmax']))  #duration max
            exit()
        else:
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
            bursts, 
            lightcurve.fluxint)
            
    
            
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
            fdbursts, 
            tophatlc.fluxint)
            
    
            
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
