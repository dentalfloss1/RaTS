import glob
import datetime
# Used to create unique filenames
starttime = datetime.datetime.now().strftime('%y%m%d%Hh%Mm%Ss')
from casatools import msmetadata, ms 
import sys
import argparse
import os




def main_proc_loop(targetobs):
    """Takes as input a ms, outputs a modified list that will become the obs file"""
   
    def get_integrationtime(s):
        """Takes as input an index to start looking for the integration time. Recursively increases this index until 
        integration time is sucessfully retrieved from the ms without error and returns this time"""

        print(s) # just to let the user know something is happening
        msobj = ms()
        msobj.open(targetobs)
        try:
            return msobj.getscansummary()[str(s)]['0']['IntegrationTime']
        except KeyError:
            return get_integrationtime(s+1)
    def correctrarange(ra):
        """Takes an input ra in degrees and returns it on the interval between 0 and 360"""
        
        if (ra <= 360) and (ra >= 0):
            return ra
        elif (ra > 360):
            print("correcting ra angle range")
            ra= ra-360
            return correctrarange(ra)
        elif (ra < 360):
            print("correcting ra angle range")
            ra = ra + 360
            return correctrarange(ra)
    def correctdecrange(dec):
        """Takes an input dec in degrees and returns it on the interval between -90 and 90"""
        
        if (dec <= 90) and (dec >= -90):
            return dec
        elif (dec > 90):
            print("correcting dec angle range")
            dec= dec-90
            return correctdecrange(dec)
        elif (dec < -90):
            print("correcting dec angle range")
            dec = dec + 90
            return correctdecrange(dec)
    def get_radec(field):
        """Takes an input field and outputs ra and dec in degrees"""

        pointdict = msmd.phasecenter(field)
        print(pointdict['m0']['value'],pointdict['m1']['value'], pointdict['refer'],pointdict['m0']['unit'])
        # There's a lot of unknowns in how the fields are written 
        # Probably best to handle if it says 'rad' or 'deg' and otherwise return None 
        if ('J2000' or 'fk5') not in pointdict['refer']:
            print('Unexpected coordinate system: please add conversion to J2000 in script')
        if 'rad' in pointdict['m0']['unit']:
            tmpra = pointdict['m0']['value']*180/np.pi
            tmpra = correctrarange(tmpra)
        elif 'deg' in pointdict['m0']['unit']:
            tmpra = pointdict['m0']['value']
            tmpra = correctrarange(tmpra)
        # Add new elif statements here for any weird units
        else:
            print("unknown units, convert to deg")

        if 'rad' in pointdict['m1']['unit']:
            tmpdec = pointdict['m1']['value']*180/np.pi
            print(tmpdec)
            tmpdec = correctdecrange(tmpdec)
        elif 'deg' in pointdict['m1']['unit']:
            tmpdec = pointdict['m1']['value']
            tmpdec = correctdecrange(tmpdec)
        # Add new elif statements here for any weird units
        else:
            print("unknown units, convert to deg")
        return tmpra, tmpdec
    
    

    rng = np.random.default_rng()
    
    try:
        # s here doesn't really matter much it just sets an index that will be incremented if the function fails, so start low
        integration_time = get_integrationtime(1)
        # The epoch commonly used in ms 
        start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)
        def writescansfile(beginepoch,epochduration,scans, scantime,scansfile):
            """Takes as input info about a long observation with gaps, writes a file detailing the scans, and returns the total time in the gaps"""

            gaps = []
            gaptime = 0
            scanlist = []
            for i in range(len(scantime)-1):
                startgap = max(scantime[i])+datetime.timedelta(seconds=round(integration_time)/2.0)
                endgap = min(scantime[i+1])+datetime.timedelta(seconds=round(integration_time)/2.0)
                startscan = min(scantime[i])+datetime.timedelta(seconds=round(integration_time)/2.0)
                endscan = max(scantime[i+1])+datetime.timedelta(seconds=round(integration_time)/2.0)
                gaps.append([startgap,endgap])
                gaptime = gaptime + (min(scantime[i+1]) - max(scantime[i])).total_seconds()
            for s,t in zip(scans,scantime):
                scanend = (max(t)+datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                scanstart = (t[0]).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                scanduration = (max(t) - t[0]).total_seconds()
                sensitivity = constant/np.sqrt(scanduration)
                scanlist.append([scanstart, scanend,  rng.normal(sensitivity, 0.08*sensitivity), tmpra, tmpdec, False])
            with open(scansfile, "w") as f:
                for t in scanlist:
                    f.write("{},{},{},{},{},{},{}\n".format(t[0], t[1], t[2], t[3], t[4], t[5], fov))
            print("wrote to "+scansfile)

            return gaptime
            
        msmd = msmetadata()
        msmd.open(msfile=targetobs)
        # Yes, this just assumes that the intent is TARGET. Maybe that can change with optional user input, but for now this seems fine
        for f in msmd.fieldsforintent('TARGET'):
            tmpra, tmpdec = get_radec(f)
            print(tmpra, tmpdec)
            scans = msmd.scansforfield(f)
            fieldname = msmd.namesforfields(f)
            print(fieldname)
            scantime = []
            
            for s in scans:
                timesinscan = msmd.timesforscan(s)
                scantime.append([start_epoch + datetime.timedelta(seconds=t) for t in timesinscan])
            # This puts an entire ms as an entry in the obsfile
            # It needs to get the right date formats and also get the gapfile name and times 
            beginepoch = min([min(time) for  time in scantime]) - datetime.timedelta(seconds=round(integration_time)/2.0)
            mjdbeginepoch = (beginepoch-start_epoch).total_seconds()/60/60/24
            epochduration = max([max(time) for time in scantime]) - beginepoch
            endepoch = beginepoch + epochduration

            scansfile = "scans"+fieldname[0]+"-"+str(mjdbeginepoch)+"-"+str(epochduration.total_seconds())+".txt"
            gaptime = writescansfile(beginepoch,epochduration,scans, scantime,scansfile)
            total_intime = epochduration.total_seconds() - gaptime

            # constant was defined previously when we read in the sample noises 
            sensitivity = constant/np.sqrt(epochduration.total_seconds() - gaptime)
            
            timelist.append([beginepoch.strftime('%Y-%m-%dT%H:%M:%S.%f+00:00'),
                        endepoch.strftime('%Y-%m-%dT%H:%M:%S.%f+00:00'), rng.normal(sensitivity, 0.08*sensitivity), tmpra, tmpdec, scansfile])
            

            for s,t in zip(scans,scantime):
                scanstart = (t[0] - datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                scanduration = ((t[-1] - t[0]) + datetime.timedelta(seconds=round(integration_time)/2.0)).total_seconds()
                scanend = (t[0]).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                sensitivity = constant/np.sqrt(scanduration)
                timelist.append([scanstart, scanend,  rng.normal(sensitivity, 0.08*sensitivity), tmpra, tmpdec, False])
            for tint in msmd.timesforfield(f):
                intstart = ((start_epoch + datetime.timedelta(seconds=tint)) - datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                intduration = round(integration_time)
                intend = ((start_epoch + datetime.timedelta(seconds=intduration) + datetime.timedelta(seconds=tint)) - datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                sensitivity = constant/np.sqrt(intduration)
                timelist.append([intstart, intend,  rng.normal(sensitivity, 0.08*sensitivity), tmpra, tmpdec, False])
            
    except RuntimeError:
        # Hopefully this doesn't happen, but if it does skipping the entire ms is probably the best option
        print("Issues with "+targetobs+". Skipping this one")

parser=argparse.ArgumentParser(
    description='''Crude script that pulls info from supplied ms or multiple ms using unix wildcards. Always check output. Requires casa6 and astropy.''',
    epilog="""Reads in input from a txt file containing ms folder names via the --obs flag""")
parser.add_argument("--obs", help="txt file containing a list of ms", required=True)
parser.add_argument("--noise", help="txt file containing duration and sensitivity in Jy: duration (seconds),sensitivity (Jy)\n", required=True)
parser.add_argument("--fov", help="radius of the fov of instrument in degrees \n", type=float, required=True)

args = parser.parse_args()
import numpy as np 
fov = args.fov
observations = np.loadtxt(args.obs, dtype='<U100')
samplenoises = np.loadtxt(args.noise, delimiter=',',dtype={'names': ('duration','sens'), 'formats': ('f8','f8')})
constant = np.mean(samplenoises['sens']*np.sqrt(samplenoises['duration']))
print(observations)
timelist = []
if observations.size<2:
    main_proc_loop(observations)  
else:
    for targetobs in observations:
        main_proc_loop(targetobs)
with open('obslist'+starttime+'.txt', 'a+') as f:
    for t in timelist:
            f.write("{},{},{},{},{},{},{}\n".format(t[0], t[1], t[2], t[3], t[4], t[5], fov))       
    f.flush()
print('wrote to '+'obslist'+starttime+'.txt')
