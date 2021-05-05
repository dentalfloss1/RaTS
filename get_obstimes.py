import glob
import datetime
starttime = datetime.datetime.now().strftime('%y%m%d%Hh%Mm%Ss')
from casatools import msmetadata, ms 
import sys
# from astropy.coordinates import SkyCoord
import argparse
import os



def main_proc_loop(targetobs):
    def get_integrationtime(s):
        print(s)
        msobj = ms()
        msobj.open(targetobs)
        try:
            return msobj.getscansummary()[str(s)]['0']['IntegrationTime']
        except KeyError:
            return get_integrationtime(s+1)
    def correctrarange(ra):
        
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
        pointdict = msmd.phasecenter(field)
        print(pointdict['m0']['value'],pointdict['m1']['value'], pointdict['refer'],pointdict['m0']['unit'])
        if ('J2000' or 'fk5') not in pointdict['refer']:
            print('Unexpected coordinate system: please add conversion to J2000 in script')
        if 'rad' in pointdict['m0']['unit']:
            tmpra = pointdict['m0']['value']*180/np.pi
            tmpra = correctrarange(tmpra)
        elif 'deg' in pointdict['m0']['unit']:
            tmpra = pointdict['m0']['value']
            tmpra = correctrarange(tmpra)
        else:
            print("unknown units, convert to deg")
        if 'rad' in pointdict['m1']['unit']:
            tmpdec = pointdict['m1']['value']*180/np.pi
            print(tmpdec)
            tmpdec = correctdecrange(tmpdec)
        elif 'deg' in pointdict['m1']['unit']:
            tmpdec = pointdict['m1']['value']
            tmpdec = correctdecrange(tmpdec)
        else:
            print("unknown units, convert to deg")
        return tmpra, tmpdec
    try:
        # s here doesn't really matter much it just sets an index that will be incremented if the function fails, so start low
        integration_time = get_integrationtime(1)
        start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)

        msmd = msmetadata()
        msmd.open(msfile=targetobs)

        timelist = []
        for f in msmd.fieldsforintent('TARGET'):
            tmpra, tmpdec = get_radec(f)
            print(tmpra, tmpdec)
            scans = msmd.scansforfield(f)
            scantime = []
            
            for s in scans:
                timesinscan = msmd.timesforscan(s)
                scantime.append([start_epoch + datetime.timedelta(seconds=t) for t in timesinscan])
            
            beginepoch = min([min(time) for  time in scantime])
            endepoch = max([max(time) for time in scantime])
            timelist.append([(beginepoch - datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00'),
                        ((endepoch - beginepoch) + datetime.timedelta(seconds=round(integration_time)/2.0)).total_seconds(), tmpra, tmpdec])
            for s,t in zip(scans,scantime):

                timelist.append([(t[0] - datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00'),
                                ((t[-1] - t[0]) + datetime.timedelta(seconds=round(integration_time)/2.0)).total_seconds(), tmpra, tmpdec])
            for tint in msmd.timesforfield(f):
                timelist.append([((start_epoch + datetime.timedelta(seconds=tint)) - datetime.timedelta(seconds=round(integration_time)/2.0)).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00'),
                                round(integration_time), tmpra, tmpdec])
            
        with open('obslist'+starttime+'.txt', 'a+') as f:
            for t in timelist:
                    f.write("{},{},{},{}\n".format(t[0], t[1], t[2], t[3]))       
            f.flush()
    except RuntimeError:
        print("Issues with "+targetobs+". Skipping this one")

parser=argparse.ArgumentParser(
    description='''Crude script that pulls info from supplied ms or multiple ms using unix wildcards. Always check output. Requires casa6 and astropy.''',
    epilog="""Reads in input from a txt file containing ms folder names via the --obs flag""")
parser.add_argument("--obs", help="txt file containing a list of ms")

args = parser.parse_args()
# just gets a timestamp for a unique filename
import numpy as np 
observations = np.loadtxt(args.obs, dtype='<U100')
print(observations)
if observations.size<2:
    main_proc_loop(observations)  

else:
    for targetobs in observations:
        
        main_proc_loop(targetobs)
print('wrote to '+'obslist'+starttime+'.txt')