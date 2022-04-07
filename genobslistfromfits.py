import numpy as np 
from astropy.io import fits
import glob
import argparse
import datetime
start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)
starttime = datetime.datetime.now().strftime('%y%m%d%Hh%Mm%Ss')
print("""This script assumes that:
1. All observations are the same duration
2. Noise can be well approximated by 1/sqrt(duration) scaling from a provided sample of noises
3. RA and DEC are in degrees in the header
4. There are no gaps in the observations""")

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

def getnoise(A, duration):
    rng = np.random.default_rng()
    sensitivity = A/np.sqrt(duration)
    return rng.normal(sensitivity, 0.08*sensitivity)


parser=argparse.ArgumentParser()
parser.add_argument("files", nargs='+', help="fits files to grab info from")
parser.add_argument("--duration", required=True, type=float, help="set singular duration for all observations in seconds (required)")
parser.add_argument("--noisefile", required=True, help="file containing a list of noises and durations to estimate noise in the format: (s, Jy) (required)")
parser.add_argument("--FOVrad", help='radius of the FOV of the images, otherwise will try to read from header')
args = parser.parse_args()
if len(args.files)==1:
     files = glob.glob(args.files)
else:
    files = args.files

samplenoises = np.loadtxt(args.noisefile, delimiter=',',dtype={'names': ('duration','sens'), 'formats': ('f8','f8')})
A = np.mean(samplenoises['sens']*np.sqrt(samplenoises['duration']))

startdates = []
ra =  []
dec = []
FOVrad = []
endtime = []
for myfile in files:
    with fits.open(myfile) as hdul:
        hdu = hdul[0]
        startdates.append(hdu.header['DATE-OBS']+'+00:00')
        ra.append(correctrarange(hdu.header['OBSRA']))
        dec.append(correctdecrange(hdu.header['OBSDEC']))
        if ("RA" in hdu.header['CTYPE1']) and ~bool(args.FOVrad):
            FOVrad.append(np.absolute(hdu.header['CDELT1']*hdu.header['CRPIX1']))
        else:
            FOVrad.append(args.FOVrad)
        filedatetime = datetime.datetime.strptime(hdu.header['DATE-OBS'],'%Y-%m-%dT%H:%M:%S.%f') 
        endfiledatetime = filedatetime + datetime.timedelta(seconds=args.duration)
        endtime.append(endfiledatetime.strftime('%Y-%m-%dT%H:%M:%S.%f+00:00'))

with open('obslist'+starttime+'.txt', 'a+') as f:
    for sd, ed, myra, mydec, FOV in zip(startdates, endtime, ra, dec, FOVrad):
            f.write("{},{},{},{},{},{},{}\n".format(sd, ed, getnoise(A, args.duration), myra, mydec, "False", FOV))       