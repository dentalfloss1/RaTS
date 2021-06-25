import configparser
import datetime
import numpy as np
import os
import glob
import argparse
import warnings
from tqdm import tqdm
from astropy import units as u 
from astropy.coordinates import SkyCoord, CartesianRepresentation
from scipy.special import binom
import scipy.interpolate as interpolate
from bitarray import bitarray
import matplotlib.pyplot as plt
from matplotlib import ticker, colors

def observing_strategy(obs_setup, det_threshold, nobs, obssens, obssig, obsinterval, obsdurations):
    """Parse observation file or set up trial mode. Return array of observation info and a regions observed"""

    rng = np.random.default_rng()
    start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)
    if obs_setup is not None:
        tstart, tdur, sens, ra, dec, gapsfile,fov  = np.loadtxt(obs_setup, unpack=True, delimiter = ',',
            dtype={'names': ('start', 'duration','sens', 'ra', 'dec','gaps', 'fov'), 'formats': ('U32','f8','f8','f8','f8','<U128','f8')})
        tstart = np.array([(datetime.datetime.strptime(t, "%Y-%m-%dT%H:%M:%S.%f+00:00") - start_epoch).total_seconds()/3600/24 for t in tstart])
        tdur = np.array([datetime.timedelta(seconds=t).total_seconds()/3600/24 for t in tdur])
        sortkey = np.argsort(tstart)
        obs = np.zeros(len(tstart),
              dtype={'names': ('start', 'duration','sens','gaps'), 'formats': ('f8','f8','f8','<U128')})
        obs['start']=tstart[sortkey]
        obs['duration']=tdur[sortkey]
        obs['sens']+=sens[sortkey]*det_threshold
        obs['gaps']=gapsfile[sortkey]
        pointing = np.array([np.array([r,d,f]) for r,d,f in zip(ra[sortkey],dec[sortkey],fov[sortkey])])
        pointFOV = pointing
    else: # Enter 'trial mode' according to specified cadence
        observations = np.zeros(nobs,dtype={'names': ('start', 'duration','sens','gaps'), 'formats': ('f8','f8','f8','<U128')})
        for i in range(nobs):
            tstart = (datetime.datetime.strptime("2019-08-08T12:50:05.0", "%Y-%m-%dT%H:%M:%S.%f") + datetime.timedelta(days=i*7) - start_epoch).total_seconds()/3600/24#) + 60*60*24*obsinterval*i)/(3600*24)    #in days at an interval of 7 days
            tdur = datetime.timedelta(days=obsdurations).total_seconds()/3600/24
            sens = rng.normal(obssens, obssig) * det_threshold # in Jy. Choice of Gaussian and its characteristics were arbitrary.
            observations['start'][i] = tstart
            observations['duration'][i] = tdur
            observations['sens'][i] = sens
            observations['gaps'][i] = "False"
        # obs = np.zeros(len(observations)
        obs = observations
            
        pointing = np.array([np.array([275.0913169, 7.185135679]) for l in observations])
        # pointing = np.array([np.array([342.7528844,-59.12614311]) for l in observations])
        tmpsc = SkyCoord(ra=275.0913169*u.degree,dec=7.185135679*u.degree,frame='fk5')
        tmpscoff1 = tmpsc.directional_offset_by(-30.00075*u.degree, (1/np.sqrt(2))*u.degree)
        tmpscoff2 = tmpsc.directional_offset_by(30.00075*u.degree, (1/np.sqrt(2))*u.degree)
        # # print()
        pointing[:15]-=[tmpsc.ra.deg - tmpscoff1.ra.deg,tmpsc.dec.deg - tmpscoff1.dec.deg]
        pointing[15:30]-=[tmpsc.ra.deg - tmpscoff2.ra.deg,tmpsc.dec.deg - tmpscoff2.dec.deg]
        # pointing[::3]=[342.5811999,-59.12369356]
        # pointing[1::3]=[342.6688451,-59.04494042]
        # print(pointing)
        # pointing[2::3]-=[-1,1]
        pointing = pointing[obs['start'].argsort()]
        obs = obs[obs['start'].argsort()] # sorts observations by date
        FOV = np.array([1.4 for l in observations]) # make FOV for all observations whatever specified here, 1.5 degrees for example
        # FOV = np.array([0.059505254 for l in observations]) # make FOV for all observations whatever specified here, 1.5 degrees for example
        # FOV[::3]=0.06105040575
        # FOV[1::3]=0.0636070677
        # pointing[0:6,0]+=2 # Make an offset between some pointings
        pointFOV = np.zeros((len(obs),3))
        pointFOV[:,0:2] += pointing
        pointFOV[:,2] += FOV
    return obs, pointFOV
        
def calculate_regions(pointFOV, observations):
    """Calculate regions based on simultaneous observing times assuming circular regions. Returns region info as structured numpy array."""
    
    uniquepoint = np.unique(pointFOV,axis=0)
    uniquesky = SkyCoord(ra=uniquepoint[:,0],dec=uniquepoint[:,1], unit='deg', frame='fk5')
    ### Next two lines set up variables for defining region properties. Region array is of maximum theoretical length assuming no more than double overlapping fovs (triple or more never get computed ##
    numrgns = len(uniquepoint) 
    regions = np.zeros(np.uint32(numrgns + binom(numrgns,2) + binom(numrgns,3)), dtype={'names': ('ra', 'dec','identity', 'area', 'timespan', 'stop', 'start'), 'formats': ('f8','f8','U32','f8', 'f8', 'f8', 'f8')})
    for i in range(numrgns): # Label the individual pointings. These regions are for example region 1 NOT 2 and 2 NOT 1 
        regions['identity'][i] = str(i)
        regions['ra'][i] = uniquepoint[i,0]
        regions['dec'][i] = uniquepoint[i,1]
        regions['area'][i] = (4*np.pi*np.sin(uniquepoint[i,2]*(np.pi/180/2))**2*(180/np.pi)**2) # Assumes single circular regions, for multiple pointings or other shapes this needs altering
        leftoff = i + 1
    obssubsection = []
    
    for p in uniquepoint:
        timeind = np.array([np.amax(np.argwhere((pointFOV[:,0] == p[0]) & (pointFOV[:,1] ==p[1]))), np.amin(np.argwhere((pointFOV[:,0] == p[0]) & (pointFOV[:,1] ==p[1])))])
        matchregion = (regions['ra']==p[0]) & (regions['dec']==p[1]) & (regions['area']==(4*np.pi*np.sin(p[2]*(np.pi/180/2))**2*(180/np.pi)**2))
        regions['timespan'][matchregion] = (observations['start'][timeind[0]] + observations['duration'][timeind[0]] - observations['start'][timeind[1]])
        regions['stop'][matchregion] = observations['start'][timeind[0]] + observations['duration'][timeind[0]]
        regions['start'][matchregion] = observations['start'][timeind[1]]
        obssubsection.append([timeind[1],timeind[0],regions['identity'][matchregion][0]])

    gamma = np.zeros((len(uniquepoint),len(uniquepoint)))
    for i in range(len(uniquesky)-1): # Label intersections: For example: 1 AND 2
        for j in range(i+1,len(uniquesky)):
            if uniquesky[i].separation(uniquesky[j]).deg <= (uniquepoint[i,2] + uniquepoint[j,2]):
                d = uniquesky[i].separation(uniquesky[j]).rad
                r1 = uniquepoint[i,2]*np.pi/180
                r2 = uniquepoint[j,2]*np.pi/180
                gamma[i,j] = np.arctan((np.cos(r2)/np.cos(r1)/np.sin(d)) - (1/np.tan(d)))
                pa = uniquesky[i].position_angle(uniquesky[j])
                # https://arxiv.org/ftp/arxiv/papers/1205/1205.1396.pdf
                # and https://en.wikipedia.org/wiki/Solid_angle#Cone,_spherical_cap,_hemisphere
                fullcone1 = 4*np.pi*np.sin(r1/2)**2 
                cutchord1 = 2*(np.arccos(np.sin(gamma[i,j])/np.sin(r1)) - np.cos(r1)*np.arccos(np.tan(gamma[i,j])/np.tan(r1))) 
                fullcone2 = 4*np.pi*np.sin(r2/2)**2
                cutchord2 = 2*(np.arccos(np.sin(gamma[i,j])/np.sin(r2)) - np.cos(r2)*np.arccos(np.tan(gamma[i,j])/np.tan(r2))) 
                centerreg = uniquesky[i].directional_offset_by(pa, gamma[i,j]*u.radian)
                regions['identity'][leftoff] = str(i)+'&'+str(j)
                regions['ra'][leftoff] = centerreg.ra.deg
                regions['dec'][leftoff] = centerreg.dec.deg
                regions['area'][leftoff] = (cutchord1 + cutchord2)*(180/np.pi)**2
                regions['start'][leftoff] = min(regions['start'][i],regions['start'][j])
                regions['stop'][leftoff] = max(regions['stop'][i],regions['stop'][j])
                regions['timespan'][leftoff] = regions['stop'][leftoff] - regions['start'][leftoff]
                leftoff+=1
    # scatterpointsra = []
    # scatterpointsdec = []
    for i in range(len(uniquesky)-2): # repeat the above, but this time for triple overlapping regions
        for j in range(i+1,len(uniquesky)-1):
            for index3 in range(j+1,len(uniquesky)):
                if ((uniquesky[i].separation(uniquesky[j]).deg <= (uniquepoint[i,2] + uniquepoint[j,2])) and 
                    (uniquesky[j].separation(uniquesky[index3]).deg <= (uniquepoint[j,2] + uniquepoint[index3,2])) and 
                    (uniquesky[i].separation(uniquesky[index3]).deg <= (uniquepoint[i,2] + uniquepoint[index3,2]))):
                    r1 = uniquepoint[i,2]*np.pi/180
                    r2 = uniquepoint[j,2]*np.pi/180
                    r3 = uniquepoint[index3,2]*np.pi/180
                    # Get coordinates of the encircled(?) spherical triangle
                    # from the triangle formed between pointing center, overlap center, and overlap nodal point.

                    angle_offset = 90*u.deg
                    halfheightr4 = np.arccos(np.cos(r2)/np.cos(gamma[i][j])) 
                    point4key = np.where(regions['identity'] == str(i)+'&'+str(j))
                    point4ra = regions['ra'][point4key]
                    point4dec = regions['dec'][point4key]
                    point4sc = SkyCoord(ra=point4ra, dec=point4dec, unit='deg',frame='fk5')
                    point4pa = point4sc.position_angle(uniquesky[j])
                    point7sc = point4sc.directional_offset_by(point4pa + angle_offset, halfheightr4)
                    # print(point7sc.separation(uniquesky[index3]).deg)
                    if point7sc.separation(uniquesky[index3]).deg > uniquepoint[index3,2]:
                        point7sc = point4sc.directional_offset_by(point4pa - angle_offset, halfheightr4)

                    halfheightr5 = np.arccos(np.cos(r3)/np.cos(gamma[j][index3]))
                    point5key = np.where(regions['identity'] == str(j)+'&'+str(index3))
                    point5ra = regions['ra'][point5key]
                    point5dec = regions['dec'][point5key]
                    point5sc = SkyCoord(ra=point5ra, dec=point5dec, unit='deg',frame='fk5')
                    point5pa = point5sc.position_angle(uniquesky[index3])
                    point8sc = point5sc.directional_offset_by(point5pa + angle_offset, halfheightr5)
                    # print(point8sc.separation(uniquesky[i]).deg)
                    if point8sc.separation(uniquesky[i]).deg > uniquepoint[i,2]:
                        point8sc = point5sc.directional_offset_by(point5pa - angle_offset, halfheightr5)

                    halfheightr6 = np.arccos(np.cos(r1)/np.cos(gamma[i][index3]))
                    point6key = np.where(regions['identity'] == str(i)+'&'+str(index3))
                    point6ra = regions['ra'][point6key]
                    point6dec = regions['dec'][point6key]
                    point6sc = SkyCoord(ra=point6ra, dec=point6dec, unit='deg',frame='fk5')
                    point6pa = point6sc.position_angle(uniquesky[i])
                    point9sc = point6sc.directional_offset_by(point6pa + angle_offset, halfheightr6)
                    # print(point9sc.separation(uniquesky[j]).deg)
                    if point9sc.separation(uniquesky[j]).deg > uniquepoint[j,2]:
                        point9sc = point6sc.directional_offset_by(point6pa - angle_offset, halfheightr6)

                    #We now get the side lengths of the encircled triangle from the coordinates. 
                    aside = point7sc.separation(point8sc).rad[0]
                    bside = point8sc.separation(point9sc).rad[0]
                    cside = point9sc.separation(point7sc).rad[0]

                    # spherical law of cosines

                    Aangle = np.arccos((np.cos(aside) - np.cos(bside)*np.cos(cside))/(np.sin(bside)*np.sin(cside)))
                    Bangle = np.arccos((np.cos(bside) - np.cos(cside)*np.cos(aside))/(np.sin(cside)*np.sin(aside)))
                    Cangle = np.arccos((np.cos(cside) - np.cos(aside)*np.cos(bside))/(np.sin(aside)*np.sin(bside)))

                    triarea = Aangle + Bangle + Cangle - np.pi

                    # We now need to get the excess area from the overlapping region not actually being a spherical triange. 
                    # we will use the triangle formed from a pointing center, a point of the overlap region, and the midpoint of
                    # the encircled circular triangle

                    gamma4 = np.arccos(np.cos(r2)/np.cos(aside/2))
                    cutchord1 = 2*(np.arccos(np.sin(gamma4)/np.sin(r2)) - np.cos(r2)*np.arccos(np.tan(gamma4)/np.tan(r2))) 
                    gamma5 = np.arccos(np.cos(r3)/np.cos(bside/2))
                    cutchord2 = 2*(np.arccos(np.sin(gamma5)/np.sin(r3)) - np.cos(r3)*np.arccos(np.tan(gamma5)/np.tan(r3))) 
                    gamma6 = np.arccos(np.cos(r1)/np.cos(cside/2))
                    cutchord3 = 2*(np.arccos(np.sin(gamma6)/np.sin(r1)) - np.cos(r1)*np.arccos(np.tan(gamma6)/np.tan(r1))) 

                    area = triarea + cutchord1 + cutchord2 + cutchord3


                    midpointra = (point7sc.ra.deg + point8sc.ra.deg + point9sc.ra.deg)/3
                    midpointdec = (point7sc.dec.deg + point8sc.dec.deg + point9sc.dec.deg)/3
                    # print(point7sc, point8sc, point9sc)
                    regions['identity'][leftoff] = str(i)+'&'+str(j)+'&'+str(index3)
                    regions['ra'][leftoff] = midpointra
                    regions['dec'][leftoff] = midpointdec
                    regions['area'][leftoff] = area*(180/np.pi)**2
                    regions['start'][leftoff] = np.amin([regions['start'][i],regions['start'][j],regions['start'][index3]])
                    regions['stop'][leftoff] = np.amax([regions['stop'][i],regions['stop'][j],regions['stop'][index3]])
                    regions['timespan'][leftoff] = regions['stop'][leftoff] - regions['start'][leftoff]
                    leftoff+=1

                    # scatterpointsra.extend([point7sc.ra,point8sc.ra,point9sc.ra])
                    # scatterpointsdec.extend([point7sc.dec,point8sc.dec,point9sc.dec])
    # from astropy.wcs import WCS
    # from astropy.io import fits
    # from astropy.utils.data import get_pkg_data_filename
    # from astropy.visualization.wcsaxes import SphericalCircle
    # # filename = get_pkg_data_filename('allsky/allsky_rosat.fits') #'/media/sarah/Elements/GRB200219A/1582287955_sdp_l0.GRB200219A_im_3.fits') #'E:\\GRB200219A\\1582287955_sdp_l0.GRB200219A_im_3.fits')
    # import matplotlib.pyplot as plt
    # ax = plt.subplot()
    # # ax.imshow(hdu.data, vmin=1.3, vmax=900, origin='lower')
    # for p in uniquepoint:
    #     ax.add_patch(SphericalCircle((p[0] * u.deg, p[1] * u.deg), p[2] * u.degree,
    #                     edgecolor='red', facecolor='none'))#,
    #                     #  transform=ax.get_transform('galactic')))
    # print("vertices: ", scatterpointsra, scatterpointsdec)
    # for i in range(len(scatterpointsra)):
    #     print(scatterpointsra[i],scatterpointsdec[i], ['s','P','X'][i])
    #     ax.scatter(scatterpointsra[i],scatterpointsdec[i], marker=['s','P','X'][i],s=100, c='black')
    # plt.show()
    # plt.close()
    print(regions)
    
    return regions[regions['identity'] != ''], obssubsection
    

def generate_pointings(n_sources, pointFOV, i, leftoff, overlapnums):
    """Simulate pointings for each simulated source. Use a monte-carlo like method to roll the dice to determine position. Return a numpy array and bitarray"""
    
    uniquepointFOV = np.unique(pointFOV, axis=0)
    rng = np.random.default_rng() 
    maxFOV = np.max(pointFOV[:,2])
    minra = min(uniquepointFOV[:,0])
    mindec = min(uniquepointFOV[:,1])
    maxra = max(uniquepointFOV[:,0])
    maxdec = max(uniquepointFOV[:,1])
    uniqueskycoord = SkyCoord(ra=uniquepointFOV[:,0]*u.deg, dec=uniquepointFOV[:,1]*u.deg, frame='fk5')
    ra_randfield = lambda n, ptfov: (rng.random(int(n))*(min(ptfov[0] + ptfov[2],360) - max(ptfov[0] - ptfov[2],0)) + max(ptfov[0] - ptfov[2],0)) * u.degree
    dec_randfield = lambda n, ptfov: (rng.random(int(n))*(min(ptfov[1] + ptfov[2],90)  - max(ptfov[1] - ptfov[2], -90)) + max(ptfov[1] - ptfov[2], -90)) * u.degree    
    rollforskycoordfield = lambda n, ptfov: SkyCoord(ra=ra_randfield(n, ptfov), dec=dec_randfield(n, ptfov), frame='fk5')
    sc = rollforskycoordfield(n_sources, uniquepointFOV[i])
    accept = np.zeros(sc.shape,dtype=bool)
    reject = np.logical_not(accept)
    cond = np.zeros(len(sc),dtype=bool)
    while True:
        cond[reject] +=((sc[reject].separation(uniqueskycoord[i]).deg) < (uniquepointFOV[i,2])) 
        if np.sum(~cond)==0:
            break
        reject = ~cond
        accept = ~reject
        print("Rejected sources: ",np.sum(reject))
        sc.data.lon[reject] = ra_randfield(len(sc[reject]),uniquepointFOV[i])
        sc.data.lat[reject] = dec_randfield(len(sc[reject]),uniquepointFOV[i])
        sc.cache.clear()

    numrgns = len(uniqueskycoord)
    print("Aggregating numbers of sources in overlap regions")
    for j in range(i+1,len(uniqueskycoord)):
        if (uniqueskycoord[i].separation(uniqueskycoord[j]).deg <= (uniquepointFOV[i,2] + uniquepointFOV[j,2])) & (i!=j):
            overlapnums['name'][leftoff - numrgns] = str(i)+'&'+str(j)
            overlapnums['sources'][leftoff - numrgns] =  np.sum((sc.separation(uniqueskycoord[i]).deg < (uniquepointFOV[i,2])) & (sc.separation(uniqueskycoord[j]).deg < (uniquepointFOV[j,2])))
            leftoff+=1
    for j in range(i+1,len(uniqueskycoord)):
        for k in range(j+1,len(uniqueskycoord)):
            if ((uniqueskycoord[i].separation(uniqueskycoord[j]).deg <= (uniquepointFOV[i,2] + uniquepointFOV[j,2])) & 
                (uniqueskycoord[j].separation(uniqueskycoord[k]).deg <= (uniquepointFOV[j,2] + uniquepointFOV[k,2])) &
                (uniqueskycoord[i].separation(uniqueskycoord[k]).deg <= (uniquepointFOV[i,2] + uniquepointFOV[k,2])) &
                (i!=j) & (j!=k)):
                overlapnums['name'][leftoff - numrgns] = str(i)+'&'+str(j)+'&'+str(k)
                overlapnums['sources'][leftoff - numrgns] =  np.sum((sc.separation(uniqueskycoord[i]).deg < (uniquepointFOV[i,2])) & 
                    (sc.separation(uniqueskycoord[j]).deg < (uniquepointFOV[j,2])) &
                    (sc.separation(uniqueskycoord[k]).deg < (uniquepointFOV[k,2])))
                leftoff+=1
    #### Put another loop here, but nested and use it to determine number of sources in triple overlap region

    return overlapnums,leftoff
        # print(i*n_sources + btarr[0:i*n_sources].count(bitarray('1')), (i*n_sources + btarr[0:i*n_sources].count(bitarray('1')) + int(targetnum)))
        # btarr[(i*n_sources + btarr[0:i*n_sources].count(bitarray('1'))):(i*n_sources + btarr[0:i*n_sources].count(bitarray('1')) + int(targetnum))] = True
    # print(btarr[0:n_sources].count(bitarray('1')))
    # return btarr

def generate_sources(n_sources, start_survey, end_survey, fl_min, fl_max, dmin, dmax, lightcurve, burstlength, burstflux):
    """Generate characteristic fluxes and characteristic durations for simulated sources. Return as numpy array"""
    
    start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)
    rng = np.random.default_rng()
    bursts = np.zeros(n_sources, dtype={'names': ('chartime', 'chardur','charflux'), 'formats': ('f8','f8','f8')}) # initialise the numpy array
    
    if not np.isnan(burstlength):
        bursts['chardur'] += burstlength
        print("setting burst length to single value")
        print(np.where(bursts['chardur']!=burstlength))
    else:
        bursts['chardur'] = 10**(rng.random(n_sources)*(np.log10(dmax) - np.log10(dmin)) + np.log10(dmin)) # random number for duration
        # bursts['chardur'] = (rng.random(n_sources)*(dmax - dmin) + dmin) # random number for duration
    # bursts['chardur'] += 500*7 + 0.01
    if not np.isnan(burstflux):
        bursts['charflux'] += burstflux
        print("setting burst flux to single value")
        print(np.where(bursts['charflux']!=burstflux))
    else:
        bursts['charflux'] = 10**(rng.random(n_sources)*(np.log10(fl_max) - np.log10(fl_min)) + np.log10(fl_min)) # random number for flux
    # bursts['charflux'] = (rng.random(n_sources)*(fl_max - fl_min) + fl_min) # random number for flux
    return bursts
    
def generate_start(bursts, potential_start, potential_end, n_sources):
    """Generate characteristic times to go with durations and fluxes. Return modified numpy array"""
    
    print("Generating starts to simulated sources")
    rng = np.random.default_rng() 
    bursts['chartime'] = rng.random(n_sources)*(potential_end - potential_start) + potential_start
    bursts.sort(axis=0,order='chartime')
    return bursts
    

def detect_bursts(obs, flux_err,  det_threshold, extra_threshold, sources, gaussiancutoff, edges, fluxint, file, dump_intermediate, write_source, pointFOV):
    """Detect simulated sources by using a series of conditionals along with the integrated flux calculation. Returns detected sources and the boolean array to get source indices"""
    
    start = datetime.datetime.now()
    uniquepointFOV = np.unique(pointFOV,axis=0)
    rng = np.random.default_rng() 
    # This bitarray is defined here because it tracks overall detections through all observations

    
    if edges[0] == 1 and edges[1] == 1: # TWO edges
        edgecondmask = lambda start_obs, end_obs: (sources['chartime'] + sources['chardur'] > start_obs) & (sources['chartime'] < end_obs) #
    elif edges[0] == 0 and edges[1] == 1: # Only ending edge, In this case, critical time is the end time
        edgecondmask = lambda start_obs, end_obs: (sources['chartime'] > start_obs)
    elif edges[0] == 1 and edges[1] == 0: # Only starting edge 
        edgecondmask = lambda start_obs, end_obs: (sources['chartime'] < end_obs)
    elif edges[0] == 0 and edges[1] == 0:
        edgecondmask = lambda start_obs, end_obs: (sources['chartime'] == sources['chartime'])
    candbitarr = bitarray(len(sources)*len(obs))
    print("Enforcing edge conditions in detection loop")
    for i in tqdm(range(len(obs))): # bitarray stores whether or not sources fall within observations
        candbitarr[i*len(sources):(i+1)*len(sources)] = bitarray(list(edgecondmask(obs['start'][i],obs['start'][i] + obs['duration'][i]))) 
 
    end = datetime.datetime.now()
    
    print('conditionals elapsed: ',end-start)
    start = datetime.datetime.now()
    candidates = bitarray(len(sources)) # True if source meets candidate criteria including det_threshold
    candidates.setall(False)
    extra_candidates = bitarray(len(sources)) # True if source meets candidate criteria plus extra threshold
    extra_candidates.setall(False)
    detallbtarr = bitarray(len(sources))
    detallbtarr.setall(True) # Bitarray that determines if source is detected in every observation. If it is, we set it to "not detected" since it is constant.
    print('Enforcing flux conditions in detection loop')
    for i in tqdm(range(len(obs))):
        if obs['gaps'][i] != 'False':
            subobs, _ = observing_strategy(obs['gaps'][i], det_threshold, 1, 1, 1, 1, 1) # We are giving the scansfile name, so the other variables are unimportant, we set them to 1 
            flux_int = np.zeros((len(sources)),dtype=np.float32)
            candind = np.array(candbitarr[i*len(sources):(i+1)*len(sources)].search(bitarray([True]))) # Turn the candbitarr into indices. Clunky, but it's the best way to do it I think.
            if candind.size == 0: # No candidates!
                detallbtarr.setall(False) # Otherwise will reject all previous candidates
                break
            F0_o = sources['charflux'][candind]
            tcrit = sources['chartime'][candind]
            tau = sources['chardur'][candind]
            subfluxint = np.zeros([len(subobs),len(candind)], dtype=np.float32)
            for j in range(len(subobs)):
                end_subobs = subobs['start'][j]+subobs['duration'][j]
                start_subobs = subobs['start'][j]
                error = np.sqrt((F0_o * flux_err)**2 + (obs['sens'][i]/det_threshold)**2) 
                F0 =rng.normal(F0_o, error)
                F0[(F0<0)] = F0[(F0<0)]*0
                subfluxint[j] = fluxint(F0, tcrit, tau, end_subobs, start_subobs)
                subfluxint[j][subfluxint[j] < 0] = 0
                # print(len(subfluxint > subobs['sens'][j]))
            flux_int[candind] = np.average(subfluxint,weights=subobs['duration']/np.sum(subobs['duration']),axis=0)               
            flux_int[candind][flux_int[candind] < 0] = 0 

        else:
            flux_int = np.zeros((len(sources)),dtype=np.float32)
            candbitarr[i*len(sources):(i+1)*len(sources)].setall(True)
            candind = np.array(candbitarr[i*len(sources):(i+1)*len(sources)].search(bitarray([True]))) # Turn the candbitarr into indices. Clunky, but it's the best way to do it I think.
            if candind.size == 0: # No candidates!
                detallbtarr.setall(False) # Otherwise will reject all previous candidates
                break
            F0_o = sources['charflux'][candind]
            tcrit = sources['chartime'][candind]
            tau = sources['chardur'][candind]
            end_obs = obs['start'][i]+obs['duration'][i]
            start_obs = obs['start'][i]
    ######## Convert from random.gauss to numpy verion  ######
    #################### ORIGINAL VERSION had an extra rng normal to sim instrument noise? #################
            # error = np.sqrt((abs(rng.normal(F0_o * flux_err, 0.05 * F0_o * flux_err)))**2 + (obs['sens'][i]/det_threshold)**2) 
            error = np.sqrt((F0_o * flux_err)**2 + (obs['sens'][i]/det_threshold)**2) 
            # error = F0_o*flux_err
            F0 =rng.normal(F0_o, error)
            # F0=F0_o
            F0[(F0<0)] = F0[(F0<0)]*0
            flux_int[candind] = fluxint(F0, tcrit, tau, end_obs, start_obs) # uses whatever class of lightcurve supplied: tophat, ered, etc      
            
        sensitivity = obs['sens'][i]
        candidates |= bitarray(list(flux_int > obs['sens'][i])) # Do a bitwise or to determine if candidates meet flux criteria. Sources rejected by edge criteria above are at zero flux anyway
        extra_sensitivity =  obs['sens'][i] * (det_threshold + extra_threshold) / det_threshold
        extra_candidates |= bitarray(list(flux_int > extra_sensitivity))
        candidates &= extra_candidates # Weird, I know. We really just use the extra detection criteria, but this is a holdover
        best_extra_sensitivity = np.amin(obs['sens']) * (det_threshold + extra_threshold) / det_threshold
        detallbtarr &= (bitarray(list(flux_int > extra_sensitivity)) | bitarray(list(flux_int > best_extra_sensitivity))) # End result is true if all are true. In other words, detected in every obs.
    candidates &= ~detallbtarr # We do a bitwise not and the result is that we only keep candidates that were not detected in every obs
    end = datetime.datetime.now()    
    print('indexing elapsed: ',end-start)
    detections = np.zeros(len(sources), dtype=bool)
    detections[candidates.search(bitarray([True]))] = True # Again, this is how we turn a bitarray into indices

    return sources[detections], detections

def statistics(fl_min, fl_max, dmin, dmax, det, all_simulated):
    """Calculate probabilities based on detections vs simulated, return a numpy array"""
    
    flux_bins = np.geomspace(fl_min, fl_max, num=int(round((np.log10(fl_max)-np.log10(fl_min))/0.05)), endpoint=True)
    dur_ints = np.geomspace(dmin, dmax, num=int(round((np.log10(dmax)-np.log10(dmin))/0.05)), endpoint=True)

    fluxes = np.array([],dtype=np.float32)
    durations = np.array([],dtype=np.float32)
    

    stats = np.zeros(((len(flux_bins)-1)*(len(dur_ints)-1), 5), dtype=np.float32)

    allhistarr, _, _ = np.histogram2d(all_simulated['chardur'], all_simulated['charflux'],bins = [dur_ints,flux_bins])
    dethistarr, _, _ = np.histogram2d(det['chardur'], det['charflux'],bins = [dur_ints,flux_bins])
    # The following ignores divide by zero and 0/0 errors and replaces nans with the smallest representable number
    # and infinities with the largest representable number. This is fine here because the probability is bounded by zero and 1 
    # so it won't explode out of control. 
    with np.errstate(divide='ignore', invalid='ignore'):
        probabilities = np.nan_to_num(dethistarr/allhistarr)
       
    durations = (dur_ints[:-1] + dur_ints[1:])/2
    fluxes = (flux_bins[:-1] + flux_bins[1:])/2   

    stats[:,0] = np.repeat(durations, len(fluxes))
    stats[:,1] = np.tile(fluxes, len(durations))
    stats[:,2] = probabilities.flatten()
    stats[:,3] = dethistarr.flatten()
    stats[:,4] = allhistarr.flatten()

    return stats

def make_mpl_plots(rgn, fl_min,fl_max,dmin,dmax,det_threshold,extra_threshold,obs,cdet,file,flux_err,toplot,gaussiancutoff,lclines,area,tsurvey,detections,confidence):
    """Use Matplotlib to make plots and if that fails dump numpy arrays. Returns an int that indicates plotting success or failure"""

    # Make histograms of observation senstivities and false detections
    fluxbins = np.geomspace(fl_min, fl_max, num=int(round((np.log10(fl_max)-np.log10(fl_min))/0.01)), endpoint=True)
    fddethist, fddetbins = np.histogram(cdet, bins=fluxbins, density=False)
    senshist, sensbins = np.histogram(obs['sens'], bins=fluxbins, density=False)
    # Calculate 99% of false detections line (from the left)
    histsum = 0
    for j in range(len(fddethist)):
        histsum += fddethist[j]
        try:
            if (histsum)/np.sum(fddethist) >= 0.99:
                print("stopped on index:",j,"exact percentage threshold:", (histsum)/np.sum(fddethist))
                vlinex = np.full(10, fluxbins[j])
                vliney = np.linspace(1, np.amax([fddethist,senshist]), num=10)
                break
        except RuntimeWarning:
            vlinex = np.full(10,fl_min)
            vliney = np.linspace(1,np.amax([fddethist,senshist]), num=10)

    # prepare variables for making surface plots
    # I don't like the way this np.log10 works, but much frustration dictates that I don't touch this because it's 
    # difficult to make work properly. 
    toplot[:,0] = np.log10(toplot[:,0])
    toplot[:,1] = np.log10(toplot[:,1])
    # An array of the gaps in the observations. This calculation is repeated for some light curves for no real reason
    gaps = np.array([],dtype=np.float32)
    for i in range(len(obs)-1):
        gaps = np.append(gaps, (obs['start'][i+1] - obs['start'][i] + obs['duration'][i]))
        # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])

    min_sens = min(obs['sens'])
    max_sens = max(obs['sens'])
    extra_thresh = max_sens / det_threshold * (extra_threshold + det_threshold)
    sens_last = obs['sens'][-1]
    sens_maxgap = obs['sens'][np.where((gaps[:] == max(gaps)))[0]+1][0]

    durmax = (obs['start'][-1] + obs['duration'][-1] - obs['start'][0])
    day1_obs = obs['duration'][0]
    max_distance = max(gaps)

    dmin=min(toplot[:,0])
    dmax=max(toplot[:,0])
    flmin=min(toplot[:,1])
    flmax=max(toplot[:,1])

    xs = np.arange(dmin, dmax, 1e-3)
    ys = np.arange(flmin, flmax, 1e-3)

    xs = xs[0:-1]
    day1_obs_x = np.empty(len(ys)) # I think these three arrays are now no longer needed, but I am hesitant to delete
    day1_obs_x.fill(day1_obs)
    
    sensmin_y = np.empty(len(xs))
    sensmin_y.fill(min_sens)
    
    sensmax_y = np.empty(len(xs))
    sensmax_y.fill(max_sens)

    X = np.linspace(dmin, dmax, num = 1001)
    Y = np.linspace(flmin, flmax, num = 1001)
    X = (X[0:-1] + X[1:])/2 # place interpolated Z values halfway between bin edges
    Y = (Y[0:-1] + Y[1:])/2

    X, Y = np.meshgrid(X, Y)
    
    Z = interpolate.griddata(toplot[:,0:2], toplot[:,2], (X, Y), method='linear')

    # do calculations for transient rate plot 
    durations = toplot[:,0]
    fluxes = toplot[:,1]
    probabilities = toplot[:,2]




    # if there is a divide by zero error, do a dummy calculation and replace infinity with the max non-infinite number
    with np.errstate(divide='ignore'):
        if detections==0:
            trial_transrate = -np.log(1-confidence)/(probabilities)/(tsurvey + durations)/area

            try: 
                ultransrates = np.nan_to_num(-np.log(1-confidence)/(probabilities)/(tsurvey + durations)/area, posinf=np.max(trial_transrate[trial_transrate < np.inf]))
                ulZrate = interpolate.griddata(toplot[:,0:2], ultransrates, (X, Y), method='linear')
                # Make plot for zero detections
                fig = plt.figure()
                # 
                # https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py
                lev_exp = np.linspace(np.floor(np.log10(ulZrate[ulZrate > 0].min())),np.ceil(np.log10(np.mean(ulZrate[ulZrate > 0]))+1), num=1000)
                levs = np.power(10, lev_exp)
                # cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())
                # levels = np.geomspace(max(np.amin(toplot[:,2]),1e-16),np.mean(toplot[:,2]),num = 1000)
                csrate = plt.contourf(10**X, 10**Y, ulZrate, levels=levs, cmap='viridis', norm=colors.LogNorm())
                cbarrate = fig.colorbar(csrate, ticks=np.geomspace(ulZrate.min(),np.ceil(np.mean(ulZrate))+1,num=10), format=ticker.StrMethodFormatter("{x:01.1e}"))
                cbarrate.set_label('Transient Rate per day per sq. deg.')
                plt.plot(10**xs, 10**np.full(xs.shape, np.log10(vlinex[0])),  color="red")
                ax = plt.gca()
                ax.set_ylabel('Characteristic Flux (Jy)')
                ax.set_xlabel('Characteristic Duration (Days)')
                ax.set_xscale('log')
                ax.set_yscale('log')
                ax.set_xlim(10**np.amin(X),10**np.amax(X))
                ax.set_ylim(10**np.amin(Y),10**np.amax(Y))
                plt.savefig('rateplot'+rgn+'.png')
                zkey = np.argsort(ultransrates[toplot[:,0] > vlinex[0]])
                xratesorted = 10**toplot[:,0][toplot[:,0] > vlinex[0]][zkey][0:10]
                yratesorted = 10**toplot[:,1][toplot[:,0] > vlinex[0]][zkey][0:10]
                zratesorted = ultransrates[toplot[:,0] > vlinex[0]][zkey][0:10]
                print("Ten lowest transient rates above the 99% false detection line:")
                print("Tau, Fpk, Rate")
                for x,y,z in zip(xratesorted, yratesorted, zratesorted):
                    print(x,y,z)

                plt.close()
            except ValueError:
                print("Issues calculating rates, skipping")
                pass

        else:
            from scipy.special import gammaincinv
            alpha = 1-confidence
            upperlimitpoisson = gammaincinv(detections+1, 1-alpha/2)
            lowerlimitpoisson = gammaincinv(detections,alpha/2.)
            lltrial_transrate, ultrial_transrate = lowerlimitpoisson/(probabilities)/(tsurvey + durations)/area, (upperlimitpoisson/(probabilities)/(tsurvey + durations)/area)
            lltransrates = np.nan_to_num(lowerlimitpoisson/(probabilities)/(tsurvey + durations)/area, posinf=np.max(lltrial_transrate[lltrial_transrate < np.inf]))
            ultransrates = np.nan_to_num(upperlimitpoisson/(probabilities)/(tsurvey + durations)/area, posinf=np.max(ultrial_transrate[ultrial_transrate < np.inf]))

            ulZrate = interpolate.griddata(toplot[:,0:2], ultransrates, (X, Y), method='linear')
            llZrate = interpolate.griddata(toplot[:,0:2], lltransrates, (X, Y), method='linear')
            # Make upper and lower limit plots 
            # https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py

            fig = plt.figure()
            #
            lev_exp = np.linspace(np.floor(np.log10(ulZrate[ulZrate>0].min())),np.ceil(np.log10(np.mean(ulZrate[ulZrate>0])+1)), num=1000)
            levs = np.power(10, lev_exp)
            # cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())
            # levels = np.geomspace(max(np.amin(toplot[:,2]),1e-16),np.mean(toplot[:,2]),num = 1000)
            csrate = plt.contourf(10**X, 10**Y, ulZrate, levels=levs, cmap='viridis', norm=colors.LogNorm())
            cbarrate = fig.colorbar(csrate, ticks=np.geomspace(ulZrate.min(),np.ceil(np.mean(ulZrate))+1,num=10), format=ticker.StrMethodFormatter("{x:01.1e}"))
            cbarrate.set_label('Transient Rate per day per sq. deg.')
            plt.plot(10**xs, 10**np.full(xs.shape, np.log10(vlinex[0])),  color="red")
            ax = plt.gca()
            ax.set_ylabel('Characteristic Flux (Jy)')
            ax.set_xlabel('Characteristic Duration (Days)')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(10**np.amin(X),10**np.amax(X))
            ax.set_ylim(10**np.amin(Y),10**np.amax(Y))
            plt.savefig('ulrateplot'+rgn+'.png')
            zkey = np.argsort(ultransrates[toplot[:,0] > vlinex[0]])
            xratesorted = 10**toplot[:,0][toplot[:,0] > vlinex[0]][zkey][0:10]
            yratesorted = 10**toplot[:,1][toplot[:,0] > vlinex[0]][zkey][0:10]
            zratesortedul = ultransrates[toplot[:,0] > vlinex[0]][zkey][0:10]
            zratesortedll = lltransrates[toplot[:,0] > vlinex[0]][zkey][0:10]

            plt.close()
            # https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py

            fig = plt.figure()
            # 
            lev_exp = np.linspace(np.floor(np.log10(llZrate[llZrate > 0].min())),np.ceil(np.log10(np.mean(llZrate[llZrate > 0]))+1), num=1000)
            levs = np.power(10, lev_exp)
            # cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())
            # levels = np.geomspace(max(np.amin(toplot[:,2]),1e-16),np.mean(toplot[:,2]),num = 1000)
            csrate = plt.contourf(10**X, 10**Y, llZrate, levels=levs, cmap='viridis', norm=colors.LogNorm())
            cbarrate = fig.colorbar(csrate, ticks=np.geomspace(llZrate[llZrate>0].min(),np.ceil(np.mean(llZrate[llZrate>0]))+1,num=10), format=ticker.StrMethodFormatter("{x:01.1e}"))
            cbarrate.set_label('Transient Rate per day per sq. deg.')
            plt.plot(10**xs, 10**np.full(xs.shape, np.log10(vlinex[0])),  color="red")
            ax = plt.gca()
            ax.set_ylabel('Characteristic Flux (Jy)')
            ax.set_xlabel('Characteristic Duration (Days)')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(10**np.amin(X),10**np.amax(X))
            ax.set_ylim(10**np.amin(Y),10**np.amax(Y))
            plt.savefig('llrateplot'+rgn+'.png')
            print("Ten lowest transient rates above the 99% false detection line:")
            print("Tau, Fpk, Rate")
            for x,y,zll,zul in zip(xratesorted, yratesorted, zratesortedll, zratesortedul):
                print(x,y,zll,'~',zul)

            plt.close()


    
    


    # print("Error below here")
    fig = plt.figure()
    cs = plt.contourf(10**X, 10**Y, Z, levels=np.linspace(0,1.0,num = int(1/0.01)+1), cmap='viridis')
    cbar = fig.colorbar(cs, ticks=np.linspace(0,1.0,num=11))
    durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices = lclines(xs, ys, durmax, max_distance, flux_err, obs)   
    plt.plot(10**xs[durmax_y_indices], durmax_y[durmax_y_indices],  color = "red")
    plt.plot(10**xs[maxdist_y_indices], maxdist_y[maxdist_y_indices],   color = "red")
    plt.plot(10**xs, 10**np.full(xs.shape, np.log10(vlinex[0])),  color="red")
    if durmax_x[0]!=' ':
        plt.plot(10**durmax_x, 10**ys,  color = "red")
    if maxdist_x[0]!=' ':    
        plt.plot(10**maxdist_x, 10**ys,  color = "red")
    if (np.amin(day1_obs_x) > np.amin(10**ys)): plt.plot(day1_obs_x, 10**ys,  color = "red")
    ax = plt.gca()
    ax.set_ylabel('Characteristic Flux (Jy)')
    ax.set_xlabel('Characteristic Duration (Days)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(10**np.amin(X),10**np.amax(X))
    ax.set_ylim(10**np.amin(Y),10**np.amax(Y))
    plt.savefig('probcont'+rgn+'.png')
    print("Saved probability contour plot to probcont"+rgn+".png")
    plt.close()
    
    plt.bar(fddetbins[0:-1][fddethist>0],fddethist[fddethist>0], width = (fddetbins[1:]-fddetbins[:-1])[fddethist>0], align='edge', alpha=0.5)
    plt.plot(vlinex, vliney, color="red")
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:1.2e}"))
    ax.set_xlabel('Actual Transient Flux (Jy)')
    ax.set_ylabel('Number of Transients')
    plt.savefig('FalseDetections'+rgn+'.png')
    print("Saved False Detection histogram to FalseDetections"+rgn+".png")
    plt.close()

    plt.bar(sensbins[0:-1][senshist>0],senshist[senshist>0], width = (sensbins[1:]-sensbins[:-1])[senshist>0], align='edge', alpha=0.5, color='gray')
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:1.2e}"))
    ax.set_xlabel('Observation Noise (Jy)')
    ax.set_ylabel('Number of Observations')
    plt.savefig('Sensitivities'+rgn+'.png')
    print("Saved observation sensitivity histogram to Sensitivites"+rgn+".png")
    plt.close()
    print("Dumping numpy arrays for you to use.")
    np.save("fddetbins"+rgn+".npy",fddetbins)
    np.save("fddethist"+rgn+".npy",fddethist)
    np.save("senshist"+rgn+".npy",senshist)
    np.save("sensbins"+rgn+".npy",sensbins)
    np.save("vlinex"+rgn+".npy",vlinex)
    np.save("vliney"+rgn+".npy",vliney)
