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
import matplotlib.pyplot as plt
from matplotlib import ticker, colors

def observing_strategy(obs_setup, det_threshold, nobs, obssens, obssig, obsinterval, obsdurations):
    """Parse observation file or set up trial mode. Return array of observation info and a regions observed"""

    rng = np.random.default_rng()
    start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)
    if obs_setup is not None:
        tstart, tend, sens, ra, dec, gapsfile,fov  = np.loadtxt(obs_setup, unpack=True, delimiter = ',',
            dtype={'names': ('start', 'end','sens', 'ra', 'dec','gaps', 'fov'), 'formats': ('U32','U32','f8','f8','f8','<U128','f8')})
        tstart = np.array([(datetime.datetime.strptime(t, "%Y-%m-%dT%H:%M:%S.%f+00:00") - start_epoch).total_seconds()/3600/24 for t in tstart])
        tend = np.array([(datetime.datetime.strptime(t, "%Y-%m-%dT%H:%M:%S.%f+00:00") - start_epoch).total_seconds()/3600/24 for t in tend])
        tdur = tend - tstart
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
        # pointing[::3]-=[tmpsc.ra.deg - tmpscoff1.ra.deg,tmpsc.dec.deg - tmpscoff1.dec.deg]
        # pointing[1::3]-=[tmpsc.ra.deg - tmpscoff2.ra.deg,tmpsc.dec.deg - tmpscoff2.dec.deg]
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
                    halfheightr4 = np.arccos(np.cos(r1)/np.cos(gamma[i][j])) 
                    point4key = np.where(regions['identity'] == str(i)+'&'+str(j))
                    point4ra = regions['ra'][point4key]
                    point4dec = regions['dec'][point4key]
                    point4sc = SkyCoord(ra=point4ra, dec=point4dec, unit='deg',frame='fk5')
                    point4pa = point4sc.position_angle(uniquesky[j])
                    point7sc = point4sc.directional_offset_by(point4pa + angle_offset, halfheightr4)
                    # print(point7sc.separation(uniquesky[index3]).deg)
                    if point7sc.separation(uniquesky[index3]).deg > uniquepoint[index3,2]:
                        point7sc = point4sc.directional_offset_by(point4pa - angle_offset, halfheightr4)

                    halfheightr5 = np.arccos(np.cos(r2)/np.cos(gamma[j][index3]))
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

    #                 scatterpointsra.extend([point7sc.ra,point8sc.ra,point9sc.ra])
    #                 scatterpointsdec.extend([point7sc.dec,point8sc.dec,point9sc.dec])
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
    # print(regions)
    # for i in range(len(uniquepoint)):
    #     for j in range(len(uniquepoint),len(regions[regions['identity'] != ''])):
    #         print(i,j)
    #         print(regions[regions['identity'] != '']['identity'][j])
    #         print(regions[regions['identity'] != '']['identity'][i])
    #         print(regions[regions['identity'] != '']['area'][j]/regions[regions['identity'] != '']['area'][i])

    # exit()
    return regions[regions['identity'] != ''], obssubsection
    

def generate_pointings(n_sources, pointFOV, i, leftoff, overlapnums):
    """Simulate pointings for each simulated source. Use a monte-carlo like method to roll the dice to determine position. Return a numpy array and """
    
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
        sc.data.lon[reject] = ra_randfield(len(sc[reject]),uniquepointFOV[i])
        sc.data.lat[reject] = dec_randfield(len(sc[reject]),uniquepointFOV[i])
        sc.cache.clear()

    numrgns = len(uniqueskycoord)
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
    # print(overlapnums)
    # print(leftoff)
    return overlapnums,leftoff

def generate_sources(n_sources, start_survey, end_survey, fl_min, fl_max, dmin, dmax, lightcurve, burstlength, burstflux):
    """Generate characteristic fluxes and characteristic durations for simulated sources. Return as numpy array"""
    
    start_epoch = datetime.datetime(1858, 11, 17, 00, 00, 00, 00)
    rng = np.random.default_rng()
    bursts = np.zeros(n_sources, dtype={'names': ('chartime', 'chardur','charflux'), 'formats': ('f8','f8','f8')}) # initialise the numpy array
    if not np.isnan(burstlength):
        bursts['chardur'] += burstlength
    else:
        bursts['chardur'] = (rng.random(n_sources)*(dmax - dmin) + dmin) # random number for duration
        # bursts['chardur'] = (rng.random(n_sources)*(dmax - dmin) + dmin) # random number for duration
    # bursts['chardur'] += 500*7 + 0.01
    if not np.isnan(burstflux):
        bursts['charflux'] += burstflux
    else:
        bursts['charflux'] = (rng.random(n_sources)*(fl_max - fl_min) + fl_min) # random number for flux
    # bursts['charflux'] = (rng.random(n_sources)*(fl_max - fl_min) + fl_min) # random number for flux
    if hasattr(lightcurve, 'docustompop'):
        exit()
    return bursts
    
def generate_start(bursts, potential_start, potential_end, n_sources):
    """Generate characteristic times to go with durations and fluxes. Return modified numpy array"""
    
    rng = np.random.default_rng() 
    bursts['chartime'] = rng.random(n_sources)*(potential_end - potential_start) + potential_start
    bursts.sort(axis=0,order='chartime')
    return bursts
    

def detect_bursts(obs, flux_err,  det_threshold, sources, fluxint):
    """Detect simulated sources by using a series of conditionals along with the integrated flux calculation. Returns detected sources and the boolean array to get source indices"""
    
    rng = np.random.default_rng() 
    if (obs['gaps']!="False").any():
        detections = np.zeros((len(obs),len(sources)),dtype=bool)
        detbool = np.zeros(len(sources)*len(obs),dtype=bool)
        for i,o in enumerate(obs):
            if o['gaps']!=False:
                subobs, _ = observing_strategy(o['gaps'], det_threshold, 1, 1, 1, 1, 1) # We are giving the scansfile name, so the other variables are unimportant, we set them to 1 
                flux_int = np.zeros((len(sources)),dtype=np.float32)
                subfluxint = np.zeros([len(subobs),len(sources)], dtype=np.float32)
                for j,s in enumerate(subobs):
                    end_subobs = s['start']+s['duration']
                    start_subobs = s['start']
                    F0_o = sources['charflux']
                    error = np.sqrt((F0_o * flux_err)**2 + (o['sens']/det_threshold)**2) 
                    F0 =rng.normal(F0_o, error)
                    F0[(F0<0)] = F0[(F0<0)]*0
                    tcrit = sources['chartime']
                    tau = sources['chardur']
                    subfluxint[j] = fluxint(F0, tcrit, tau, end_subobs, start_subobs)
                    subfluxint[j][subfluxint[j] < 0] = 0
                flux_int = np.average(subfluxint,weights=subobs['duration']/np.sum(subobs['duration']),axis=0)
                flux_int[flux_int < 0] = 0
            else:
                flux_int = np.zeros((len(sources)),dtype=np.float32)
                F0_o = sources['charflux']
                tcrit = sources['chartime']
                tau = sources['chardur']
                end_obs = o['start']+o['duration']
                start_obs = o['start']
                error = np.sqrt((F0_o * flux_err)**2 + (o['sens']/det_threshold)**2) 
                F0 =rng.normal(F0_o, error)
                F0[(F0<0)] = F0[(F0<0)]*0
                flux_int = fluxint(F0, tcrit, tau, end_obs, start_obs)
                flux_int[flux_int < 0] = 0
            sensitivity = o['sens']
            detections[i] = (flux_int > sensitivity)
        constant = np.all(detections==True, axis=0)
        detectany = np.any(detections==True,axis=0)
        detections = detectany & np.logical_not(constant)
            
    else:  # vectorized
        sensitivity = np.tile(obs['sens'],len(sources))
        F0_o = np.repeat(sources['charflux'],len(obs))
        error = np.sqrt((F0_o * flux_err)**2 + (sensitivity/det_threshold)**2) 
        # error = F0_o*flux_err
        F0 =rng.normal(F0_o, error)
        # F0=F0_o
        F0[(F0<0)] = F0[(F0<0)]*0
        t0 = np.repeat(sources['chartime'],len(obs))
        tau0 = np.repeat(sources['chardur'],len(obs))
       #  print(tau0.min(),tau0.max(),obs['duration'].min(),obs['duration'].max())
       #  print("^durations")
        end_obs = np.tile(obs['start']+obs['duration'],len(sources))
        start_obs = np.tile(obs['start'],len(sources))
        flux_int = fluxint(F0, t0, tau0, end_obs, start_obs) # uses whatever class of lightcurve supplied: tophat, ered, etc      
       #  for f,t,d,t1,t2,fi in zip(F0,t0,tau0,start_obs,end_obs,flux_int):
       #      print(f,t,d,t1,t2,fi)
        flux_int[flux_int < 0] = 0
        detections = (flux_int > sensitivity).reshape(len(sources),len(obs)).transpose()
        constant = np.all(detections==True, axis=0)
        detectany = np.any(detections==True,axis=0)
        nondetect = np.all(detections==False,axis=0)
       #  print(np.sum(constant), "constant sources")
       #  print(np.sum(detectany), "total detections")
       #  print(np.sum(nondetect), "total undetected")
        detections = detectany & np.logical_not(constant)

    return sources[detections], detections

def statistics(fl_min, fl_max, dmin, dmax, det, all_simulated):
    """Calculate probabilities based on detections vs simulated, return a numpy array"""

    start = datetime.datetime.now()
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

    end = datetime.datetime.now()
    return stats

def make_mpl_plots(rgn, fl_min,fl_max,dmin,dmax,det_threshold,extra_threshold,obs,cdet,file,flux_err,toplot,gaussiancutoff,lclines,area,tsurvey,detections,confidence,filename):
    """Use Matplotlib to make plots and if that fails dump numpy arrays. Returns an int that indicates plotting success or failure"""
    fddethist = None
    fddetbins = None

    start = datetime.datetime.now()
    # Make histograms of observation senstivities and false detections
    fluxbins = np.geomspace(fl_min, fl_max, num=int(round((np.log10(fl_max)-np.log10(fl_min))/0.01)), endpoint=True)
    if cdet:
        fddethist, fddetbins = cdet 
    senshist, sensbins = np.histogram(obs['sens'], bins=fluxbins, density=False)
    # Calculate 99% of false detections line (from the left)
    histsum = 0
    if cdet:
        for j in range(len(fddetbins)):
            histsum += fddethist[j]
            try:
                if (histsum)/np.sum(fddethist) >= 0.99:
                    print("stopped on index:",j,"exact percentage threshold:", (histsum)/np.sum(fddethist))
                    vlinex = np.full(10, fddetbins[j])
                    vliney = np.linspace(1, np.amax([fddethist.max(),senshist.max()]), num=10)
                    break
            except RuntimeWarning:
                vlinex = np.full(10,fl_min)
                vliney = np.linspace(1,np.amax([fddethist.max(),senshist.max()]), num=10)
    else:
        vlinex = np.full(10,fl_min)
        vliney = np.linspace(1,np.amax(senshist), num=10)
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
    mindurationobs = np.amin(obs['duration'])
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
            ultransrates = -np.log(1-confidence)/(probabilities)/(tsurvey + 10**durations)/area
            # try: 
            # ultransrates = np.nan_to_num(-np.log(1-confidence)/(probabilities)/(tsurvey + durations)/area, posinf=np.max(trial_transrate[trial_transrate < np.inf]))
            ulZrate = interpolate.griddata(toplot[:,0:2], ultransrates, (X, Y), method='linear')
            # figsc = plt.figure()
            # plt.scatter(10**X[-2,:],ulZrate[-2,:], s=3)
            # ax = plt.gca()
            # ax.set_xscale('log')
            # ax.set_yscale('log')
            # ax.set_xlabel('Characteristic Flux (Jy)')
            # ax.set_ylabel('Transient Rate per day per sq. deg.')
            # ax.set_title('Transient Rate at '+str(10**Y[-2,0])+' Jy')
            # plt.savefig(filename+'ratescatter.png')
            # Make plot for zero detections
            fig = plt.figure()
            # 
            # https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py
            
            lev_exp = np.linspace(np.log10(ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].min()),np.log10(ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].max()), num=1000)
            levs = np.power(10, lev_exp)
            # cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())
            # levels = np.geomspace(max(np.amin(toplot[:,2]),1e-16),np.mean(toplot[:,2]),num = 1000)
            ulcsrate = plt.contourf(10**X, 10**Y, ulZrate, levels=levs, cmap='viridis', norm=colors.LogNorm(), vmin=ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].min(), vmax=ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].max())
            ulrateticks = np.geomspace(ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].min(),ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].max(),num=10)
            cbarrate = fig.colorbar(ulcsrate, ticks=ulrateticks, format=ticker.StrMethodFormatter("{x:01.1e}"))
            cbarrate.set_label('Transient Rate per day per sq. deg.')
            plt.plot(10**xs, 10**np.full(xs.shape, np.log10(vlinex[0])),  color="red")
            plt.plot(10**np.full(ys.shape, np.log10(5/60/24)),10**ys,  color="red")
            ax = plt.gca()
            ax.set_ylabel('Characteristic Flux (Jy)')
            ax.set_xlabel('Characteristic Duration (Days)')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(10**np.amin(X),10**np.amax(X))
            ax.set_ylim(10**np.amin(Y),10**np.amax(Y))
            plt.savefig(filename+'rateplot'+rgn+'.png')
            zkey = np.argsort(ultransrates[toplot[:,0] > vlinex[0]])
            xratesorted = 10**toplot[:,0][toplot[:,0] > vlinex[0]][zkey][0:10]
            yratesorted = 10**toplot[:,1][toplot[:,0] > vlinex[0]][zkey][0:10]
            zratesorted = ultransrates[toplot[:,0] > vlinex[0]][zkey][0:10]
            print("Ten lowest transient rates above the 99% false detection line:")
            print("Tau, Fpk, Rate")
            for x,y,z in zip(xratesorted, yratesorted, zratesorted):
                print(x,y,z)

            plt.close()
            # except ValueError:
            #     print("Issues calculating rates, skipping")
            #     pass

        else:
            from scipy.special import gammaincinv
            alpha = 1-confidence
            upperlimitpoisson = gammaincinv(detections+1, 1-alpha/2)
            lowerlimitpoisson = gammaincinv(detections,alpha/2.)
            lltrial_transrate, ultrial_transrate = lowerlimitpoisson/(probabilities)/(tsurvey + 10**durations)/area, (upperlimitpoisson/(probabilities)/(tsurvey + 10**durations)/area)
            lltransrates = lowerlimitpoisson/(probabilities)/(tsurvey + 10**durations)/area
            ultransrates = upperlimitpoisson/(probabilities)/(tsurvey + 10**durations)/area

            ulZrate = interpolate.griddata(toplot[:,0:2], ultransrates, (X, Y), method='linear')
            llZrate = interpolate.griddata(toplot[:,0:2], lltransrates, (X, Y), method='linear')
            # Make upper and lower limit plots 
            # https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py

            fig = plt.figure()
            #
            lev_exp = np.linspace(np.log10(ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].min()),np.log10(ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].max()), num=1000)
            levs = np.power(10, lev_exp)
            # cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())
            # levels = np.geomspace(max(np.amin(toplot[:,2]),1e-16),np.mean(toplot[:,2]),num = 1000)
            ulcsrate = plt.contourf(10**X, 10**Y, ulZrate, levels=levs, cmap='viridis', norm=colors.LogNorm(), vmin=ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].min(), vmax=ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].max())
            ulrateticks = np.geomspace(ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].min(),ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].max(),num=10)
            cbarrate = fig.colorbar(ulcsrate, ticks=ulrateticks, format=ticker.StrMethodFormatter("{x:01.1e}"))
            cbarrate.set_label('Transient Rate per day per sq. deg.')
            plt.plot(10**xs, 10**np.full(xs.shape, np.log10(vlinex[0])),  color="red")
            ax = plt.gca()
            ax.set_ylabel('Characteristic Flux (Jy)')
            ax.set_xlabel('Characteristic Duration (Days)')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(10**np.amin(X),10**np.amax(X))
            ax.set_ylim(10**np.amin(Y),10**np.amax(Y))
            plt.savefig(filename+'ulrateplot'+rgn+'.png')
            zkey = np.argsort(ultransrates[toplot[:,0] > vlinex[0]])
            xratesorted = 10**toplot[:,0][toplot[:,0] > vlinex[0]][zkey][0:10]
            yratesorted = 10**toplot[:,1][toplot[:,0] > vlinex[0]][zkey][0:10]
            zratesortedul = ultransrates[toplot[:,0] > vlinex[0]][zkey][0:10]
            zratesortedll = lltransrates[toplot[:,0] > vlinex[0]][zkey][0:10]

            plt.close()
            # https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py

            fig = plt.figure()
            # 
            lev_exp = np.linspace(np.log10(llZrate[(llZrate > 0) & (llZrate != np.inf)].min()),np.log10(llZrate[(llZrate > 0) & (llZrate != np.inf)].max()), num=1000)
            levs = np.power(10, lev_exp)
            # cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())
            # levels = np.geomspace(max(np.amin(toplot[:,2]),1e-16),np.mean(toplot[:,2]),num = 1000)
            llcsrate = plt.contourf(10**X, 10**Y, llZrate, levels=levs, cmap='viridis', norm=colors.LogNorm(), vmin=llZrate[(llZrate > 0) & (llZrate != np.inf)].min(), vmax=llZrate[(llZrate > 0) & (llZrate != np.inf)].max())
            llrateticks = np.geomspace(llZrate[(llZrate > 0) & (llZrate != np.inf)].min(),llZrate[(llZrate > 0) & (llZrate != np.inf)].max(),num=10)
            cbarrate = fig.colorbar(llcsrate, ticks=llrateticks, format=ticker.StrMethodFormatter("{x:01.1e}"))
            cbarrate.set_label('Transient Rate per day per sq. deg.')
            plt.plot(10**xs, 10**np.full(xs.shape, np.log10(vlinex[0])),  color="red")
            plt.plot(10**np.full(ys.shape, np.log10(5/60/24)),10**ys,  color="red")
            ax = plt.gca()
            ax.set_ylabel('Characteristic Flux (Jy)')
            ax.set_xlabel('Characteristic Duration (Days)')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(10**np.amin(X),10**np.amax(X))
            ax.set_ylim(10**np.amin(Y),10**np.amax(Y))
            plt.savefig(filename+'llrateplot'+rgn+'.png')
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
    if (mindurationobs > np.amin(10**ys)): plt.plot(np.full(ys.shape, mindurationobs), 10**ys,  color = "red")
    ax = plt.gca()
    ax.set_ylabel('Characteristic Flux (Jy)')
    ax.set_xlabel('Characteristic Duration (Days)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(10**np.amin(X),10**np.amax(X))
    ax.set_ylim(10**np.amin(Y),10**np.amax(Y))
    plt.savefig(filename+'probcont'+rgn+'.png')
    print("Saved probability contour plot to probcont"+rgn+".png")
    plt.close()
    if cdet: 
        plt.bar(fddetbins[0:-1][fddethist>0],fddethist[fddethist>0], width = (fddetbins[1:]-fddetbins[:-1])[fddethist>0], align='edge', alpha=0.5, color='gray')
        plt.plot(vlinex, vliney, color="red")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:1.2e}"))
        ax.set_xlabel('Actual Transient Flux (Jy)')
        ax.set_ylabel('Number of Transients')
        plt.savefig(filename+'FalseDetections'+rgn+'.png')
        print("Saved False Detection histogram to FalseDetections"+rgn+".png")
        plt.close()

    plt.bar(sensbins[0:-1][senshist>0],senshist[senshist>0], width = (sensbins[1:]-sensbins[:-1])[senshist>0], align='edge', alpha=0.5, color='gray')
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:1.2e}"))
    ax.set_xlabel('Observation Noise (Jy)')
    ax.set_ylabel('Number of Observations')
    plt.savefig(filename+'Sensitivities'+rgn+'.png')
    print("Saved observation sensitivity histogram to Sensitivites"+rgn+".png")
    plt.close()
    print("Dumping numpy arrays for you to use.")
    now = (datetime.datetime.now() - datetime.datetime(1858, 11, 17, 00, 00, 00, 00)).total_seconds()/60/60/24
    if detections==0:
        np.savez_compressed(filename+"myrun"+str(now).replace('.','_')+".npz", 
            fddetbins=fddetbins, 
            fddethist=fddethist, 
            senshist=senshist, 
            sensbins=sensbins, 
            X=X,
            Y=Y,
            Z=Z,
            ulZrate=ulZrate,
            ulcsrate=ulcsrate,
            ulrateticks=ulrateticks,
            vlinex=vlinex,
            xs=xs,
            ys=ys,
            durmax_x=durmax_x,
            durmax_y=durmax_y,
            durmax_y_indices=durmax_y_indices,
            maxdist_x=maxdist_x,
            maxdist_y=maxdist_y,
            maxdist_y_indices=maxdist_y_indices,
            day1_obs_x=day1_obs_x)
    else:
        np.savez_compressed(filename+"myrun"+str(now).replace('.','_')+".npz", 
            fddetbins=fddetbins, 
            fddethist=fddethist, 
            senshist=senshist, 
            sensbins=sensbins, 
            X=X,
            Y=Y,
            Z=Z,
            ulZrate=ulZrate,
            llZrate=llZrate,
            ulcsrate=ulcsrate,
            llcsrate=llcsrate,
            ulrateticks=ulrateticks,
            llrateticks=llrateticks,
            vlinex=vlinex,
            xs=xs,
            ys=ys,
            durmax_x=durmax_x,
            durmax_y=durmax_y,
            durmax_y_indices=durmax_y_indices,
            maxdist_x=maxdist_x,
            maxdist_y=maxdist_y,
            maxdist_y_indices=maxdist_y_indices,
            day1_obs_x=day1_obs_x)
    end = datetime.datetime.now()
    print('Plotting elapsed: ',end-start)

