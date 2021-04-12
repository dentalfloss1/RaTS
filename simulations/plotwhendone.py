import configparser
import random
import time
import datetime
import numpy as np
import os
import argparse
import math
import warnings
from bokeh.plotting import figure, show, output_file
from bokeh.models import LinearColorMapper, SingleIntervalTicker, ColorBar, Title
from bokeh.io import export_png
import scipy.interpolate as interpolate
from tqdm import tqdm
from astropy import units as u 
from astropy.coordinates import SkyCoord
from scipy.special import binom
import ered
warnings.simplefilter("error", RuntimeWarning)

ered_cls = ered()
obs_setup = 'obsnoise21012012h32m50s.txt'
tstart, tdur, sens, ra, dec  = np.loadtxt(obs_setup, unpack=True, delimiter = ',',
            dtype={'names': ('start', 'duration','sens', 'ra', 'dec'), 'formats': ('U32','f8','f8','f8','f8')})
        
tstart = np.array([time.mktime(datetime.datetime.strptime(t, "%Y-%m-%dT%H:%M:%S.%f+00:00").timetuple())/(3600*24) for t in tstart])
sortkey = np.argsort(tstart)
obs = np.zeros((len(tstart),3))
obs[:,0]+=tstart[sortkey]
obs[:,1]+=tdur[sortkey]/3600/24 # convert into days
obs[:,2]+=sens[sortkey]
pointing = np.array([np.array([r,d]) for r,d in zip(ra[sortkey],dec[sortkey])])


file='plotout'
extra_threshold = 3
det_threshold = 5
flux_err = 1e-1
toplot = np.loadtxt('obsnoise21012012h32m50s.txt', skiprows=1)
gaussiancutoff =  0.1
lclines = ered.lines


lightcurve = ''
toplot[:,0] = np.log10(toplot[:,0])
toplot[:,1] = np.log10(toplot[:,1])

gaps = np.array([],dtype=np.float32)
for i in range(len(obs)-1):
    gaps = np.append(gaps, obs[i+1,0] - obs[i,0] + obs[i,1])
    # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
min_sens = min(obs[:,2])
max_sens = max(obs[:,2])
extra_thresh = max_sens / det_threshold * (extra_threshold + det_threshold)
sens_last = obs[-1,2]
sens_maxgap = obs[np.where((gaps[:] == max(gaps)))[0]+1 ,2][0]

durmax = obs[-1,0] + obs[-1,1] - obs[0,0]
day1_obs = obs[0,1]
max_distance = max(gaps)

dmin=min(toplot[:,0])
dmax=max(toplot[:,0])
flmin=min(toplot[:,1])
flmax=max(toplot[:,1])

xs = np.arange(dmin, dmax, 1e-3)
ys = np.arange(flmin, flmax, 1e-3)

day1_obs_x = np.empty(len(ys))
day1_obs_x.fill(day1_obs)

sensmin_y = np.empty(len(xs))
sensmin_y.fill(min_sens)

sensmax_y = np.empty(len(xs))
sensmax_y.fill(max_sens)

extra_y = np.empty(len(xs))
extra_y.fill(extra_thresh)

X = np.linspace(dmin, dmax, num = 1000)
Y = np.linspace(flmin, flmax, num = 1000)

X, Y = np.meshgrid(X, Y)

Z = interpolate.griddata(toplot[:,0:2], toplot[:,2], (X, Y), method='linear')
p = figure(title="Probability Contour Plot",tooltips = [("X", "$X"), ("Y", "$Y"), ("value", "@image")], x_axis_type = "log", y_axis_type = "log")
p.x_range.range_padding = p.y_range.range_padding = 0
color_mapper = LinearColorMapper(palette="Viridis256",low = 0.0, high = 1.0)
color_bar = ColorBar(color_mapper=color_mapper, ticker=SingleIntervalTicker(interval = 0.1), label_standoff=12, border_line_color=None, location=(0,0))
p.image(image=[Z], x=np.amin(10**xs), y=np.amin(10**ys), dw=(np.amax(10**xs)-np.amin(10**xs)), dh=(np.amax(10**ys)-np.amin(10**ys)),palette="Viridis256")
durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices = lclines(xs, ys, durmax, max_distance, flux_err, obs)   
p.line(10**xs[durmax_y_indices], durmax_y[durmax_y_indices],  line_width=2, line_color = "red")
p.line(10**xs[maxdist_y_indices], maxdist_y[maxdist_y_indices],  line_width=2, line_color = "red")
if durmax_x[0]!=' ':
    p.line(10**durmax_x, 10**ys,   line_width=2, line_color = "red")
if maxdist_x[0]!=' ':    
    p.line(10**maxdist_x, 10**ys,  line_width=2, line_color = "red")
if (np.amin(day1_obs_x) > np.amin(10**ys)): p.line(day1_obs_x, 10**ys,  line_width=2, line_color = "red")
if (sensmin_y[0] > np.amin(10**ys)): p.line(10**xs, sensmin_y,  line_width=2, line_color = "red")
if (sensmax_y[0] > np.amin(10**ys)): p.line(10**xs, sensmax_y,  line_width=2, line_color = "red")
if (extra_y[0] > np.amin(10**ys)): p.line(10**xs, extra_y,  line_width=2, line_color = "red")
p.add_layout(color_bar, 'right')
p.add_layout(Title(text="Duration (days)", align="center"), "below")
p.add_layout(Title(text="Transient Flux Density (Jy)", align="center"), "left")
p.toolbar.logo = None
p.toolbar_location = None
p.toolbar.active_drag = None
p.toolbar.active_scroll = None
p.toolbar.active_tap = None

output_file(file + "_ProbContour.html", title = "Probability Contour")
export_png(p, filename=file + "_ProbContour.png")
show(p)