from bokeh.plotting import figure, show, output_file
from bokeh.models import LinearColorMapper, SingleIntervalTicker, ColorBar, Title, Range1d, ColumnDataSource
import random
from bokeh.io import export_png
import numpy as np 
import warnings
import configparser
# import matplotlib.pyplot as plt
from argparse import ArgumentParser
from datetime import datetime,timedelta
from scipy.special import binom
from scipy.special import gammaincinv
warnings.simplefilter("error", RuntimeWarning)

# Parse command line input
parser = ArgumentParser()
parser.add_argument('obsfile', help='Supply a file formatted with: "noise,filename"')
args = parser.parse_args()

noise, filenames = np.loadtxt(args.obsfile, delimiter = ',', unpack=True, dtype=str)

# Read surfdens.ini settings into variables

params = configparser.ConfigParser()
params.read('surfdens.ini')
conf_lev = np.float(params['STATISTICAL']['confidence'])
extract_rad = np.float(params['DATA']['extract_rad'])
sigtonoise = np.float(params['STATISTICAL']['sigtonoise'])
tsnap = np.float(params['DATA']['minint'])/60./60./24.
num_skyrgns = int(params['DATA']['num_skyrgns'])
detections = int(params['DATA']['detections'])

if detections > 0:
    alpha = 1-conf_lev
    upperlimitpoisson = gammaincinv(detections+1, 1-alpha/2)
    lowerlimitpoisson = gammaincinv(detections,alpha/2.)



def trad_nondet_surf(num_ims, num_skyrgns):
    if detections==0:
        return  np.log(1.0-conf_lev)/(-(num_ims-num_skyrgns)*np.pi*extract_rad**2)
    else:
        return lowerlimitpoisson/((num_ims-num_skyrgns)*np.pi*extract_rad**2), upperlimitpoisson/((num_ims-num_skyrgns)*np.pi*extract_rad**2)


noise = np.array(noise)
noise = noise.astype(np.float)

noise_levs = np.array([max(noise)],dtype = np.float)

for i in range(1,int(len(noise)-num_skyrgns)):
    noise_levs = np.append(noise_levs,max(noise[noise < noise_levs[i-1]])) 
noise_levs_rates = np.array([trad_nondet_surf(j, num_skyrgns) for j in range(len(noise),num_skyrgns, -1)])


print(noise_levs_rates)
exit()


plot = figure(title=" ", x_axis_type = "log", y_axis_type = "log" )
plot.cross(x=sigtonoise*max(noise), y=trad_nondet_surf(len(noise), num_skyrgns), size=20, color="#386CB0")
plot.dot(x=sigtonoise*noise_levs, y=noise_levs_rates, size=20, color="#386CB0")
# plot.y_range = Range1d(np.min(trad_nondet_surf(len(noise), num_skyrgns)),np.max(trad_nondet_surf(len(noise), num_skyrgns)))
# plot.x_range = Range1d(np.min(sigtonoise*max(noise)),np.max(sigtonoise*max(noise)))
for i in range(0,30,5):
    g = float(i)/10
    # print(g)
    omega = np.pi*extract_rad**2
    sens_range = np.linspace(0.00001,5,num=100)    
    nstar = -(np.log(1.0-conf_lev)/omega)*(sens_range/sigtonoise)**(-g)*(1/np.sum(np.power(noise,-g)))

    num_dens = nstar*np.power((sens_range/(10.*np.average(noise))),g)
    plot.line(sens_range, nstar)
    # plot.line(sens_range, num_dens)
    # y = -(np.log(1-conf_lev)/omega)*(10.*np.average(noise)/sigtonoise)**(-g)*(1/np.sum(np.power(noise,g)))
    # plot.line(noise, -(np.log(1.0-conf_lev)/omega)*(10.*np.average(noise)/sigtonoise)**(-g)*(1/np.sum(np.power(noise,g))),    line_width=2, line_color = "red")
plot.add_layout(Title(text="Noise (Jy/BM)", align="center"), "below")
plot.add_layout(Title(text="Transient Surface Density (1/deg^2))", align="center"), "left")
plot.toolbar.logo = None
plot.toolbar_location = None
plot.toolbar.active_drag = None
plot.toolbar.active_scroll = None
plot.toolbar.active_tap = None
output_file("tradsurfdens.html", title = "Traditional Transient Surface Density")
show(plot)






    



