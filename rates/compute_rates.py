from bokeh.plotting import figure, show, output_file
from bokeh.models import LinearColorMapper, SingleIntervalTicker, ColorBar, Title, Range1d
import random
from bokeh.io import export_png
import numpy as np 
import warnings
import configparser

from argparse import ArgumentParser
from datetime import datetime,timedelta
from scipy.special import binom
warnings.simplefilter("error", RuntimeWarning)

# Parse command line input
parser = ArgumentParser()
parser.add_argument('obsfile', help='supply a list of observations')
args = parser.parse_args()

observations = np.loadtxt(args.obsfile,dtype={'names': ('dateobs', 'duration'), 'formats': ('U20','f8')})

# Read config.ini settings into variables

params = configparser.ConfigParser()
params.read('config.ini')

conf_lev = np.float(params['STATISTICAL']['confidence'])
extract_rad = np.float(params['DATA']['extract_rad'])
sigtonoise = np.float(params['STATISTICAL']['sigtonoise'])
tsnap = np.float(params['DATA']['minint'])/60./60./24.
num_skyrgns = np.float(params['DATA']['num_skyrgns'])

# THE PARAMETERS COMMENTED OUT BELOW ARE FOR REPLICATING THE ORIGINAL PAPER ###############################################
# extract_rad = 5*np.sqrt(13./21/np.pi)
# num_skyrgns = 4
# The parameter below needs to change when considering a real, general setup. This is just replicating what the paper does
# sampletimescales = np.array([15, 30, 45, 60, 75, 90, 105, 7*24*60, 14*24*60, 30*24*60, 61*24*60])/60.0/24.0
# tsurvey = ((datetime.fromisoformat(observations[-1][0].replace('Z','+00:00')) + timedelta(minutes=observations[-1][1]*observations[-1][2])) - datetime.fromisoformat(observations[0][0].replace('Z','+00:00'))).total_seconds()/60./60./24.
# tsurvey = ((datetime.fromisoformat(observations[-1][0].replace('Z','+00:00')) + timedelta(minutes=observations[-1][1]*observations[-1][2])) - datetime.fromisoformat(observations[0][0].replace('Z','+00:00'))).total_seconds()/60./60./24.
# #########################################################################################################################
tsurvey = ((datetime.fromisoformat(observations[-1][0]) + timedelta(seconds=observations[-1][1])) - datetime.fromisoformat(observations[0][0])).total_seconds()/60./60./24.
sampletimescales = np.geomspace(tsnap, tsurvey, num=10)

# tgap = np.zeros((len(observations)-1,))
# totalimpersnap = np.array([obsdays//s for s in sampletimescales])
# npairs = np.array([binom(t,2) for t in totalimpersnap[totalimpersnap>2]])



# for s in sampletimescales:
    # i=0
    # while i<len(observations):
        # snap1 = datetime.fromisoformat(observations[i][0].replace('Z','+00:00'))
        # snapend = datetime.fromisoformat(observations[i][0].replace('Z','+00:00')) + timedelta(days=s)
        # setcond = 0
        # for j in range(1,len(observations)):
            # snapdiff = (datetime.fromisoformat(observations[j][0].replace('Z','+00:00'))) - snapend
            # print(snapdiff.total_seconds())
            # if snapdiff.total_seconds() > 0:
                # setcond = 1 # we need this to know if there are no additional iterations to go
                # tgap.append(snapdiff.total_seconds()/60./60./24.)
                # i=j
        
        # if setcond==0:
            # i = len(observations)
def npairsperT(T):
    T[T<tsnap] = tsnap
    timescalearr = []
    npairarr = []
    print(len(T))
    for t in T:
        # print(t/sampletimescales )
        # print(sampletimescales[(t//sampletimescales >= 1)][-1] )
        timescalearr.append(sampletimescales[(t//sampletimescales >= 1)][-1] )
    totalimpersnap = np.array([obsdays//t for t in timescalearr])
    return np.array([binom(totalimpersnap[i],2)*timescalearr[i] for i in range(len(totalimpersnap))])
    
def trad_nondet_surf(num_ims, num_skyrgns, num_dets):
    if num_dets>0:
        return num_dets/((num_ims-num_skyrgns)*np.pi*extract_rad**2)
    else:
        return  np.log(1.0-conf_lev)/(-(num_ims-num_skyrgns)*np.pi*extract_rad**2)
        #return  np.log(1.0-conf_lev)/(-2275)


def prob_gaps(tdur): # eqn 3.12 in Dario's thesis
    prob = np.zeros(tdur.shape,dtype=float)
    for i in range(len(tdur)):
        tau = tdur[i]
        numerator = tgap-tau
        numerator[numerator<0]=0
        prob[i] =  np.sum(numerator)/tsurvey
    return prob

def transrate(tdur): # eqn 3.15 in Dario's thesis

    omega = np.pi*extract_rad**2
    #The constant converts to 1/sky
    return -41252.96*np.log(1-conf_lev)/num_skyrgns/omega/(np.sum(tgap)+tdur)/(1-prob_gaps(tdur))
    
    
def transrateuncorr(T): # eqn 3.15 in Dario's thesis

    omega = np.pi*extract_rad**2
    #The constant converts to 1/sky
    print(npairsperT(T))
    return -41252.96*np.log(1-conf_lev)/num_skyrgns/omega/npairsperT(T)/tsnap
    
    
tdur = np.geomspace(start= 10**(-3), stop=10**(3), num=100)

print(transrateuncorr(tdur[npairsperT(tdur)>1]))
# exit()

# transrate(num)

rateplot = figure(title=" ", x_axis_type = "log", y_axis_type = "log" )

rateplot.cross(x=tdur, y=transrate(tdur), size=15, color="#386CB0", legend_label="Gap Corrected")

# print(transrate(tdur)/(1-prob_gaps(tdur)))
rateplot.diamond(x=tdur[npairsperT(tdur)>1], y=transrateuncorr(tdur[npairsperT(tdur)>1]), size=15, color="#b07c38", legend_label="Uncorrected")
rateplot.add_layout(Title(text="Duration (days)", align="center"), "below")
rateplot.add_layout(Title(text="Transient Rate (per sky, per day)", align="center"), "left")
rateplot.toolbar.logo = None
rateplot.toolbar_location = None
rateplot.toolbar.active_drag = None
rateplot.toolbar.active_scroll = None
rateplot.toolbar.active_tap = None
output_file("rateplot.html", title = "Transient Rate")
#export_png(p, filename=file + "_ProbContour.png")
show(rateplot)





noise, filenames = np.loadtxt("observations.txt", delimiter = ',', unpack=True, dtype=str)
# noise = [random.uniform(0.025,0.125) for l in range(151)]
#print(noise)
noise = np.array(noise)
noise = noise.astype(np.float)

noise_levs = np.array([max(noise)],dtype = np.float)
# print(len(noise))
for i in range(1,int(len(noise)-num_skyrgns)):
    noise_levs = np.append(noise_levs,max(noise[noise < noise_levs[i-1]])) 
noise_levs_rates = np.array([trad_nondet_surf(j, num_skyrgns, 0) for j in range(len(noise),num_skyrgns, -1)])

plot = figure(title=" ", x_axis_type = "log", y_axis_type = "log" )
plot.cross(x=8.0*max(noise), y=trad_nondet_surf(len(noise), num_skyrgns, 0), size=20, color="#386CB0")
plot.dot(x=8.0*noise_levs, y=noise_levs_rates, size=20, color="#386CB0")
plot.y_range = Range1d(np.min(trad_nondet_surf(len(noise)), num_skyrgns, 0),np.max(trad_nondet_surf(len(noise)))
plot.x_range = Range1d(np.min(8.0*max(noise)),np.max(8.0*max(noise)))
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






    
