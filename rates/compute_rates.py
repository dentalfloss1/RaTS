from bokeh.plotting import figure, show, output_file
from bokeh.models import LinearColorMapper, SingleIntervalTicker, ColorBar, Title, Range1d
import random
from bokeh.io import export_png
import numpy as np 
import warnings
warnings.simplefilter("error", RuntimeWarning)

conf_lev = 0.95
# extract_rad = 1.39
extract_rad = 5*np.sqrt(13./21/np.pi)
sigtonoise = 8
# num_skyrgns = 1
tsnap = 15/60/24. # 15 minute snapshot into seconds
num_skyrgns = 4
# The parameter below needs to change when considering a real, general setup. This is just replicating what the paper does
sampletimescales = np.array([15, 30, 45, 60, 75, 90, 105, 7*24*60, 14*24*60, 30*24*60, 61*24*60])/60.0/24.0
tsurvey = max(sampletimescales)
tgap = sampletimescales


def trad_nondet_surf(num_ims, num_skyrgns, num_dets):
    if num_dets>0:
        return num_dets/((num_ims-num_skyrgns)*np.pi*extract_rad**2)
    else:
        return  np.log(1.0-conf_lev)/(-(num_ims-num_skyrgns)*np.pi*extract_rad**2)
        #return  np.log(1.0-conf_lev)/(-2275)


def prob_gaps(tdur, tgap): # eqn 3.12 in Dario's thesis
    prob = np.array([],dtype=float)
    for tau in tdur:
        numerator = tgap-tau
        numerator[numerator<0]=0
        prob = np.append(prob, np.sum(numerator/np.sum(tgap)))
    return prob

def transrate(tdur): # eqn 3.15 in Dario's thesis

    omega = np.pi*extract_rad**2
    #The constant converts to 1/sky
    return -41252.96*np.log(1-conf_lev)/num_skyrgns/omega/(np.sum(tgap)+tdur)/(1-prob_gaps(tdur,tgap))
    
tdur = np.geomspace(start= 10**(-3), stop=10**(3), num=100)

# transrate(num)


rateplot = figure(title=" ", x_axis_type = "log", y_axis_type = "log" )
rateplot.cross(x=tdur, y=transrate(tdur), size=15, color="#386CB0", legend_label="Gap Corrected")
rateplot.diamond(x=tdur, y=np.multiply(transrate(tdur),(1-prob_gaps(tdur,tgap))), size=15, color="#b07c38", legend_label="Uncorrected")
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
print(len(noise))
for i in range(1,int(len(noise)-num_skyrgns)):
    noise_levs = np.append(noise_levs,max(noise[noise < noise_levs[i-1]])) 
noise_levs_rates = np.array([trad_nondet_surf(j, num_skyrgns, 0) for j in range(len(noise),num_skyrgns, -1)])

plot = figure(title=" ", x_axis_type = "log", y_axis_type = "log" )
plot.cross(x=8.0*max(noise), y=trad_nondet_surf(len(noise), num_skyrgns, 0), size=20, color="#386CB0")
plot.dot(x=8.0*noise_levs, y=noise_levs_rates, size=20, color="#386CB0")
plot.y_range = Range1d(1e-5,1e-1)
plot.x_range = Range1d(1e-4,5e-1)
for i in range(0,30,5):
    g = float(i)/10
    print(g)
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






    
