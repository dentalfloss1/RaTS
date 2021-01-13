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
parser.add_argument('obsfile', help='supply a list of observations')
args = parser.parse_args()

observations = np.loadtxt(args.obsfile,dtype={'names': ('dateobs', 'duration'), 'formats': ('U32','f8')})

if observations.size>1:
    observations = observations[np.argsort([datetime.fromisoformat(o) for o in observations['dateobs']])]
    
# Read config.ini settings into variables

params = configparser.ConfigParser()
params.read('config.ini')
conf_lev = np.float(params['STATISTICAL']['confidence'])
extract_rad = np.float(params['DATA']['extract_rad'])
sigtonoise = np.float(params['STATISTICAL']['sigtonoise'])
tsnap = np.float(params['DATA']['minint'])/60./60./24.
num_skyrgns = int(params['DATA']['num_skyrgns'])
detections = int(params['DATA']['detections'])


if observations.size>1:
    start_survey = min([datetime.fromisoformat(o) for o in observations['dateobs']])
    stop_survey = max([datetime.fromisoformat(o) for o in observations['dateobs']])
    stop_survey_ind = np.argmax([datetime.fromisoformat(o) for o in observations['dateobs']])
    tsurvey = ((stop_survey + timedelta(seconds=observations[stop_survey_ind][1])) - start_survey).total_seconds()/60/60/24
elif observations.size==1:
    start_survey = datetime.fromisoformat(str(observations['dateobs'])) 
    stop_survey = datetime.fromisoformat(str(observations['dateobs']))
    tsurvey = timedelta(seconds=np.float(observations['duration'])).total_seconds()/60/60/24
else:
    print("must have at least one observation")
    exit()
onsourcetime = np.sum(observations['duration'])/60/60/24
sampletimescales = np.geomspace(tsnap, tsurvey, num=10)
# #########################################################################################################################
# THE PARAMETERS COMMENTED OUT BELOW ARE FOR REPLICATING THE ORIGINAL PAPER ###############################################
extract_rad = 5*np.sqrt(13./21/np.pi)
num_skyrgns = 4
sampletimescales = np.array([15, 30, 45, 60, 75, 90, 105, 7*24*60, 14*24*60, 30*24*60, 61*24*60])/60.0/24.0
tsurvey = ((datetime.fromisoformat(observations[-1][0]) + timedelta(minutes=observations[-1][1])) - datetime.fromisoformat(observations[0][0])).total_seconds()/60./60./24.
# tsurvey = ((datetime.fromisoformat(observations[-1][0]) + timedelta(minutes=observations[-1][1])) - datetime.fromisoformat(observations[0][0])).total_seconds()/60./60./24.
onsourcetime = np.sum(observations['duration'])/60/24
observations['duration'] = observations['duration']/60/24
tsurvey = tsurvey*1.01 # perhaps a bit off with my estimation? 
detections=1
# #########################################################################################################################

if observations.size > 1:
    tgap = np.zeros((len(observations)-1,))
    tgap = np.array([(datetime.fromisoformat(observations[i+1][0]) - (datetime.fromisoformat(observations[i][0]) + timedelta(seconds=observations[i][1]))).total_seconds()/60/60/24 for i in range(len(observations)-1)])
elif observations.size == 1:
    tgap = 1e-12 # some small number


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
            
def upperlimitpoisson(alpha, k):
    return gammaincinv(k+1, 1-alpha/2)
def lowerlimitpoisson(alpha, k):
    return gammaincinv(k,alpha/2.)
def npairsperT(T):
    # print(T)
    timescalearr = []
    npairarr = []
    # print(len(T))
    # totalimpersnap = np.array([int(round(((onsourcetime/t)-1))) for t in T],dtype='i4')
    # totalimpersnap = np.array([binom(onsourcetime/t,2) for t in T],dtype='i4')
    # imhist = totalimpersnap
    imhist = np.zeros(len(T),dtype='i4')
    # imhistcmp = np.zeros(len(T),dtype='i4')
    for i in range(len(T)):
        t = T[i]
        # startbin = datetime.fromisoformat(observations['dateobs'][0])
        # stopbin = datetime.fromisoformat(observations['dateobs'][-1]) + timedelta(days=observations['duration'][-1]) + timedelta(days=t)
        # totalbins = int(round((stopbin-startbin).total_seconds()/timedelta(days=t).total_seconds()))
        # imhist[i]=totalbins
        # print(t,tgap)
        # if t<np.all(tgap):
            # imhist[i] = int(round(onsourcetime/t))
            # print('cond1')
        # else:
        # maximages = int(round(onsourcetime/tsnap))
        # if t<=tsnap:
            # imhist[i] = maximages
        # else:
            # for o in observations:
                # if (t<=(timedelta(days=o['duration']).total_seconds/60/60/24)):
                    # imhist[i]+=int(max(1,(timedelta(days=o['duration']).total_seconds/60/60/24)/t))
                # elif (t>(timedelta(days=o['duration']).total_seconds/60/60/24)):
                    
        startbin = datetime.fromisoformat(observations['dateobs'][0])
        stopbin = datetime.fromisoformat(observations['dateobs'][-1]) + timedelta(days=observations['duration'][-1]) + timedelta(days=t)
        totalbins = int(round((stopbin-startbin).total_seconds()/timedelta(days=t).total_seconds()))
        # print(totalbins)
        if t<min(observations['duration']):
            imhist[i] = int(round(onsourcetime/t))
        else:
            for j in range(totalbins):
                # datelist = [datetime.fromisoformat(o) for o in observations['dateobs']]
                # durlist = [timedelta(days=o) for o in observations['duration']]
                localbinL = (startbin + j*timedelta(days=t))
                localbinR = (startbin + (j+1)*timedelta(days=t))
                # print(localbinL,localbinR)
                # tfmaskarray = [(localbinL<=d) and (localbinR>=d) for d in datelist]
                # print(observations)
                # print(observations['dateobs'])
                # print(observations['dateobs'].shape)
                obsdates = [datetime.fromisoformat(d) for d in observations['dateobs']]
                obsdurs = [timedelta(days=d) for d in observations['duration']]
                # print(np.sum([localbinL<=dat+dur<localbinR for (dat,dur) in zip(obsdates,obsdurs)]))
                # input("presskey")
                # if np.sum([localbinL<=dat+dur<localbinR for (dat,dur) in zip(obsdates,obsdurs)])!=0:
                    # print(np.sum([localbinL<=dat+dur<localbinR for (dat,dur) in zip(obsdates,obsdurs)]))
                # imhist[i] = imhist[i] -  (np.sum([localbinL<=dat+dur<localbinR for (dat,dur) in zip(obsdates,obsdurs)])>1) 
                # if t>5e-02:
                    # print((np.sum([localbinL<=dat+dur<localbinR for (dat,dur) in zip(obsdates,obsdurs)])>1) )
                    # print(imhistcmp[i], imhist[i], int(np.sum([localbinL<=dat+dur<localbinR for (dat,dur) in zip(obsdates,obsdurs)])>1) )
                
                for (date, dur) in zip(observations['dateobs'], observations['duration']):
                    # print(date,dur)
                    startobs = datetime.fromisoformat(date)
                    endobs = datetime.fromisoformat(date) + timedelta(days=dur)
                    # print(min(localbinR,endobs) - max(localbinL, startobs))
                    # input("presskey")

                    leftcond = localbinL<=endobs
                    # print(localbinL,endobs)
                    # input("presskey")
                    # print("------------------")
                    observedinbin = max(0,(min(localbinR,endobs) - max(localbinL, startobs)).total_seconds())/60/60/24
                    if (max(localbinL, startobs) < min(localbinR,endobs)) and observedinbin>tsnap:
                        imhist[i]+=1
                        # if t>(104/60/24):
                            # print(max(localbinL, startobs), min(localbinR,endobs), observedinbin)
                        # imhistcmp[i]+=1
                        break
                        # print((min(localbinR,endobs) - max(localbinL, startobs)))
                        # print(max(localbinL, startobs))
                        # print(min(localbinR,endobs))

                        # imhistflt[i] += observedinbin/t
                        # print(max(0,(min(localbinR,endobs) - max(localbinL, startobs)).total_seconds()))
                        # print(t)
                        # print(imhistflt[i])
                        # input("presskey")
                        # print(j,localbinL,startobs,max(localbinL, startobs))
                        # print((datetime.fromisoformat(o['dateobs'])+timedelta(days=o['duration'])-timedelta(days=tsnap)),t,timedelta(days=o['duration']).total_seconds())
                        # print(imhistflt)
                        # input("presskey")
            # imhist[i] = int(round(imhistflt[i]))
                # anydetections = int(np.any([(localbinL<=(datetime.fromisoformat(o['dateobs'])+timedelta(days=o['duration'])-timedelta(days=tsnap)) and (localbinR>=datetime.fromisoformat(o['dateobs']))) for o in observations]))
                # print(anydetections)
                # print(localbinL, localbinR)
                # print([d for d in observations['dateobs']])
                # print([(localbinL<=(datetime.fromisoformat(o['dateobs'])+timedelta(days=o['duration'])) - timedelta(days=tsnap)) for o in observations])
                # print([(localbinR>=datetime.fromisoformat(o['dateobs'])) for o in observations])
                # input("presskey")
                # tfmask = np.multiply(tfmaskarray, [max(1,int(round(d.total_seconds()/60/60/24/t))) for d in durlist])
                # print([max(1,int(round(d.total_seconds()/60/60/24/t))) for d in durlist])
                # print(tfmask)
                # if np.any(tfmask!=(tfmask*0)):
                    # print(tfmask)
                # import pdb
                # pdb.set_trace()
                # if i >4:
                    # print(tfmask)
                # input("press key")
                # imhistflt[i]+=anydetections
                # print(imhist[i])
                    # input("press key")
                # imhist[i] = np.sum(tfmask)
                # print(totalbins)
                
            # for j in range(int(round(tsurvey/t))):
                # startbin = datetime.fromisoformat(observations['dateobs'][0]) + j*timedelta(days=t)
                # stopbin = datetime.fromisoformat(observations['dateobs'][0]) + (j+1)*timedelta(days=t)
                # for o in observations:
                    # if (startbin<=datetime.fromisoformat(o[0])) and (stopbin>=datetime.fromisoformat(o[0])):
                        # imhist[i] += int(round(o[1]/t))
                        
            # print(imhist)
            # imhist = np.zeros(int(round(tsurvey/t)))
            # datebins = [(datetime.fromisoformat(observations['dateobs'][0]) + i*timedelta(days=t)) for i in range(len(imhist)+1)]
            # print(datebins
            # maskedarray = np.array([d<datebins[3] for d in datebins])
            # print(datebins[np.where(maskedarray)[0]])
            # print(datebins[np.wheremaskedarray)[0][:]])
            # print(datebins[maskedarray])
        # print(sampletimescales)
        # print(imhistflt)
        # print(imhist)
        # exit()
            # for i in range(len(round(tsurvey/t))):
                # lower = datetime.fromisoformat(observations['dateobs'][i])
                # upper = 
            # imhist[observations['dateobs']
    # print(T)
    # print(imhist)
    # print(imhistcmp)
    imhist-=1
    # print([int(t) for t in totalimpersnap])
    # exit()
    

    return imhist
    
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
        # print(tgap, tau)
        # input("press key")
        if numerator.size>1:
            numerator[numerator<0]=0
        else:
            numerator = 0
        # print(numerator)
        # print(np.sum(numerator)/tsurvey)
        prob[i] =  np.sum(numerator)/tsurvey
    # print(prob)
    return prob

def transrate(T): # eqn 3.15 in Dario's thesis

    omega = np.pi*extract_rad**2
    #The constant converts to 1/sky
    if detections==0:
        return -41252.96*np.log(1-conf_lev)/num_skyrgns/omega/(tsurvey+T)/(1-prob_gaps(T))
    else:
        return 4152.96*lowerlimitpoisson(1-conf_lev, detections)/num_skyrgns/omega/(tsurvey+T)/(1-prob_gaps(T)), 4152.96*upperlimitpoisson(1-conf_lev,detections)/num_skyrgns/omega/(tsurvey+T)/(1-prob_gaps(T))
    
    
def transrateuncorr(T): # eqn 3.15 in Dario's thesis

    omega = np.pi*extract_rad**2
    #The constant converts to 1/sky
    # print(npairsperT(T))
    if detections==0:
        return -41252.96*np.log(1-conf_lev)/num_skyrgns/omega/npairsperT(T)/sampletimescales[npairsperT(sampletimescales)>1]
    else:
        return 4152.96*lowerlimitpoisson(1-conf_lev, detections)/omega/num_skyrgns/npairsperT(T)/sampletimescales[npairsperT(sampletimescales)>1], 4152.96*upperlimitpoisson(1-conf_lev,detections)/omega/num_skyrgns/npairsperT(T)/sampletimescales[npairsperT(sampletimescales)>1]
    
    
# tdur = np.geomspace(start= tsnap/10, stop=tsurvey*100, num=100)
tdur = sampletimescales
# print(tdur)

# fig = plt.figure()
# plt.scatter(x=tdur, y=transrate(tdur), marker='X', label='Gap Corrected')
# plt.scatter(x=tdur[npairsperT(tdur)>1], y=transrateuncorr(tdur[npairsperT(tdur)>1]), marker='D', label='Uncorrected')
# ax = plt.gca()
# ax.legend()
# plt.xscale('log')
# plt.yscale('log')

# plt.show()
# print(tdur)
# print(transrate(tdur))
if detections == 0:
    rateplot = figure(title=" ", x_axis_type = "log", y_axis_type = "log" )
    rateplot.cross(x=tdur, y=transrate(tdur), size=15, color="#386CB0", legend_label="Gap Corrected")
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
else:
# 
    rateplot = figure(title=" " , x_axis_type = "log", y_axis_type = "log")
    # rateplot.cross(x=tdur, y=transrate(tdur), size=15, color="#386CB0", legend_label="Gap Corrected")
    uncorrlower, uncorrupper = transrateuncorr(tdur[npairsperT(tdur)>1])
    corrlower, corrupper = transrate(tdur)
    # print(tdur[npairsperT(tdur)>1])
    # print(uncorrlower,uncorrupper)
    # data = {'tdurcorr': tdur,
    # 'tduruncorr': tdur[npairsperT(tdur)>1],
    # 'uncorrlower': uncorrlower,
    # 'uncorrupper': uncorrupper,
    # 'corrlower': corrlower,
     # 'corrupper': corrupper,
     # 'width': tdur*5e-02}

    # source = ColumnDataSource(data=data)    
    # rateplot.vbar_stack(['uncorrupper', 'corrupper'],width='width', x=['tduruncorr','tdurcorr'], color=['blue', 'red'], source=source)
    rateplot.vbar(x=tdur+5e-2*tdur, width=5e-2*tdur,bottom=uncorrlower, top=uncorrupper, color="#b07c38", legend_label="Uncorrected")
    rateplot.vbar(x=tdur[npairsperT(tdur)>1]-5e-2*tdur[npairsperT(tdur)>1], width=5e-2*tdur[npairsperT(tdur)>1],bottom=corrlower, top=corrupper, color="#386CB0", legend_label="Corrected")
    # rateplot.vbar(bottom=stack('uncorrlower'), top=stack('uncorrupper'), x=stack('tduruncorr'), width=5e-2*tdur, color='#b07c38', source=source, name='uncorrlower')
    # rateplot.vbar(bottom=stack('uncorrlower'), top=stack('corrupper'), x=stack('tduruncorr'), width=5e-2*tdur, color='#b07c38', source=source, name='corrlower')
    # rateplot.vbar(bottom=stack('2016'), top=stack('2016', '2017'), x=10, width=0.9, color='red',  source=source, name='2017')
    # rateplot.vbar_stack(x=['tduruncorr','tdurcorr'], width=5e-2*tdur,bottom=['uncorrlower','corrlower'], top=['uncorrupper','corrupper'], color=["#386CB0","#b07c38"], source=source)
    rateplot.y_range = Range1d(np.min(np.concatenate((uncorrlower,corrlower)))*0.9,np.max(np.concatenate((corrupper,uncorrupper)))*1.1)
    # rateplot.y_range = Range1d(1e-3,1e4)
    # rateplot.y_range = Range1d(1e-3,1e4)
    rateplot.x_range = Range1d(np.min(tdur)*0.9,np.max(tdur)*1.1)
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

exit()
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
# plot.y_range = Range1d(np.min(trad_nondet_surf(len(noise), num_skyrgns, 0)),np.max(trad_nondet_surf(len(noise), num_skyrgns, 0)))
# plot.x_range = Range1d(np.min(8.0*max(noise)),np.max(8.0*max(noise)))
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






    
