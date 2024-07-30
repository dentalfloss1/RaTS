#!python
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib import colors, ticker
from tqdm import tqdm
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--flux",required=True,type=float, help="Flux at which to show transient rate")
parser.add_argument("--simfile",required=True,type=str, help="npz file from the transient simulations")
parser.add_argument("--outfile",required=True,type=str, help="outfile basename")
args = parser.parse_args()
plt.rcParams.update({
    "text.usetex": True})
files = []
font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 22}

matplotlib.rc('font', **font)
matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 

# files.append('outputmyrun59820_683119518384.npz')
# files.append('outputmyrun59820_68719579131.npz')
# files.append('outputmyrun59820_69226022373.npz'  )
# filenames = ['rgn1.png','rgn2.png','overlap.png']
f = args.simfile
fname = args.outfile



fig, axs = plt.subplots(1,2, sharex=True, sharey=False,figsize = (5,7.5))
fig.set_figwidth(20)
# fig.set_tight_layout(pad=)
counter = 0
for a in tqdm(axs):
    curfile = np.load(f,allow_pickle=True)
    fddetbins = curfile['fddetbins']
    fddethist = curfile['fddethist']
    senshist = curfile['senshist']
    X = curfile['X']
    Y = curfile['Y']
    Z = curfile['Z']
    ulZrate = curfile['ulZrate']
    vlinex = curfile['vlinex']
    xs = curfile['xs']
    ys = curfile['ys']
    durmax_x = curfile['durmax_x']
    durmax_y = curfile['durmax_y']
    durmax_y_indices = curfile['durmax_y_indices']
    maxdist_x = curfile['maxdist_x']
    maxdist_y_indices = curfile['maxdist_y_indices']
    day1_obs_x = curfile['day1_obs_x']
    ndecimal = int(np.round(-np.log10(args.flux)))+1
    whereindex = np.where(np.round(10**Y[:,0],decimals=ndecimal) == args.flux)
    index = whereindex[0][0]
    if counter==0:
        # https://matplotlib.org/stable/gallery/images_contours_and_fields/contourf_log.html#sphx-glr-gallery-images-contours-and-fields-contourf-log-py
        lev_exp = np.linspace(np.log10(ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].min()),np.log10(ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].max()), num=1000)
        # lev_exp = np.linspace(np.floor(np.log10(ulZrate[ulZrate > 0].min())),np.ceil(np.log10(np.mean(ulZrate[ulZrate > 0]))+1), num=1000)
        levs = np.power(10, lev_exp)
        # cs = ax.contourf(X, Y, z, levs, norm=colors.LogNorm())
        # levels = np.geomspace(max(np.amin(toplot[:,2]),1e-16),np.mean(toplot[:,2]),num = 1000)
        ulcsrate = a.contourf(10**X, 10**Y, ulZrate, levels=levs, cmap='viridis', norm=colors.LogNorm(), vmin=ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].min(), vmax=ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].max())
        ulrateticks = np.geomspace(ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].min(),ulZrate[(ulZrate > 0) & (ulZrate != np.inf)].max(),num=10)
        # ulrateticks = np.geomspace(ulZrate.min(),np.ceil(np.mean(ulZrate))+1,num=10)
        
        a.plot(10**xs, 10**np.full(xs.shape, np.log10(vlinex[0])),  color="red")
        # a.plot(10**np.full(ys.shape, np.log10(5/60/24)),10**ys,  color="red")
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_ylabel('Characteristic Flux (Jy)')
        a.set_xlim((10**np.amin(X),10**np.amax(X)))
        a.set_ylim(10**np.amin(Y),10**np.amax(Y))
        a.plot(np.geomspace(10**np.amin(X),10**np.amax(X),len(10**Y[index,:])),10**Y[index,:], ls='solid',color='red')
        cbarrate = fig.colorbar(ulcsrate, ticks=ulrateticks, format=ticker.StrMethodFormatter("{x:01.1e}"), ax=a)
        cbarrate.set_label('Transient Rate per day per sq. deg.')
        a.set_title('Transient Rate Upper Limits')
    elif counter==1:
        a.scatter(10**X[index,:], ulZrate[index,:] )
        a.set_ylim(1e-5,1e-2)
        a.set_xscale('log')
        a.set_yscale('log')
        a.set_title('Transient Rate at '+str(np.round(10**Y[index,0], decimals=3))+' Jy')
        a.set_ylabel('Transient Rate per day per sq. deg.')

    a.set_xlabel('Characteristic Duration (Days)')
    counter += 1 


plt.tight_layout()
plt.savefig(fname)
plt.close()
