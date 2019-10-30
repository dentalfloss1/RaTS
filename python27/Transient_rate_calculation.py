# USAGE:
# python transient_rate_cleaned.py [options] input_file

from random import random
import math
import os
import sys
import numpy
import getopt
import pyfits
import optparse
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib import rc
import pylab
rc('text', usetex=True)




usage = "USAGE: python transient_rate_cleaned.py [options] input_file"
parser = optparse.OptionParser()
parser.add_option("-s", "--single", dest = "single", help = "Calculates the transient rate upper limit for each field separately.", default=False, action="store_true")
parser.add_option("-c", "--combined", dest = "combined", help = "Calculates the transient rate upper limit for the combined images.", default=False, action="store_true")


(opts, args) = parser.parse_args()
images_noise = []
for i in range (0, len(args)):
	images_noise.append(args[i])

## exponent of the power law distribution of sources in flux: N(F) = k * (F/Fstar)**(-alpha)
alpha = [3.5, 3.0, 2.5, 2.0, 1.5, 1.0]

FOV = 15.48		## field of view analysed in deg**2

snapshot_time = 11./60/24/365	## amount of observing time per snapshot, in years
det_threshold = 8	## detection threshold used in the source finding step (in sigmas)

conf_lim = 0.05		## confidence level used in estimating the upper limit for k
my_x = []
my_y = []

def initialise():
	os.system('touch all_images')
	b = open('all_images', 'r+')
	for i in range (len(images_noise)):
		a = open(images_noise[i], 'r')
		a.seek(0)
		single_images = a.read().split('\n')
		for line in single_images:
			print >> b, line
	b.close()
	
	b = open('all_images', 'r')
	surv_per_field = per_field(b)
	unique_time_diff = surv_per_field[0]
	num_pairs = surv_per_field[1]

	for l in range(len(alpha)):
		selecting = read_select(b)
		values = calculation(selecting, alpha[l], unique_time_diff, num_pairs)
	a.close()
	b.close()
	os.system('rm all_images')

	if opts.single:
		for i in range (len(images_noise)):
			a = open(images_noise[i], 'r')
			a.seek(0)
			surv_per_field = per_field(a)
			unique_time_diff = surv_per_field[0]
			num_pairs = surv_per_field[1]

			for l in range(len(alpha)):
				selecting = read_select(a)
				values = calculation(selecting, alpha[l], unique_time_diff, num_pairs)
			a.close()

	my_x = values[0]
	my_y = values[1]
	old_xy = values[2]
	Fstar_new = values[3]
	tr_rate_perim = values[4]
	tr_rate_timescale = values[5]

	plot_comparison(my_x, my_y, old_xy, Fstar_new, tr_rate_perim, alpha, unique_time_diff, tr_rate_timescale)

def read_select(a):
	a.seek(0)
	noises = []
	images_ID = []
	bands = []
	discard_combined = []

	t = a.readlines()
	for line in t:
		image = line.split(' ')
		obsids = image[1].split('/')
		ids = obsids[-1].split('.')
		band = ids[0].split('_')

		noises.append(float(image[0]))		## reading the noise in each image
		discard_combined.append(band[0]+'_'+band[1])

	selecting = average_and_select(discard_combined, noises)
	return selecting

def calculation(selecting, alpha_value, unique_time_diff, pairs_tot):
	unique_noise = selecting[0]
	survivals = selecting[1]
	unique_noise_NEW = unique_noise[:]
	
	## calculating the total amount of square degrees observed and total observing time
	tot_sq_deg = (len(survivals)-4) * FOV
	tot_int_time = (len(survivals)-4) * snapshot_time
	
##############			NEW METHOD			#############	
	summation = 0.
	weight_noise = 0.

	## Calculating the total noise (sum and weighted sum)
	for i in range(len(unique_noise)):
		summation += unique_noise[i]
		weight_noise += (8 * unique_noise[i])**(1. - alpha_value)
	
	## measuring Fstar
	Fstar = 0.5

	## transient rate upper limit in # deg**(-2) and in # deg**(-2) s**(-1)
	tr_rate_deg_NEW = math.log(1./conf_lim) * (Fstar)**(1. - alpha_value) / FOV / weight_noise
	tr_rate_deg_time_NEW = math.log(1./conf_lim) * (Fstar)**(1. - alpha_value) / FOV / weight_noise / tot_int_time
	
	my_x.append(Fstar)
	my_y.append(tr_rate_deg_NEW)

##############			OLD METHOD (also new...)			#############
## transient rate upper limit in # deg**(-2)
	F0 = 0.
	for i in range(len(unique_noise)):
		if unique_noise[i] > F0:
			F0 = unique_noise[i]
			det = 8 * unique_noise[i]
		else: continue

	tr_rate_deg_OLD = - math.log(conf_lim) / tot_sq_deg
	tr_rate_deg_time_OLD= math.log(1./conf_lim) / tot_sq_deg / tot_int_time

	old_xy = [det, tr_rate_deg_OLD]
	
	Fstar_new = []
	tr_rate_perim = []
	for l in range(len(unique_noise_NEW)):
		max_noise = 0.
		for i in range(len(unique_noise_NEW)):
			if unique_noise_NEW[i] > max_noise:
				max_noise = unique_noise_NEW[i]
			else: continue
		det = 8. * max_noise
		Fstar_new.append(det)
		tr_rate_perim.append(math.log(1./conf_lim) / FOV / (len(unique_noise_NEW) - 4))
		unique_noise_NEW.remove(max_noise)
		if len(unique_noise_NEW) == 4: break
		else: continue

##############			NEW METHOD COMP		#############	
## New method but noise is just the highest in the dataset
	
	tr_rate_comp = math.log(1./conf_lim) * (det)**(1. - alpha_value) / FOV / weight_noise

##############			TIMESCALES IN PLACE		#############	
	tr_rate_timescale = []
	for i in range(len(unique_time_diff)):
		tr_rate_timescale.append(math.log(1./conf_lim) / FOV / pairs_tot[i])

	return my_x, my_y, old_xy, Fstar_new, tr_rate_perim, tr_rate_timescale


def get_DMJDs(survivals):
	snapshot_names = []
	for line in survivals:
		snap = line.split('_')
		snapshot_names.append(snap[0])

	DMJDs = []
	for line in snapshot_names:
		hdulist = pyfits.open('noises/'+line+'_SAP000_band3.fits') # open a FITS file
		hdr = hdulist[0].header
		date_obs = hdr['DATE-OBS']
		date_time = date_obs.split('T')
		split_date = date_time[0].split('-')
		split_hour = date_time[1].split(':')
		year = int(split_date[0])
		month = int(split_date[1])
		day = int(split_date[2])
		hour = int(split_hour[0])
		minute = int(split_hour[1])
		second = float(split_hour[2])

		if month == 1: M = 0
		elif month == 2: M = 31
		elif month == 3: M = 59
		elif month == 4: M = 90
		elif month == 5: M = 120
		elif month == 6: M = 150
		elif month == 7: M = 181
		elif month == 8: M = 212
		elif month == 9: M = 243
		elif month == 10: M = 273
		elif month == 11: M = 304
		elif month == 12: M = 334

		D = M + day
		## Day number
		
		sec_frac = second /60.
		min = minute + sec_frac
		min_frac = min / 60.
		h = hour + min_frac
		h_frac = h / 24.

		DMJD = D + h_frac
		DMJDs.append(DMJD)
		## DMJD = 0.0 on January 1st 2014

	unique_DMJDs = []
	[unique_DMJDs.append(item) for item in DMJDs if item not in unique_DMJDs]

	return unique_DMJDs, snapshot_names
	
def coalesce(input_data, bin_size):
	bin_start = input_data[0]
	list = []
	array = []
	for i in range(len(input_data)):
		if bin_start <= input_data[i] < bin_start + bin_size:
			list.append(input_data[i])
		else:
			array.append(list)
			bin_start = input_data[i]
			list=[input_data[i]]
	if list != []:
		array.append(list)
		
	return array

def number_of_pairs(timestamps, bin_size, timescale, tolerance):
	bins = coalesce(timestamps, bin_size)	
	pair_count = 0
	true_bins=[]

	for i in bins:
		if i != []:
			true_bins.append(i)
	for bin_number, my_bin in enumerate(true_bins[:-1]):
		if (timescale - tolerance) < (true_bins[bin_number+1][0] - my_bin[-1]) < (timescale + tolerance) or (timescale - tolerance) < (true_bins[bin_number+1][-1] - my_bin[0]) < (timescale + tolerance) or true_bins[bin_number+1][0] - my_bin[-1] <= timescale < true_bins[bin_number+1][-1] - my_bin[0]:
			pair_count += 1

	return pair_count


def per_field(a):
	## Getting the list of snapshots (different frequencies of the same snapshot merged together)
	selecting = read_select(a)
	survivals = selecting[1]

	## Getting the star time of each snapshot
	getDMJDs = get_DMJDs(survivals)
	DMJDs = getDMJDs[0]
	snapshot_names = getDMJDs[1]
	
	## Ordering the stat times
	index = range(len(DMJDs))
	index.sort(key = DMJDs.__getitem__)
	DMJDs[:] = [DMJDs[i] for i in index]
	snapshot_names[:] = [snapshot_names[i] for i in index]

	total_pairs_1snap = number_of_pairs(DMJDs, 0.01, 0.0104166666667, 0.01)
	total_pairs_2snap = number_of_pairs(DMJDs, 0.02, 0.0208333333449, 0.01)
	total_pairs_3snap = number_of_pairs(DMJDs, 0.03, 0.03125, 0.01)
	total_pairs_4snap = number_of_pairs(DMJDs, 0.04, 0.0416666666667, 0.01)
	total_pairs_5snap = number_of_pairs(DMJDs, 0.05, 0.0520833333449, 0.01)
	total_pairs_6snap = number_of_pairs(DMJDs, 0.06, 0.0625, 0.01)
	total_pairs_7snap = number_of_pairs(DMJDs, 0.07, 0.0729166666667, 0.01)
	total_pairs_1week = number_of_pairs(DMJDs, 1., 7., 2.)
	total_pairs_2weeks = number_of_pairs(DMJDs, 9., 14., 2.)
	total_pairs_1month = number_of_pairs(DMJDs, 16., 30., 12.)
	total_pairs_2months = number_of_pairs(DMJDs, 55., 60., 12.)
	
	num_pairs = []
	unique_time_diff = [
	0.0104166666667, #1snapshot, 15minutes
	0.0208333333449, #2snapshot, 30minutes
	0.03125, #3snapshot, 45minutes
	0.0416666666667, #4snapshot, 60minutes
	0.0520833333449, #5snapshot, 75minutes
	0.0625, #6snapshot, 90minutes
	0.0729166666667, #7snapshot, 105minutes
	7., #7days
	14., #14days
	30., #30days
	60.#60days
	]
	
	num_pairs.append(total_pairs_1snap)
	num_pairs.append(total_pairs_2snap)
	num_pairs.append(total_pairs_3snap)
	num_pairs.append(total_pairs_4snap)
	num_pairs.append(total_pairs_5snap)
	num_pairs.append(total_pairs_6snap)
	num_pairs.append(total_pairs_7snap)
	num_pairs.append(total_pairs_1week)
	num_pairs.append(total_pairs_2weeks)
	num_pairs.append(total_pairs_1month)
	num_pairs.append(total_pairs_2months)

	return unique_time_diff, num_pairs
	
def average_and_select(discard_combined, noises):
	## Selecting only unique snapshots (one snapshot = 1 obs; more bands are 1 snaphot!)
	IDslist = []
	[IDslist.append(item) for item in discard_combined if item not in IDslist]

	## Cutting out concatenated snapshots
	survivals = []
	for element in IDslist:
		if opts.combined:
			if element[0] == 'S':
				survivals.append(element)
		else:
			if element[0] != 'S':
				survivals.append(element)

		
	## Averaging noises for the same snapshot
	unique_noise = []
	for i in range(len(survivals)):
		noise_i = []
		for j in range(len(discard_combined)):
			if survivals[i] == discard_combined[j]:
				noise_i.append(noises[j])
			else: continue
		sum_noise = 0.
		for l in range(len(noise_i)):
			sum_noise += noise_i[l]
		unique_noise.append(sum_noise/len(noise_i))
	
	return unique_noise, survivals



# PLOTTING EVERYTHING NOW
def plot_comparison(my_x, my_y, old_xy, Fstar_new, tr_rate_perim, alpha, unique_time_diff, tr_rate_timescale):
	comp_freq_flND = [5.5, 2500.]
	comp_freq_rateND = [0.000075, 0.000000095]
	comp_freq_flTD = [0.014, 0.0021]
	comp_freq_rateTD = [0.013, 0.12]

# With Gal-Yam 2006 - GRB afterglows only
#	NOcomp_freq_flND = [0.008, 0.006, 0.040, 0.00009, 0.001, 0.01, 0.07, 3., 0.000500, 0.00037, 0.0002]
#	NOcomp_freq_rateND = [0.032, 0.0015, 0.004, 6., 1.,0.3, 0.003, 0.0009, 17., 0.6, 3.]

# Without Gal-Yam 2006 - GRB afterglows only
	NOcomp_freq_flND = [0.008, 0.040, 0.00009, 0.001, 0.01, 0.07, 3., 0.000500, 0.00037, 0.0002]
	NOcomp_freq_rateND = [0.032, 0.004, 6., 1.,0.3, 0.003, 0.0009, 17., 0.6, 3.]


	corr_flux_min_ND = [5.6, 1520.]
	corr_flux_min_TD = [0.047, 0.0049]
#	NOcorr_flux_min_ND = [0.038, 0.029, 0.19, 0.00102, 0.008, 0.083, 0.334, 14., 0.0057, 0.0042, 0.0023]
	NOcorr_flux_min_ND = [0.038, 0.19, 0.00102, 0.008, 0.083, 0.334, 14., 0.0057, 0.0042, 0.0023]


	corr_flux_min_VDL_ND = [5.3, 6129.]
	corr_flux_min_VDL_TD = [0.0016, 0.00046]
#	NOcorr_flux_min_VDL_ND = [0.00047, 0.00035, 0.0023, 0.0000011, 0.000021, 0.00021, 0.0041, 0.18, 0.000006, 0.0000046, 0.0000025]
	NOcorr_flux_min_VDL_ND = [0.00047, 0.0023, 0.0000011, 0.000021, 0.00021, 0.0041, 0.18, 0.000006, 0.0000046, 0.0000025]

	corr_flux_REFREP_ND = [5.8, 608.]
	corr_flux_REFREP_TD = [0.44, 0.010]
#	NOcorr_flux_REFREP_ND = [0.70, 0.52, 3.48, 0.092, 0.43, 4.3, 6.10, 261., 0.53, 0.38, 0.20]
	NOcorr_flux_REFREP_ND = [0.70, 3.48, 0.092, 0.43, 4.3, 6.10, 261., 0.53, 0.38, 0.20]



	xlabel = 'Flux Density S [Jy]'
	ylabel = 'Snapshot Surface Density [deg$^{-2}$]'

	if opts.combined:
		plotname = 'combined_upper_limit'
	else:
		plotname = 'snap_upper_limit'

	plt.figure()

	pylab.ylabel(r'{Snapshot Surface Density [deg$^{-2}$]', {'color':'k', 'fontsize':16})
	pylab.xlabel(r'{Flux Density S [Jy]', {'color':'k', 'fontsize':16})
	pylab.xticks(fontsize=18)
	pylab.yticks(fontsize=18)

	xvals = []
	yvals = []

	for i in NOcomp_freq_flND:
		xvals.append(i)
	for i in comp_freq_flND:
		xvals.append(i)
	for i in comp_freq_flTD:
		xvals.append(i)	
	for i in corr_flux_min_VDL_ND:
		xvals.append(i)
	for i in corr_flux_min_VDL_TD:
		xvals.append(i)
	for i in NOcorr_flux_min_VDL_ND:
		xvals.append(i)
	for i in corr_flux_min_ND:
		xvals.append(i)
	for i in corr_flux_min_TD:
		xvals.append(i)
	for i in NOcorr_flux_min_ND:
		xvals.append(i)
	for i in corr_flux_REFREP_ND:
		xvals.append(i)
	for i in corr_flux_REFREP_TD:
		xvals.append(i)
	for i in NOcorr_flux_REFREP_ND:
		xvals.append(i)
	for i in NOcomp_freq_rateND:
		yvals.append(i)
	for i in comp_freq_rateND:
		yvals.append(i)
	for i in comp_freq_rateTD:
		yvals.append(i)

	F = numpy.arange(min(xvals)*0.3, max(xvals)*2., 0.1)
	g=[]
	for i in range(len(alpha)):
		g.append(my_y[i]*(F / my_x[i])**(1 - alpha[i]))

	ymin=min(yvals)*0.3
	ymax=max(yvals)*2.
	xmin=min(xvals)*0.5
	xmax=max(xvals)*1.5
	ax = plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)

	for i in range(len(alpha)):
		ax.plot(F, g[i], color='g')

	ax.plot(F, g[0], color='g', label='model dependent')
	ax.plot(Fstar_new, tr_rate_perim,".", color='grey', label='model independent')

	yerr_up1=[0]*len(NOcomp_freq_flND)
	yerr_up2=[0]*len(comp_freq_flND)
	NOcomp_freq_rateND = numpy.double(NOcomp_freq_rateND)
	comp_freq_rateND = numpy.double(comp_freq_rateND)
	yerr_down1 = NOcomp_freq_rateND*0.6
	yerr_down2 = comp_freq_rateND*0.6
	lowerlimits1 = [5]*len(NOcomp_freq_flND)
	lowerlimits2 = [5]*len(comp_freq_flND)
	ax.errorbar(NOcomp_freq_flND, NOcomp_freq_rateND, yerr=[yerr_down1, yerr_up1], lolims=lowerlimits1, ecolor='k', fmt=None)
	ax.errorbar(comp_freq_flND, comp_freq_rateND, yerr=[yerr_down2, yerr_up2], lolims=lowerlimits2, ecolor='k', fmt=None)

	ax.plot(NOcomp_freq_flND, NOcomp_freq_rateND, 'o', mew=2., mec='red', mfc='None')
	ax.plot(comp_freq_flND, comp_freq_rateND, 'ro', label='Transient upper limits')

	ax.plot(comp_freq_flTD, comp_freq_rateTD, 'bD', label='Transient detection')

	ax.plot(my_x, my_y, '.', color='k')

	ax.text(1000, 0.0015, r'$\gamma$ = 0.0',fontsize=10)
	ax.text(1000, 0.00003, r'$\gamma$ = 0.5', fontsize=10)
	ax.text(1000, 0.0000006, r'$\gamma$ = 1.0', fontsize=10)
	ax.text(120, 0.0000002, r'$\gamma$ = 1.5', fontsize=10)
	ax.text(13, 0.0000002, r'$\gamma$ = 2.0', fontsize=10)
	ax.text(1.6, 0.0000002, r'$\gamma$ = 2.5', fontsize=10)
	ax.text(500, 10, r'Original', fontsize=13)
#	ax.text(0.003, 0.0008, r'\textbf{(}', fontsize=23)
#	ax.text(0.0075, 0.0008, r'\textbf{)}', fontsize=23)
		
	ax.errorbar(old_xy[0], old_xy[1], yerr=[[old_xy[1]*0.6], [0.]], lolims=5, ecolor='k', fmt=None)
	ax.plot(old_xy[0], old_xy[1], '*', color='r', ms=10)
	
	ax.tick_params(axis='both', direction='out')
#	ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
#	ax.get_yaxis().tick_left()
#	ax.legend(numpoints=1)

	plt.savefig(plotname+'_empty.eps',dpi=400)
#	plt.savefig(plotname+'_empty.png')
	plt.close()

	plt.figure()

	pylab.ylabel(r'{Snapshot Surface Density [deg$^{-2}$]', {'color':'k', 'fontsize':16})
	pylab.xlabel(r'{Flux Density S [Jy]', {'color':'k', 'fontsize':16})
	pylab.xticks(fontsize=18)
	pylab.yticks(fontsize=18)

	ax = plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)

	for i in range(len(alpha)):
		ax.plot(F, g[i], color='g')

	ax.plot(F, g[0], color='g', label='model dependent')
	ax.plot(Fstar_new, tr_rate_perim,".", color='grey', label='model independent')

	yerr_up1=[0]*len(NOcorr_flux_min_ND)
	yerr_up2=[0]*len(corr_flux_min_ND)
	NOcomp_freq_rateND = numpy.double(NOcomp_freq_rateND)
	comp_freq_rateND = numpy.double(comp_freq_rateND)
	yerr_down1 = NOcomp_freq_rateND*0.6
	yerr_down2 = comp_freq_rateND*0.6
	lowerlimits1 = [5]*len(NOcomp_freq_flND)
	lowerlimits2 = [5]*len(comp_freq_flND)
	ax.errorbar(NOcorr_flux_min_ND, NOcomp_freq_rateND, yerr=[yerr_down1, yerr_up1], lolims=lowerlimits1, ecolor='k', fmt=None)
	ax.errorbar(corr_flux_min_ND, comp_freq_rateND, yerr=[yerr_down2, yerr_up2], lolims=lowerlimits2, ecolor='k', fmt=None)

	ax.plot(NOcorr_flux_min_ND, NOcomp_freq_rateND, 'o', mew=2., mec='red', mfc='None')
	ax.plot(corr_flux_min_ND, comp_freq_rateND, 'ro', label='Transient upper limits')

	ax.plot(corr_flux_min_TD, comp_freq_rateTD, 'bD', label='Transient detection')

	ax.plot(my_x, my_y, '.', color='k')
	
	ax.text(1000, 0.0015, r'$\gamma$ = 0.0',fontsize=10)
	ax.text(1000, 0.00003, r'$\gamma$ = 0.5', fontsize=10)
	ax.text(1000, 0.0000006, r'$\gamma$ = 1.0', fontsize=10)
	ax.text(120, 0.0000002, r'$\gamma$ = 1.5', fontsize=10)
	ax.text(13, 0.0000002, r'$\gamma$ = 2.0', fontsize=10)
	ax.text(1.6, 0.0000002, r'$\gamma$ = 2.5', fontsize=10)
	ax.text(10, 10, r'Incoherent sources:', fontsize=13)
	ax.text(10, 5, r'slope $-$0.7', fontsize=13)
#	ax.text(0.015, 0.0008, r'\textbf{(}', fontsize=23)
#	ax.text(0.035, 0.0008, r'\textbf{)}', fontsize=23)

	ax.errorbar(old_xy[0], old_xy[1], yerr=[[old_xy[1]*0.6], [0.]], lolims=5, ecolor='k', fmt=None)
	plt.plot(old_xy[0], old_xy[1], '*', color='r', ms=10)

	ax.tick_params(axis='both', direction='out')
#	ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
#	ax.get_yaxis().tick_left()
#	ax.legend(numpoints=1)

	plt.savefig(plotname+'_corrFlux_empty.eps',dpi=400)
#	plt.savefig(plotname+'_corrFlux_empty.png')
	plt.close()

	plt.figure()

	pylab.ylabel(r'{Snapshot Surface Density [deg$^{-2}$]', {'color':'k', 'fontsize':16})
	pylab.xlabel(r'{Flux Density S [Jy]', {'color':'k', 'fontsize':16})
	pylab.xticks(fontsize=18)
	pylab.yticks(fontsize=18)

	ax = plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)

	for i in range(len(alpha)):
		ax.plot(F, g[i], color='g')

	ax.plot(F, g[0], color='g', label='model dependent')
	ax.plot(Fstar_new, tr_rate_perim,".", color='grey', label='model independent')


	yerr_up1=[0]*len(NOcorr_flux_min_VDL_ND)
	yerr_up2=[0]*len(corr_flux_min_VDL_ND)
	NOcomp_freq_rateND = numpy.double(NOcomp_freq_rateND)
	comp_freq_rateND = numpy.double(comp_freq_rateND)
	yerr_down1 = NOcomp_freq_rateND*0.6
	yerr_down2 = comp_freq_rateND*0.6
	lowerlimits1 = [5]*len(NOcomp_freq_flND)
	lowerlimits2 = [5]*len(comp_freq_flND)
	ax.errorbar(NOcorr_flux_min_VDL_ND, NOcomp_freq_rateND, yerr=[yerr_down1, yerr_up1], lolims=lowerlimits1, ecolor='k', fmt=None)
	ax.errorbar(corr_flux_min_VDL_ND, comp_freq_rateND, yerr=[yerr_down2, yerr_up2], lolims=lowerlimits2, ecolor='k', fmt=None)


	ax.plot(NOcorr_flux_min_VDL_ND, NOcomp_freq_rateND, 'o', mew=2., mec='red', mfc='None')
	ax.plot(corr_flux_min_VDL_ND, comp_freq_rateND, 'ro', label='Transient upper limits')

	ax.plot(corr_flux_min_VDL_TD, comp_freq_rateTD, 'bD', label='Transient detection')

	ax.plot(my_x, my_y, '.', color='k')
	
	ax.text(1000, 0.0015, r'$\gamma$ = 0.0',fontsize=10)
	ax.text(1000, 0.00003, r'$\gamma$ = 0.5', fontsize=10)
	ax.text(1000, 0.0000006, r'$\gamma$ = 1.0', fontsize=10)
	ax.text(120, 0.0000002, r'$\gamma$ = 1.5', fontsize=10)
	ax.text(13, 0.0000002, r'$\gamma$ = 2.0', fontsize=10)
	ax.text(1.6, 0.0000002, r'$\gamma$ = 2.5', fontsize=10)
	ax.text(0.5, 10, r'Expanding syncrotron bubble', fontsize=13)
#	ax.text(0.00018, 0.0008, r'\textbf{(}', fontsize=23)
#	ax.text(0.00042, 0.0008, r'\textbf{)}', fontsize=23)


	ax.errorbar(old_xy[0], old_xy[1], yerr=[[old_xy[1]*0.6], [0.]], lolims=5, ecolor='k', fmt=None)
	plt.plot(old_xy[0], old_xy[1], '*', color='r', ms=10)

	ax.tick_params(axis='both', direction='out')
#	ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
#	ax.get_yaxis().tick_left()
#	ax.legend(numpoints=1)

	plt.savefig(plotname+'_corrFlux_VDL_empty.eps',dpi=400)
#	plt.savefig(plotname+'_corrFlux_VDL_empty.png')
	plt.close()



	plt.figure()

	pylab.ylabel(r'{Snapshot Surface Density [deg$^{-2}$]', {'color':'k', 'fontsize':16})
	pylab.xlabel(r'{Flux Density S [Jy]', {'color':'k', 'fontsize':16})
	pylab.xticks(fontsize=18)
	pylab.yticks(fontsize=18)

	ymin=min(yvals)*0.3
	ymax=max(yvals)*2.
	xmin=min(xvals)*0.5
	xmax=max(xvals)*1.5

	ax = plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)

	for i in range(len(alpha)):
		ax.plot(F, g[i], color='g')

	ax.plot(F, g[0], color='g', label='model dependent')
	ax.plot(Fstar_new, tr_rate_perim,".", color='grey', label='model independent')

	yerr_up1=[0]*len(NOcorr_flux_REFREP_ND)
	yerr_up2=[0]*len(corr_flux_REFREP_ND)
	NOcomp_freq_rateND = numpy.double(NOcomp_freq_rateND)
	comp_freq_rateND = numpy.double(comp_freq_rateND)
	yerr_down1 = NOcomp_freq_rateND*0.6
	yerr_down2 = comp_freq_rateND*0.6
	lowerlimits1 = [5]*len(NOcorr_flux_REFREP_ND)
	lowerlimits2 = [5]*len(corr_flux_REFREP_ND)
	ax.errorbar(NOcorr_flux_REFREP_ND, NOcomp_freq_rateND, yerr=[yerr_down1, yerr_up1], lolims=lowerlimits1, ecolor='k', fmt=None)
	ax.errorbar(corr_flux_REFREP_ND, comp_freq_rateND, yerr=[yerr_down2, yerr_up2], lolims=lowerlimits2, ecolor='k', fmt=None)

	ax.plot(NOcorr_flux_REFREP_ND, NOcomp_freq_rateND, 'o', mew=2., mec='red', mfc='None')
	ax.plot(corr_flux_REFREP_ND, comp_freq_rateND, 'ro', label='Transient upper limits')
	ax.plot(corr_flux_REFREP_TD, comp_freq_rateTD, 'bD', label='Transient detection')

	ax.plot(my_x, my_y, '.', color='k')
	
	ax.text(1000, 0.0015, r'$\gamma$ = 0.0',fontsize=10)
	ax.text(1000, 0.00003, r'$\gamma$ = 0.5', fontsize=10)
	ax.text(1000, 0.0000006, r'$\gamma$ = 1.0', fontsize=10)
	ax.text(120, 0.0000002, r'$\gamma$ = 1.5', fontsize=10)
	ax.text(13, 0.0000002, r'$\gamma$ = 2.0', fontsize=10)
	ax.text(1.6, 0.0000002, r'$\gamma$ = 2.5', fontsize=10)
	ax.text(20, 10, r'Coherent sources:', fontsize=13)
	ax.text(20, 5, r'slope $-$2', fontsize=13)
#	ax.text(0.32, 0.0008, r'\textbf{(}', fontsize=23)
#	ax.text(0.62, 0.0008, r'\textbf{)}', fontsize=23)

	ax.errorbar(old_xy[0], old_xy[1], yerr=[[old_xy[1]*0.6], [0.]], lolims=5, ecolor='k', fmt=None)
	plt.plot(old_xy[0], old_xy[1], '*', color='r', ms=10)

	ax.tick_params(axis='both', direction='out')
#	ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
#	ax.get_yaxis().tick_left()
#	ax.legend(numpoints=1)

	plt.savefig(plotname+'_corrFlux_REFREPORT.eps',dpi=400)
	plt.close()


	xlabel = 'Timescale [days]'
	ylabel = 'Snapshot Surface Density [deg$^{-2}$]'

	plt.figure()
	pylab.ylabel(r'{Snapshot Surface Density [deg$^{-2}$]', {'color':'k', 'fontsize':16})
	pylab.xlabel(r'{Timescale [days]', {'color':'k', 'fontsize':16})
	pylab.xticks(fontsize=18)
	pylab.yticks(fontsize=18)

	xvals = []
	yvals = []

	## NDs
	x1 = numpy.arange(10./60./24., 365., 0.1)	#comp_freq
	x3 = numpy.arange(4.3, 45.3, 0.1)	#NOcomp_freq
	x5 = numpy.arange(81., 15.*365.25, 0.1)	#NOcomp_freq
	x11 = numpy.arange(1./24./60., 90., 0.1)	#NOcomp_freq

	##TDs
	x12 = numpy.arange(1., 365., 0.1)	#comp_freq
	x13 = numpy.arange(1., 3.*30., 0.1)	#comp_freq
	x14 = numpy.arange(20./24./60., 2.*7., 0.1)	#NOcomp_freq

	for i in comp_freq_rateND:
		yvals.append(i)
	for i in comp_freq_rateTD:
		yvals.append(i)
	for i in NOcomp_freq_rateND:
		yvals.append(i)

	for i in tr_rate_timescale:
		yvals.append(i)

	ymin = min(yvals) * 0.3
	ymax = max(yvals) * 2.
	xmin = 0.001
	xmax = 15. * 365.25

	ax = plt.gca()
	ax.set_yscale('log')
	ax.set_xscale('log')
	ax.set_xlim(xmin, xmax)
	ax.set_ylim(ymin, ymax)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	
	f1 = numpy.arange(min(unique_time_diff), max(unique_time_diff), 0.1)
	h1 = [old_xy[1]] * len(f1)
	ax.plot(f1, h1, color='grey', linewidth=3)

	g1 = [comp_freq_rateND[0]] * len(x1)
	g3 = [NOcomp_freq_rateND[0]] * len(x3)
	g5 = [NOcomp_freq_rateND[1]] * len(x5)
	g11 = [NOcomp_freq_rateND[7]] * len(x11)
	g12 = [comp_freq_rateTD[0]] * len(x12)
	g13 = [comp_freq_rateTD[1]] * len(x13)
	g14 = [NOcomp_freq_rateND[8]] * len(x14)

	ax.plot(x1, g1, color='r')
	ax.plot(x3, g3, color='r')
	ax.plot(x5, g5, color='r')
	ax.plot(x11, g11, color='r')
	ax.plot(x12, g12, color='b')
	ax.plot(x13, g13, color='b')
	ax.plot(x14, g14, color='r')

	timescales_ND = [5./60./24., 365.25, 30., 30., 1., 1., 26./60./24., 365.25, 2*30.]
	tr_rate_time_ND = [comp_freq_rateND[1], NOcomp_freq_rateND[2], NOcomp_freq_rateND[3], NOcomp_freq_rateND[4], NOcomp_freq_rateND[5], NOcomp_freq_rateND[6], comp_freq_rateND[0], comp_freq_rateND[0], NOcomp_freq_rateND[9]]

	yerr_up=[0]*len(tr_rate_time_ND)
	tr_rate_time_ND = numpy.double(tr_rate_time_ND)
	yerr_down = tr_rate_time_ND*0.6
	lowerlimits = [5]*len(tr_rate_time_ND)

	ax.errorbar(timescales_ND, tr_rate_time_ND, yerr=[yerr_down, yerr_up], lolims=lowerlimits, ecolor='k', fmt=None)
	ax.plot(timescales_ND, tr_rate_time_ND, 'ro', label='Transient upper limits')

#	ax.text(0.0018, 0.0012, r'1.4', fontsize=10)
	ax.text(0.0022, 0.00000013, r'0.074', fontsize=10)
	ax.text(0.52, 0.0007, r'1.4', fontsize=10)
	ax.text(0.52, 0.0025, r'1.4', fontsize=10)
	ax.text(16, 0.85, r'3.1', fontsize=10)
	ax.text(16, 0.25, r'3.1', fontsize=10)
	ax.text(35, 4, r'4.8/8.4', fontsize=10)
	ax.text(210, 8, r'4.8/8.4', fontsize=10)
	ax.text(1.5, 0.000045, r'0.15', fontsize=10)
	ax.text(0.2, 10, r'4.8', fontsize=10)
	ax.text(0.02, 0.7, r'4.8/8.4', fontsize=10)
	ax.text(1, 0.07, r'0.325', fontsize=10)
	ax.text(1, 0.0076, r'0.84', fontsize=10)
	ax.text(100, 0.0022, r'1.4', fontsize=10)
	ax.text(5, 0.017, r'1.4/4.8/8.4', fontsize=10)

	ax.plot(unique_time_diff, tr_rate_timescale,"o", color='grey', label='timescales')

	ax.tick_params(axis='both', direction='out')
#	ax.get_xaxis().tick_bottom()   # remove unneeded ticks 
#	ax.get_yaxis().tick_left()
#	ax.legend(numpoints=1)

	plt.savefig(plotname+'_timescale_moreBin.eps', dpi=400)
	plt.close()


if __name__ == "__main__":
	initialise()
	print 'done'
	exit()
