configfilestring="""[INITIAL PARAMETERS]
; integer number of sources to be simulated
srcperbin = 100
; Minimum simulated flux, in same units as the flux in the observations file ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
fl_min = 5e-5
; Maximum simulated flux ,in same units as the flux in the observations file  
fl_max = 5  
; Fractional error in the simulated flux ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
flux_err = 1e-1  
; Minimum simulated duration, in same units as the duration in the observations file  
dmin = 5e-2 
; Maximum simulated duration, in same units as the duration in the observations file     
dmax = 5e3
; detection threshold, to be multiplied by the noise in the images
det_threshold = 5 
; integer, extra detection threshold
extra_threshold = 0
; Name to be used for the output files
file = myrun
; Must be present as a python file 
lightcurvetype = tophat
; only used for choppedgaussian
gaussiancutoff = 0.1
; confidence levels for limits
confidence = 95
; Upper limit on detections. 3 is 95% confidence
detections = 0

; The following is only used if no observations are provided in an observations file
[SIM]
; Number of observations to be simulated 
nobs =  46
; Average sensitivity of simulated observations in Jy 
obssens = 21.7e-6
; std dev of simulated sens in Jy 
obssig = 4.6e-6
; interval between simulated observations in days 
obsinterval = 7
; duration of observations in days 
obsdurations =  0.009"""
