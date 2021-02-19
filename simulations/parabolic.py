import numpy as np
# from tqdm import tqdm
# from scipy.integrate import quad
# from scipy.integrate import trapz
class parabolic:
    """parabolic lightcurve class"""
    def __init__(self):
        self.edges=[1,1] # 1 is a definite edge, tophat is the default and has a definite beginning and end. Therefore it is [1,1]
        
    def earliest_crit_time(self, start_survey, tau): # For lightcurves with definite edges, tcrit MUST be the beginning
        return start_survey - tau
    
    def latest_crit_time(self, end_survey, tau):
        return end_survey 
    
    def lightcurve(self, t, F0, tau, tcrit):
        
        return -(F0/((tau/2.0)**2))*(t + tau/2.0 - tcrit)**2 + F0
    
    def fluxint(self, F0, tcrit, tau, end_obs, start_obs):
        """Return the integrated flux"""
        # print("Begin obs calc")
        tstart = np.maximum(tcrit , start_obs) 
        tend = np.minimum(tcrit + tau, end_obs) 
        # I = np.array([],dtype=np.float64)
        # Ierr = np.array([],dtype=np.float64)
        # fluxint = np.array([], dtype=np.float)
        # for i in tqdm(range(len(tcrit))):
            # tmp = quad(self.lightcurve, tstart[i], tend[i], args=(F0[i], tau[i], tcrit[i]))
            # I = np.append(I,tmp[0]/(end_obs-start_obs))
            # Ierr = np.append(Ierr,tmp[1])
            # fluxint = np.append(fluxint,(F0[i]*(tend[i]-tstart[i]) - ((F0[i]/(3.0*(tau[i]/2.0)**2.0))*((tend[i]-tcrit[i])**3.0-(tstart[i]-tcrit[i])**3.0)))/(end_obs-start_obs))
            # if fluxint[i] < 0 :
                # print("Negative Value: ",fluxint[i])
        # lin = F0*(tend-tstart)
        # numer = (F0*(tend-tstart) - (F0/(3.0*np.power((tau/2.0),2.0))*(np.power((tend - tau/2.0 - tcrit),3.0)-np.power((tstart - tau/2.0 - tcrit),3.0))))
        # denom = (end_obs-start_obs)
        fluxint = (F0*(tend-tstart) - (F0/(3.0*np.power((tau/2.0),2.0))*(np.power((tend - tau/2.0 - tcrit),3.0)-np.power((tstart - tau/2.0 - tcrit),3.0))))/(end_obs-start_obs)
        # negvals = np.where(fluxint<0)[0]
        # negvals_num = np.where(numer<0)[0]
        # negvals_denom = np.where(denom<0)[0]
        # negvals_lin = np.where(lin<0)[0]
        # for i in range(len(negvals)):
            # print(fluxint[negvals[i]])
        # print(end_obs, start_obs)
        # print(np.amin(fluxint))
        # print(np.amin(lin))
        # print(tend[np.where(np.amin(lin)==lin)[0]],tstart[np.where(np.amin(lin)==lin)[0]],F0[np.where(np.amin(lin)==lin)[0]]*(tend[np.where(np.amin(lin)==lin)[0]] - tstart[np.where(np.amin(lin)==lin)[0]]),F0[np.where(np.amin(lin)==lin)[0]]  )
                #print(tcrit[negvals_lin[i]],tau[negvals_lin[i]], end_obs, start_obs, end_obs-start_obs, tstart[negvals_lin[i]],tend[negvals_lin[i]])
        # print("Max Error is: ",np.amax(Ierr))
        # print("Max diff is", np.amax(fluxint - I))
        return fluxint
    
    def lines(self, xs, ys, durmax, max_distance, flux_err, obs):
        gaps = np.array([],dtype=np.float32)
        for i in range(len(obs)-1):
            gaps = np.append(gaps, obs['start'][i+1] - obs['start'][i] + obs['duration'][i,1])
            # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
        min_sens = min(obs['sens'])
        max_sens = max(obs['sens'])
        sens_last = obs['sens'][-1]
        sens_maxgap = obs['sens'][np.where((gaps[:] == max(gaps)))[0]+1][0]
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        day1_obs = obs['duration'][0]
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        maxdist_x = np.empty(len(ys))
        maxdist_x.fill(np.log10(max_distance)) 
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs /  (day1_obs - (1.0/(3.0*np.power(np.power(10,x)/2.0,2.0)))*(np.power(durmax + np.power(10,x), 3.0) - np.power(durmax - day1_obs + np.power(10.0, x), 3.0))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y = np.append(maxdist_y, (1. + flux_err) * sens_maxgap * day1_obs /  (day1_obs - (1.0/(3*np.power(np.power(10,x)/2.0,2.0)))*(np.power(max_distance + day1_obs, 3.0) - np.power(max_distance, 3.0))))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)    
        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) &  (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         