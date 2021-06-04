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
        
        return -(F0/((tau/2.0)**2))*(t - tau/2.0 - tcrit)**2 + F0
    
    def fluxint(self, F0, tcrit, tau, end_obs, start_obs):
        """Return the integrated flux"""
        
        tstart = np.maximum(tcrit , start_obs) 
        tend = np.minimum(tcrit + tau, end_obs) 
        fluxint = (F0*(tend-tstart) - (F0*(np.power((tend - tau/2.0 - tcrit),3.0)-np.power((tstart - tau/2.0 - tcrit),3.0))/(3.0*np.power((tau/2.0),2.0))))/(end_obs-start_obs)
        return fluxint
    
    def lines(self, xs, ys, durmax, max_distance, flux_err, obs):
        gaps = np.zeros(len(obs)-1,dtype=np.float32)
        for i in range(len(obs)-1):
            gaps[i] =  obs['start'][i+1] - obs['start'][i] + obs['duration'][i]
            # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
        min_sens = min(obs['sens'])
        max_sens = max(obs['sens'])
        sens_last = obs['sens'][-1]
        sens_maxgap = obs['sens'][np.where((gaps[:] == max(gaps)))[0]+1][0]
        dur_maxgap = obs['duration'][np.where((gaps[:] == max(gaps)))[0]+1][0]
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        day1_obs = obs['duration'][0]
        lastday_obs = obs['duration'][-1]
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        # maxdist_x = np.empty(len(ys))
        # maxdist_x.fill(np.log10(max_distance)) 
        # durmax_x = ' '
        maxdist_x = ' '
        durmax_y = np.zeros(xs.shape,dtype=np.float64)
        maxdist_y = np.zeros(xs.shape,dtype=np.float64)

# F0 = S_obs*(T_end - T_start)/((tend-tstart) - (np.power((tend - tau/2.0 - tcrit),3.0)-np.power((tstart - tau/2.0 - tcrit),3.0))/(3.0*np.power((tau/2.0),2.0)))
        for i in range(len(xs)):
            x = xs[i] 
            x = np.power(10,x)
            try:
                durmax_y[i] = (1 + flux_err) * sens_last * (lastday_obs) /(x - (np.power((durmax - lastday + x/2.0 ),3.0)-np.power((durmax - lastday_obs - x/2.0 ),3.0))/(3.0*np.power((x/2.0),2.0)))
                # durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs /  (day1_obs - (1.0/(3.0*np.power(x/2.0,2.0)))*(np.power(durmax + x, 3.0) - np.power(durmax - day1_obs + np.power(10.0, x), 3.0))))
            except:
                durmax_y[i] = np.inf
                # durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y[i] = (1 + flux_err) * sens_maxgap * dur_maxgap / (( x - max_distance) - (np.power((x/2.0),3.0)-np.power((max_distance - x/2.0),3.0))/(3.0*np.power((x/2.0),2.0)))
                # maxdist_y = np.append(maxdist_y, (1. + flux_err) * sens_maxgap * day1_obs /  (day1_obs - (1.0/(3*np.power(x/2.0,2.0)))*(np.power(max_distance + day1_obs, 3.0) - np.power(max_distance, 3.0))))
            except:
                maxdist_y[i] = np.inf
                # maxdist_y = np.append(maxdist_y, np.inf)    

        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) &  (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         