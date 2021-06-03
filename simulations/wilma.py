import numpy as np
class wilma:
    """wilma lightcurve class"""
    def __init__(self):
        self.edges=[0,1] # 1 is a definite edge, tophat is the default and has a definite beginning and end. Therefore it is [1,1]
        
    def earliest_crit_time(self, start_survey, tau):
        return start_survey 
    
    def latest_crit_time(self, end_survey, tau):
        return end_survey + tau
    
    def fluxint(self, F0, tcrit, tau, end_obs, start_obs):
        """Return the integrated flux"""
        tstart = start_obs - tcrit# Burst really starts, so we always start at the beginning of the observation
        tend = np.minimum(end_obs,tcrit) - tcrit 
    
        return np.multiply(F0, np.multiply(tau, np.divide(np.exp(np.divide(tend,tau)) - np.exp(np.divide(tstart,tau)) , (end_obs-start_obs))))
        
    def lines(self, xs, ys, durmax, max_distance, flux_err, obs):
        gaps = np.array([],dtype=np.float32)
        for i in range(len(obs)-1):
            gaps = np.append(gaps, obs['start'][i+1] - obs['start'][i] + obs['duration'][i])
            # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
        min_sens = min(obs['sens'])
        max_sens = max(obs['sens'])
        sens_last = obs['sens'][-1]
        sens_maxgap = obs['sens'][np.where((gaps[:] == max(gaps)))[0]-1][0]
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        day1_obs = obs['duration'][0]
        obspregap = obs['duration'][np.where((gaps[:] == max(gaps)))[0]-1][0]
        sens_first = obs['sens'][0]
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs / np.power(10,x) / (np.exp(-(durmax - day1_obs) /  np.power(10,x)) - np.exp(-((durmax) / np.power(10,x)))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, (((1. + flux_err) * sens_maxgap * obspregap) /  np.power(10,x))   / (np.exp(-(max_distance / np.power(10,x))) - np.exp(-(max_distance + obspregap) / np.power(10,x))))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)    
        durmax_x = ' '
        maxdist_x = ' '
        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) &  (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         