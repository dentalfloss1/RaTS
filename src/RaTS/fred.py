import numpy as np
class fred:
    """fred lightcurve class"""
    def __init__(self):
        self.edges=[1,0] # 1 is a definite edge, tophat is the default and has a definite beginning and end. Therefore it is [1,1]
        
    def earliest_crit_time(self, start_survey, tau):
        return start_survey - tau
    
    def latest_crit_time(self, end_survey, tau):
        return end_survey
    
    def fluxint(self, F0, tcrit, tau, end_obs, start_obs):
        """Return the integrated flux"""
        tstart = np.maximum(tcrit, start_obs) - tcrit
        tend = end_obs - tcrit # Burst is never really "Off" so, how long is it on? It's on from tcrit until the end of the universe.
        exp1 = -np.divide(tstart,tau)
        exp1res = np.nan_to_num(np.exp(exp1))
        exp2 = -np.divide(tend,tau)
        exp2res = np.nan_to_num(np.exp(exp2))
        return np.multiply(F0, np.multiply(tau, np.divide(exp1res - exp2res, (end_obs-start_obs))))
        
    def lines(self, xs, ys, durmax, max_distance, flux_err, obs):
        gaps = np.zeros(len(obs)-1,dtype=np.float32)
        for i in range(len(obs)-1):
            gaps[i] =  obs['start'][i+1] - obs['start'][i] + obs['duration'][i]
            # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
        min_sens = min(obs['sens'])
        max_sens = max(obs['sens'])
        sens_last = obs['sens'][-1]
        sens_maxgap = obs['sens'][np.where((gaps[:] == max(gaps)))[0]+1][0]
        durmax_y = np.zeros(xs.shape,dtype=np.float64)
        maxdist_y = np.zeros(xs.shape,dtype=np.float64)
        lastdayobs = obs['duration'][-1]
        duration_maxgap = obs['duration'][np.where((gaps[:] == max(gaps)))[0]+1][0]
        for i in range(len(xs)):
            x = xs[i]
            try:
                durmax_y[i] =  (1. + flux_err) * sens_last * lastdayobs / np.power(10,x) / (np.exp(-(durmax - lastdayobs ) /  np.power(10,x)) - np.exp(-((durmax) / np.power(10,x))))
            except:
                durmax_y[i] =  np.inf
            try:
                maxdist_y[i] =   (((1. + flux_err) * sens_maxgap * duration_maxgap) /  np.power(10,x))   / (np.exp(-(max_distance / np.power(10,x))) - np.exp(-(max_distance + duration_maxgap) / np.power(10,x)))
            except:
                maxdist_y[i] =  np.inf
        durmax_x = ' '
        maxdist_x = ' '
        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) &  (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         