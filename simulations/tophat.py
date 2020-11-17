import numpy as np
class tophat:
    """tophat lightcurve class"""
    def __init__(self):
        self.edges=[1,1] # 1 is a definite edge, tophat is the default and has a definite beginning and end. Therefore it is [1,1]
        
    def earliest_crit_time(self, start_survey, tau):
        return start_survey - tau
    
    def latest_crit_time(self, end_survey, tau):
        return end_survey
    
    def fluxint(self, F0, tcrit, tau, end_obs, start_obs):
        """Return the integrated flux"""
        tstart = np.maximum(tcrit, start_obs) - tcrit
        tend = np.minimum(tcrit + tau, end_obs) - tcrit
        return np.multiply(F0, np.divide((tend - tstart), (end_obs-start_obs)))
    
    def lines(self, xs, ys, durmax, max_distance, flux_err, obs):
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        maxdist_x = np.empty(len(ys))
        maxdist_x.fill(np.log10(max_distance))
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            durmax_y = np.append(durmax_y, np.inf)
            maxdist_y = np.append(maxdist_y, np.inf)
        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) &  (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         