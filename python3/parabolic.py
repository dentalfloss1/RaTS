import numpy as np
class parabolic:
    """parabolic lightcurve class"""
    def __init__(self):
        self.edges=[1,1] # 1 is a definite edge, tophat is the default and has a definite beginning and end. Therefore it is [1,1]
        
    def earliest_crit_time(self, start_survey, tau):
        return start_survey - (tau/2.0)
    
    def latest_crit_time(self, end_survey, tau):
        return end_survey + (tau/2.0)
    
    def fluxint(self, F0, tcrit, tau, end_obs, start_obs):
        """Return the integrated flux"""
        tstart = np.maximum(tcrit - tau/2.0, start_obs) - (tcrit - tau/2.0)
        tend = np.minimum(tcrit + tau/2.0, end_obs) - (tcrit + tau/2.0)
        return (F0*(tend-tstart)/(end_obs-start_obs)) - ((F0/(3.0*(tau/2)**2(end_obs-start_obs)))*((tend - tau/2)**3 - (tstart - tau/2)**3))
    
    def lines(self, xs, ys, durmax, max_distance, flux_err, obs):
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        maxdist_x = np.empty(len(ys))
        maxdist_x.fill(np.log10(max_distance)) 
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs /  (day1_obs - (1.0/(3*np.power(np.power(10,x)/2.0,2.0)))*(np.power(durmax + np.power(10,x), 3.0) - np.power(durmax - day1_obs + np.power(10.0, x), 3.0))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y = np.append(maxdist_y, (1. + flux_err) * sens_maxgap * day1_obs /  (day1_obs - (1.0/(3*np.power(np.power(10,x)/2.0,2.0)))*(np.power(max_distance + day1_obs, 3.0) - np.power(max_distance, 3.0))))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)    
        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) &  (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         