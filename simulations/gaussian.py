import numpy as np
import sys
from scipy.stats import norm
from scipy.special import erf
class gaussian:
    """gaussian lightcurve class"""
    def __init__(self):
        self.edges=[0,0] # 1 is a definite edge, tophat is the default and has a definite beginning and end. Therefore it is [1,1]
        
    def earliest_crit_time(self, start_survey, tau):
        return start_survey - tau
    
    def latest_crit_time(self, end_survey, tau):
        return end_survey + tau
    
    def fluxint(self, F0, tcrit, tau, end_obs, start_obs):
        """Return the integrated flux"""
                  
        return (F0*tau*np.sqrt(np.pi/2.0)*(erf((end_obs-tcrit)/(tau*np.sqrt(2)))-erf((start_obs-tcrit)/(tau*np.sqrt(2)))))/(end_obs-start_obs)
        
    def gausscdf(self, x, t):
        return (((x*np.sqrt(2.0*np.pi)))*norm.cdf(t, loc = 0, scale = x)) 
            
    def lines(self, xs, ys, durmax, max_distance, flux_err, obs):
        gaps = np.array([],dtype=np.float32)
        for i in range(len(obs)-1):
            gaps = np.append(gaps, obs[i+1,0] - obs[i,0] + obs[i,1])
            # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
        min_sens = min(obs[:,2])
        max_sens = max(obs[:,2])
        sens_last = obs[-1,2]
        sens_maxgap = obs[np.where((gaps[:] == max(gaps)))[0]+1 ,2][0]
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        day1_obs = obs[0,1]
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        durmax_y = np.array([],dtype=np.float64)
        maxdist_y = np.array([],dtype=np.float64)
        for x in xs:
            # print(self.gausscdf(np.power(10,x),durmax + np.power(10,x)))
            try:
                # print(self.gausscdf(np.power(10,x),durmax + np.power(10,x)))
                durmax_y = np.append(durmax_y, ((1. + flux_err) * sens_last * day1_obs  ) / (self.gausscdf(np.power(10,x),durmax + np.power(10,x)) - self.gausscdf(np.power(10,x), durmax - day1_obs + np.power(10,x))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, ((1. + flux_err) * sens_last * day1_obs ) / (self.gausscdf(np.power(10,x),max_distance + day1_obs ) - self.gausscdf(np.power(10,x), max_distance)))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)
       

        maxdist_x = np.empty(len(ys))
        maxdist_x.fill(np.log10(max_distance))        
        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) & (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         