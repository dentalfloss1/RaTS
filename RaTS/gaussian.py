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
        
        return (F0*tau*np.sqrt(np.pi/8.0)*(erf(np.sqrt(2)*(end_obs-tcrit)/(tau))-erf(np.sqrt(2)*(start_obs-tcrit)/(tau))))/(end_obs-start_obs)
        
    def gausscdf(self, x, t):
        return x*np.sqrt(np.pi/2)*erf(t/x/np.sqrt(2))
            
    def lines(self, xs, ys, durmax, max_distance, flux_err, obs):
        gaps = np.array([],dtype=np.float32)
        for i in range(len(obs)-1):
            gaps = np.append(gaps, obs['start'][i+1] - obs['start'][i] + obs['duration'][i])
            # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
        min_sens = min(obs['sens'])
        max_sens = max(obs['sens'])
        sens_last = obs['sens'][-1]
        day1_obs = obs['duration'][0]
        # durmax_x = np.empty(len(ys))
        # durmax_x.fill(np.log10(durmax))
        durmax_y = np.zeros(xs.shape,dtype=np.float64)
        maxdist_y = np.zeros(xs.shape,dtype=np.float64)
        sens_maxgapbefore = obs['sens'][np.where((gaps[:] == max(gaps)))[0]-1][0]
        sens_maxgapafter = obs['sens'][np.where((gaps[:] == max(gaps)))[0]+1][0]
        sens_maxgap = min(sens_maxgapbefore, sens_maxgapafter)
        sens_argmin = np.argmin(np.array([sens_maxgapbefore, sens_maxgapafter]))
        sens_maxdur = obs['duration'][sens_argmin]
        sens_last = obs['sens'][-1]
        sens_first = obs['sens'][0]
        sens_maxtime = max(sens_last, sens_first)
        maxargobs = np.argmin(np.array([sens_last, sens_first]))
        maxtime_obs = obs['duration'][maxargobs]
        for i in range(len(xs)):
            x = np.power(10,xs[i])
            try:
                durmax_y[i] = (1.+flux_err)*sens_maxtime*(maxtime_obs)/x/np.sqrt(np.pi/8.0)/(erf(np.sqrt(2)*(-durmax/2 + maxtime_obs)/(x))-erf(np.sqrt(2)*(-durmax/2)/(x)))
            except RuntimeWarning:
                durmax_y[i]=np.inf
            try:
                maxdist_y[i] = (1.+flux_err)*sens_maxgap*(sens_maxdur)/x/np.sqrt(np.pi/8.0)/(erf(np.sqrt(2)*(-max_distance/2)/(x))-erf(np.sqrt(2)*(-(sens_maxdur + max_distance/2))/(x)))
            except RuntimeWarning:
                maxdist_y[i] = np.inf
        #     # print(self.gausscdf(np.power(10,x),durmax + np.power(10,x)))

        #     try:
        #         # print(self.gausscdf(np.power(10,x),durmax + np.power(10,x)))
        #         durmax_y = np.append(durmax_y, ((1. + flux_err) * sens_last * day1_obs  ) / (self.gausscdf(np.power(10,x),durmax + np.power(10,x)) - self.gausscdf(np.power(10,x), durmax - day1_obs + np.power(10,x))))
        #     except:
        #         durmax_y = np.append(durmax_y, np.inf)
        #     try:
        #         maxdist_y =  np.append(maxdist_y, ((1. + flux_err) * sens_maxgap * day1_obs ) / (self.gausscdf(np.power(10,x),max_distance + day1_obs ) - self.gausscdf(np.power(10,x), max_distance)))
        #     except:
        #         maxdist_y = np.append(maxdist_y, np.inf)
        durmax_x = ' '
        maxdist_x = ' '
        # maxdist_x = np.empty(len(ys))
        # maxdist_x.fill(np.log10(max_distance))        
        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) & (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         