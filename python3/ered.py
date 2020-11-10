import numpy as np
import sys
from fred import fred
from wilma import wilma
class ered:
    """exponential rise exponential decay lightcurve class"""
    def __init__(self):
        self.edges=[0,0] # 1 is a definite edge, tophat is the default and has a definite beginning and end. Therefore it is [1,1]
        
    def earliest_crit_time(self, start_survey, tau):
        return start_survey - tau
    
    def latest_crit_time(self, end_survey, tau):
        return end_survey + tau
    
    def fluxint(self, F0, tcrit, tau, end_obs, start_obs):
        """Return the integrated flux"""
        fredcls = fred()
        wilmacls = wilma()
        scen1 = np.where(tcrit < start_obs)[0]
        scen2 = np.where(tcrit > end_obs)[0]
        scen3 = np.where((tcrit <= end_obs) & (tcrit >= start_obs))[0]
        fluxint = np.array([],dtype=np.float64)
        fluxint1 = fredcls.fluxint(F0[scen1], tcrit[scen1], tau[scen1]/2.0, end_obs, start_obs)
        fluxint2 = wilmacls.fluxint(F0[scen2], tcrit[scen2], tau[scen2]/2.0, end_obs, start_obs)
        fluxint3 = wilmacls.fluxint(F0[scen3], tcrit[scen3], tau[scen3]/2.0, end_obs, start_obs) + fredcls.fluxint(F0[scen3], tcrit[scen3], tau[scen3]/2.0, end_obs, start_obs)
        fluxint = np.zeros(len(fluxint1) + len(fluxint2) + len(fluxint3))
        fluxint[0:len(fluxint1)] = fluxint[0:len(fluxint1)] + fluxint1
        fluxint[len(fluxint1):(len(fluxint2)+len(fluxint1))] = fluxint[len(fluxint1):(len(fluxint2)+len(fluxint1))] + fluxint2
        fluxint[(len(fluxint2)+len(fluxint1)):(len(fluxint2)+len(fluxint1)+len(fluxint3))] = fluxint[(len(fluxint2)+len(fluxint1)):(len(fluxint2)+len(fluxint1)+len(fluxint3))] + fluxint3
        return fluxint
        
        # argofexponential1 = (start_obs - tcrit)/(tau/2.0)
        # exponential1 = np.exp(argofexponential1)
        # argofexponential2 = -(end_obs-tcrit)/(tau/2.0)
        # try: 
            # exponential2 = np.exp(argofexponential2)
        # except:
            # print(np.amax(-(end_obs-tcrit)))
            # print(end_obs)
            # print(tcrit)
            # print(start_obs)
            # print(np.amax(argofexponential1))
            # sys.exit(1)
        # return ((F0*(tau/2.0))/(end_obs-start_obs))*(2.0-exponential1-exponential2)
        
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
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs / np.power(10,x) / (2.0 -  np.exp(-(durmax - day1_obs + np.power(10,x)) /  np.power(10,x)) - np.exp(-((durmax + np.power(10,x)) / np.power(10,x)))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, (((1. + flux_err) * sens_maxgap * day1_obs) /  np.power(10,x))   / (2.0 - np.exp(-(max_distance / np.power(10,x))) - np.exp(-(max_distance + day1_obs) / np.power(10,x))))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)
        durmax_x = np.empty(len(ys))
        durmax_x.fill(np.log10(durmax))
        maxdist_x = np.empty(len(ys))
        maxdist_x.fill(np.log10(max_distance))        
        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) & (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         