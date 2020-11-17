import numpy as np
import sys
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
        scen1 = np.where(tcrit < start_obs)[0]
        scen2 = np.where(tcrit > end_obs)[0]
        scen3 = np.where((tcrit <= end_obs) & (tcrit >= start_obs))[0]
        flux_int = np.zeros(tcrit.shape[0],dtype=np.float64)
               
        
        tstart1 = np.maximum(tcrit[scen1], start_obs) - tcrit[scen1]
        tend1 = end_obs - tcrit[scen1]
        flux_int[scen1] = np.multiply(F0[scen1], np.multiply(tau[scen1], np.divide(np.exp(-np.divide(tstart1,tau[scen1])) - np.exp(-np.divide(tend1,tau[scen1])), (end_obs-start_obs))))
        
        tstart2 = start_obs - tcrit[scen2]# Burst really starts, so we always start at the beginning of the observation
        tend2 = np.minimum(end_obs,tcrit[scen2]) - tcrit[scen2] 
        flux_int[scen2] = np.multiply(F0[scen2], np.multiply(tau[scen2],  np.divide(np.exp(np.divide(tend2,tau[scen2])) - np.exp(np.divide(tstart2,tau[scen2])) , (end_obs-start_obs))))

        tstart3f = np.maximum(tcrit[scen3], start_obs) - tcrit[scen3]
        tend3f = end_obs - tcrit[scen3]
        tstart3w = start_obs - tcrit[scen3]# Burst really starts, so we always start at the beginning of the observation
        tend3w = np.minimum(end_obs,tcrit[scen3]) - tcrit[scen3] 
        flux_int[scen3] = np.multiply(F0[scen3], np.multiply(tau[scen3], np.divide(np.exp(-np.divide(tstart3f,tau[scen3])) - np.exp(-np.divide(tend3f,tau[scen3])), (end_obs-start_obs)))) + np.multiply(F0[scen3], np.multiply(tau[scen3],  np.divide(np.exp(np.divide(tend3w,tau[scen3])) - np.exp(np.divide(tstart3w,tau[scen3])) , (end_obs-start_obs))))
        # fluxint = np.zeros(len(fluxint1) + len(fluxint2) + len(fluxint3))
        # fluxint[0:len(fluxint1)] = fluxint[0:len(fluxint1)] + fluxint1
        # fluxint[len(fluxint1):(len(fluxint2)+len(fluxint1))] = fluxint[len(fluxint1):(len(fluxint2)+len(fluxint1))] + fluxint2
        # fluxint[(len(fluxint2)+len(fluxint1)):(len(fluxint2)+len(fluxint1)+len(fluxint3))] = fluxint[(len(fluxint2)+len(fluxint1)):(len(fluxint2)+len(fluxint1)+len(fluxint3))] + fluxint3
       # fluxint = np.append(0*fluxint1,fluxint2)

        # print(fluxint.shape, fluxint1.shape, fluxint2.shape, fluxint3.shape, fluxint1.shape[0] + fluxint2.shape[0] + fluxint3.shape[0])
        return flux_int
        
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
        for x in xs:
            try:
                durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs / np.power(10,x) / (np.exp(-(durmax - day1_obs + np.power(10,x)) /  np.power(10,x)) - np.exp(-((durmax + np.power(10,x)) / np.power(10,x)))))
            except:
                durmax_y = np.append(durmax_y, np.inf)
            try:
                maxdist_y =  np.append(maxdist_y, (((1. + flux_err) * sens_maxgap * day1_obs) /  np.power(10,x))   / (np.exp(-(max_distance / np.power(10,x))) - np.exp(-(max_distance + day1_obs) / np.power(10,x))))
            except:
                maxdist_y = np.append(maxdist_y, np.inf)    
        durmax_x = ' '
        maxdist_x = ' '
        durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) &  (durmax_y > np.amin(10**ys)))[0]
        maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]


       # gaps = np.array([],dtype=np.float32)
        # for i in range(len(obs)-1):
            # gaps = np.append(gaps, obs[i+1,0] - obs[i,0] + obs[i,1])
            # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
        # min_sens = min(obs[:,2])
        # max_sens = max(obs[:,2])
        # sens_last = obs[-1,2]
        # sens_maxgap = obs[np.where((gaps[:] == max(gaps)))[0]+1 ,2][0]
        # durmax_y = np.array([],dtype=np.float64)
        # maxdist_y = np.array([],dtype=np.float64)
        # day1_obs = obs[0,1]
        # durmax_y = np.array([],dtype=np.float64)
        # maxdist_y = np.array([],dtype=np.float64)
        # for x in xs:
            # try:
                # durmax_y = np.append(durmax_y, (1. + flux_err) * sens_last * day1_obs / np.power(10,x) / (-np.exp((durmax - day1_obs + np.power(10,x)) /  np.power(10,x)) + np.exp(-((durmax + np.power(10,x)) / np.power(10,x)))))
            # except:
                # durmax_y = np.append(durmax_y, np.inf)
            # try:
                # maxdist_y =  np.append(maxdist_y, (((1. + flux_err) * sens_maxgap * day1_obs) /  np.power(10,x))   / (-np.exp((max_distance / np.power(10,x))) + np.exp(-(max_distance + day1_obs) / np.power(10,x))))
            # except:
                # maxdist_y = np.append(maxdist_y, np.inf)
        # durmax_x = np.empty(len(ys))
        # durmax_x.fill(np.log10(durmax))
        # maxdist_x = np.empty(len(ys))
        # maxdist_x.fill(np.log10(max_distance))        
        # durmax_y_indices = np.where((durmax_y < np.amax(10**ys)) & (durmax_y > np.amin(10**ys)))[0]
        # maxdist_y_indices = np.where((maxdist_y < np.amax(10**ys)) & (maxdist_y > np.amin(10**ys)))[0]
        #fredcls = fred()
        #durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices = fredcls.lines(xs, ys, durmax, max_distance, flux_err, obs)   
        return  durmax_x, maxdist_x, durmax_y, maxdist_y, durmax_y_indices, maxdist_y_indices
         