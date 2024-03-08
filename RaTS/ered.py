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
        

        # Q: Why is this part so complicated? 
        # A: Because of overflows, underflows, and all kinds of computational issues when tcrit falls outside of the observation. Also for debugging.
        # Q: Why is there a tau/2? 
        # A: Good question. Ok, so how do we want to define tau? We already decided in definitions of earliest_crit_time and latest_crit_time to 
        #    define tau as something that is equivalent to the tau in the fred and wilma cases. Therefore, because we have a wilma glued to a fred for this light curve,
        #    we set half of tau to be the exponential decay factor. 

        tau = tau/2
        
        # In order to avoid computational issues, we create masks for three different possible scenarios:
        # 1. The critical time of the transient is before the start of the observation and looks like a fred in the observation
        # 2. The critical time is after the end of the observation and looks like a wilma
        # 3. The critical time is during the observation and the wilma and fred portions must be calculated. 
        scen1 = np.zeros(len(tcrit),dtype=bool)
        scen2 = np.zeros(len(tcrit),dtype=bool)
        scen3 = np.zeros(len(tcrit),dtype=bool)
        scen1 += (tcrit < start_obs)
        scen2 += (tcrit > end_obs)
        scen3 += (tcrit <= end_obs) & (tcrit >= start_obs)


        flux_int = np.zeros(tcrit.shape[0],dtype=np.float64)
               
        # Scenario one 
        tstart1 = np.maximum(tcrit[scen1], start_obs) - tcrit[scen1]
        tend1 = end_obs - tcrit[scen1]
        flux_int[scen1] = F0[scen1]*tau[scen1]*(np.exp(-tstart1/tau[scen1]) - np.exp(-tend1/tau[scen1]))/(end_obs - start_obs)
        
        # Scenario two
        tstart2 = start_obs - tcrit[scen2]
        tend2 = np.minimum(end_obs,tcrit[scen2]) - tcrit[scen2] 
        flux_int[scen2] = F0[scen2]*tau[scen2]*(np.exp(tend2/tau[scen2]) - np.exp(tstart2/tau[scen2]))/(end_obs - start_obs)

        # Scenario three 
        tstart3f = np.maximum(tcrit[scen3], start_obs) - tcrit[scen3]
        tend3f = end_obs - tcrit[scen3]
        tstart3w = start_obs - tcrit[scen3]
        tend3w = np.minimum(end_obs,tcrit[scen3]) - tcrit[scen3] 
        flux_int[scen3] = F0[scen3]*tau[scen3]*(np.exp(-tstart3f/tau[scen3]) - np.exp(-tend3f/tau[scen3]) + np.exp(tend3w/tau[scen3]) - np.exp(tstart3w/tau[scen3]))/(end_obs - start_obs)

        return flux_int

        
    def lines(self, xs, ys, durmax, max_distance, flux_err, obs):
        gaps = np.zeros(len(obs)-1,dtype=np.float32)
        for i in range(len(obs)-1):
            gaps[i] =  obs['start'][i+1] - obs['start'][i] + obs['duration'][i]
            # gaps = np.append(gaps, obs[i+1,0] - obs[i,0])
        min_sens = min(obs['sens'])
        max_sens = max(obs['sens'])
        
        sens_maxgapbefore = obs['sens'][np.where((gaps[:] == max(gaps)))[0]-1][0]
        sens_maxgapafter = obs['sens'][np.where((gaps[:] == max(gaps)))[0]+1][0]
        sens_maxgap = min(sens_maxgapbefore, sens_maxgapafter)
        sens_argmin = np.argmin(np.array([sens_maxgapbefore, sens_maxgapafter]))
        gapobs = obs['duration'][sens_argmin]
        # sens_maxgap = obs['sens'][np.where((gaps[:] == max(gaps)))[0]+1][0]
        durmax_y = np.zeros(xs.shape,dtype=np.float64)
        maxdist_y = np.zeros(xs.shape,dtype=np.float64)
        day1_obs = obs['duration'][0]
        lastday_obs = obs['duration'][-1]
        sens_last = obs['sens'][-1]
        sens_first = obs['sens'][0]
        sens_maxtime = max(sens_last, sens_first)
        maxargobs = np.argmin(np.array([sens_last, sens_first]))
        maxtime_obs = obs['duration'][maxargobs]
        print(maxtime_obs)
        print(durmax)
        print(sens_maxtime)
        for i in range(len(xs)):
            x = xs[i]
            try:
                durmax_y[i] =  (1. + flux_err) * maxtime_obs * sens_maxtime / (np.power(10,x)/2) / (np.exp(-(durmax  / np.power(10,x))) - np.exp(-(durmax + 2*maxtime_obs) /  np.power(10,x)))
            except:
                durmax_y[i] = np.inf
            try:
                # maxdist_y[i] =  (((1. + flux_err) * sens_maxgap * day1_obs) /  np.power(10,x))   / (np.exp(-(max_distance / np.power(10,x))) - np.exp(-(max_distance + day1_obs) / np.power(10,x)))
                maxdist_y[i] =  (((1. + flux_err) * sens_maxgap * gapobs) /  (np.power(10,x)/2)   / (np.exp(-(max_distance / np.power(10,x))) - np.exp(-(max_distance + 2*gapobs) / np.power(10,x))))
            except:
                maxdist_y[i] = np.inf    
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
         