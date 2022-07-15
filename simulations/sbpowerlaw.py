import numpy as np
import sys
from scipy.special import hyp2f1
from fractions import Fraction
from decimal import Decimal
class sbpowerlaw:
    """smoothly broken power law lightcurve class"""
    # class variables
    edges=[0,0] # 1 is a definite edge, tophat is the default and has a definite beginning and end. Therefore it is [1,1]
    fractionalcut = 100
    def __init__(self, alpha1=-0.8, alpha2=2.1, Delta=0.1):
        # sample indices. alpha1 must be negative
        # alpha2 must be positive. 0 < Delta < 1 Convention differs
        # from Mooley et al. 2018 https://arxiv.org/abs/1810.12927
        self.alpha1 = alpha1
        self.alpha2 = alpha2
        self.Delta = Delta # smoothness AKA (1/s) 

    def earliest_crit_time(self, start_survey, tau):       
        return start_survey 
    
    def latest_crit_time(self, end_survey, tau):
        # This function is ultimately not important
        return end_survey 
    
    def indefiniteint(self, x, A, xb):
        alpha1 = self.alpha1
        alpha2 = self.alpha2
        Delta = self.Delta

        try:
            p1 = 0.5**((alpha1-alpha2)/Delta)
            if bool(xb.shape):
                p2 = np.zeros(xb.shape)
                for i in range(len(xb)):
                    p2pre = Fraction.from_float(-x/xb[i])**(Fraction(-alpha1))
                    if type(p2pre) is complex:
                        p2[i] = -np.absolute(p2pre)
                    else:
                        p2[i] = p2pre
                p3 = -(x/xb)**(1/Delta)
                retval = np.zeros(p2.shape)
                for i in range(len(p2)):
                    numerator = -(Decimal(A[i])*Decimal(x)*Decimal(p1)*Decimal(p2[i])*Decimal(np.nan_to_num(hyp2f1((alpha2-alpha1)/Delta, Delta-alpha1*Delta, -alpha1*Delta+Delta+1,p3[i] ))))      
                    denom =Decimal(alpha1-1) 
                    retval[i] = np.nan_to_num(np.float64(numerator) / np.float64(denom) )
            else:
                p2pre = Fraction.from_float(-x/xb)**(Fraction(-alpha1))
                if type(p2pre) is complex:
                    p2 = -np.absolute(p2pre)
                else:
                    p2 = p2pre
                p3 = -(x/xb)**(1/Delta)
                numerator = -(Decimal(A)*Decimal(x)*Decimal(p1)*Decimal(p2)*Decimal(np.nan_to_num(hyp2f1((alpha2-alpha1)/Delta, Delta-alpha1*Delta, -alpha1*Delta+Delta+1,p3 ))))      
                denom =Decimal(alpha1-1) 
                retval = np.nan_to_num(np.float64(numerator) / np.float64(denom) )

 
            return retval
        except Exception as e:
            print(e)
            print(numerator)
            print(denom)
    
    def fluxint(self, F0, tcrit, tau, end_obs, start_obs):
        """Return the integrated flux"""
        # break time is determined by the duration. duration
        # is set to be the duration for which the flux is above
        # 10% of the peak flux as calculated by the individual power
        # law components

        alpha1= self.alpha1
        alpha2 = self.alpha2
        Delta = self.Delta
        xb = tau / ((1/sbpowerlaw.fractionalcut)**(1/-alpha1) - (1/sbpowerlaw.fractionalcut)**(1/alpha2))
        A = F0
        intflux =  self.indefiniteint(end_obs,A, xb)- self.indefiniteint(start_obs, A, xb) 
        return intflux
            
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
                alpha1= self.alpha1
                alpha2 = self.alpha2
                Delta = self.Delta
                xb = x / ((1/sbpowerlaw.fractionalcut)**(1/-alpha1) - (1/sbpowerlaw.fractionalcut)**(1/alpha2))
                durmax_y[i] = (1.+flux_err)*sens_maxtime*(maxtime_obs)/self.indefiniteint(maxtime_obs,1, xb)     
            except RuntimeWarning:
                durmax_y[i]=np.inf
            try:
                alpha1= self.alpha1
                alpha2 = self.alpha2
                Delta = self.Delta
                xb = x / ((1/sbpowerlaw.fractionalcut)**(1/-alpha1) - (1/sbpowerlaw.fractionalcut)**(1/alpha2))
                maxdist_y[i] = (1.+flux_err)*sens_maxgap*(sens_maxdur)/self.indefiniteint(sens_maxdur,1, xb)
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
         
