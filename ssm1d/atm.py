import numpy as np
from scipy.interpolate import interp1d
def get_custom_atm(par,vres1=np.arange(0,3e4,1e2),vres2=np.arange(0,3e4,1e2)):
    # scalar inputs: mid-tropospheric RH (RHmid) and temeprature (Tmid) and boolean option for uniform RH
    # Vertical resolution can also be specified to facilitate RFM vs analytical model. 
    #############################################################
    def get_optimal_RH(Tsfc2trp, Ts=300, Tmid=260, Ttrp=200, RHs=0.75, RHmid=0.4, RHtrp=0.75):
        print('Ts',Ts)
        print('Tmid',Tmid)
        print('Ttrp',Ttrp)
        print('RHs',RHs)
        print('RHmid',RHmid)
        print('RHtrp',RHtrp)
        # Compute alpha_left and alpha_right based on target RH values
        alpha_left  = -np.log(RHtrp / RHmid)/(Ttrp - Tmid) ** 2
        alpha_right = -np.log(RHs / RHmid)/(Ts - Tmid) ** 2 
        # Compute RH values for the provided temperature range
        RH_opt = [
            RHmid * np.exp(-(temp - Tmid) ** 2 * alpha_left) if temp < Tmid else RHmid * np.exp(-(temp - Tmid) ** 2 * alpha_right)
            for temp in Tsfc2trp
        ]
        return {'RH':RH_opt,'alpha_lt':alpha_left, 'alpha_gt':alpha_right}
    #############################################################
    def get_esat_over_l(par,T):
        import math
        # SVP over liquid from DAM (Pa)
        xi = 1 # factor to change the SVP
        pstar = xi*par.ptrip * (T/par.Ttrip)**((par.cpv-par.cvl)/par.rgasv) * math.exp( (par.E0v - (par.cvv-par.cvl)*par.Ttrip) / par.rgasv * (1/par.Ttrip - 1/T) )
        return pstar
    #############################################################
    # Lapse Rate and Temperature (Troposphere)
    Gamma = np.ones([len(par.z)],dtype='float')*par.Gamma
    T  = par.Ts - par.Gamma*par.z
    # stratospheric mask
    mask = np.where(T<par.Ttrp)
    # Lapse Rate and Temperature (Stratosphere)
    T  = np.where(T<par.Ttrp,par.Ttrp,T)
    Gamma[mask] = 0
    # Identify the height of the tropopause
    ktrp = int(np.amin(np.where(T==par.Ttrp)))
    ztrp = par.z[ktrp]
    # Custom RH
    if par.uniform is True:
        RH       = np.ones([len(par.z)])*par.RHmid
        RH[mask] = 0
        alpha_lt = 0
        alpha_gt = 0
    else:
        RH  = np.ones([len(par.z)])*par.RHs
        RH[mask] = 0 # stratospheric mask
        foo = get_optimal_RH(T[0:(ktrp+1)], par.Ts, par.Tmid, par.Ttrp, par.RHs, par.RHmid, par.RHtrp)
        RH[0:(ktrp+1)] = foo['RH']
        alpha_lt       = foo['alpha_lt']
        alpha_gt       = foo['alpha_gt']
    
    # Solve for environmental pressure and density 
    rho  = np.zeros_like(par.z)
    p    = np.zeros_like(par.z)
    p[0] = par.ps + RH[0]*get_esat_over_l(par,par.Ts) # 10^5(dry)+pH2O
    arg  = 0
    # initialize molar fraction (relative to total)
    xN2  = np.zeros_like(par.z)
    xO2  = np.zeros_like(par.z)
    xH2O = np.zeros_like(par.z)
    xCO2 = np.zeros_like(par.z)
    # initialize mass fraction (relative to total)
    wN2  = np.zeros_like(par.z)
    wO2  = np.zeros_like(par.z)
    wH2O = np.zeros_like(par.z)
    wCO2 = np.zeros_like(par.z)
    # initialize specific gas constant
    Rtot  = np.zeros_like(par.z)
    # initialize mean molecular mass
    Mave  = np.zeros_like(par.z)
    # initialize specific heat of air
    cp    = np.zeros_like(par.z) # J/kg/K
    cpmol = np.zeros_like(par.z) # J/mol/K
    for k in range(len(par.z)):
        if k<(len(par.z)-1):
            dz   = par.z[k+1]-par.z[k]
        else:
            dz   = par.z[1]-par.z[0]
        # molar mixing ratio of H2O
        if k<=ktrp: # tropospheric value set by Clausius-Clapeyron
            pH2O    = RH[k]*get_esat_over_l(par,T[k])
            xH2O[k] = pH2O/p[k]
        else: # stratospheric mixing ratio fixed to tropopause value
            xH2O[k] = xH2O[ktrp]
            # calculate the "implied" relative humidity
            RH[k]   = xH2O[k]*p[k]/get_esat_over_l(par,T[k])
        # Compute dry-air scaling factor (applies to N2, O2, and CO2)
        scaling_factor = 1 - xH2O[k]  # The fraction of dry air remaining
        # Apply fixed dry-air composition (keeping N2, O2, and CO2 proportions constant in dry-air)
        xN2[k]  = scaling_factor * par.xN2             # 78%
        xCO2[k] = scaling_factor * par.xCO2            # 395 ppm
        xO2[k]  = scaling_factor * par.xO2 - xCO2[k]   # ~22% (note: robbing from over-estimated O2)
        tol = 1e-6 # tolerance for precision in mole fractions
        if abs(par.xN2+par.xO2-1)>tol:
            print('Error: sum of dry-air abundances of N2,O2 is non-unity. Tolerance exceeded')
        if abs(xN2[k] + xO2[k] + xH2O[k] + xCO2[k] - 1) > tol:
            print(f'Error: sum of molar fractions is non-unity. Tolerance exceeded.')
            #print(abs(xN2[k] + xO2[k] + xH2O[k] + xCO2[k] - 1))
        # mean molecular weight of air
        Mave[k]  = xN2[k]*par.MN2 + xO2[k]*par.MO2 + xH2O[k]*par.MH2O + xCO2[k]*par.MCO2
        # molar specific heat of air (J/mol/K)
        cpmol[k] = xN2[k]*par.cpmolN2 + xO2[k]*par.cpmolO2 + xH2O[k]*par.cpmolH2O + xCO2[k]*par.cpmolCO2
        # mass mixing ratios of N2, O2, H2O, CO2
        wN2[k]  = xN2[k] *par.MN2/Mave[k]
        wO2[k]  = xO2[k] *par.MO2/Mave[k]
        wH2O[k] = xH2O[k]*par.MH2O/Mave[k]
        wCO2[k] = xCO2[k]*par.MCO2/Mave[k]
        if abs(wN2[k] + wO2[k] + wH2O[k] + wCO2[k] - 1) > tol:
            print(f'Error: sum of mixing ratios is non-unity. Tolerance exceeded.')
        # specific heat of air (J/kg/K)
        cp[k] = wN2[k]*par.cpN2 + wO2[k]*par.cpO2 + wH2O[k]*par.cpH2O + wCO2[k]*par.cpCO2
        # specific gas constant of air
        Rtot[k]  = wN2[k]*par.RN2 + wO2[k]*par.RO2 + wH2O[k]*par.RH2O + wCO2[k]*par.RCO2
        # solve for total density of air
        rho[k] = p[k]/(Rtot[k]*T[k])
        # solve for exponential term
        arg    = -par.ggr/Rtot[k]*dz/T[k]
        # solve for pressure at next level
        if k<(len(par.z)-1):
            p[k+1] = p[k]*np.exp(arg)
    #############################################################
    # Export fields in their desired vertical resolution (vres1,vres2)
    def interpolate(var, vres):
        # func = interp1d(par.z,var,kind)(input of function)
        return interp1d(par.z, var, kind='cubic')(vres)
    T1, T2         = interpolate(T, vres1), interpolate(T, vres2)
    Gamma1, Gamma2 = interpolate(Gamma, vres1), interpolate(Gamma, vres2)
    p1, p2         = interpolate(p, vres1), interpolate(p, vres2)
    rho1, rho2     = interpolate(rho, vres1), interpolate(rho, vres2)
    RH1, RH2       = interpolate(RH, vres1), interpolate(RH, vres2)
    xN2_1, xN2_2   = interpolate(xN2, vres1), interpolate(xN2, vres2)
    xO2_1, xO2_2   = interpolate(xO2, vres1), interpolate(xO2, vres2)
    xH2O_1, xH2O_2 = interpolate(xH2O, vres1), interpolate(xH2O, vres2)
    xCO2_1, xCO2_2 = interpolate(xCO2, vres1), interpolate(xCO2, vres2)
    cp1, cp2       = interpolate(cp, vres1), interpolate(cp, vres2)
    cpmol1, cpmol2 = interpolate(cpmol, vres1), interpolate(cpmol, vres2)
    #############################################################
    dat1 = {'T': T1, 'p': p1, 'Gamma': Gamma1, 'rho': rho1, 'z': vres1, 'RH': RH1,
            'xN2': xN2_1, 'xO2': xO2_1, 'xH2O': xH2O_1, 'xCO2': xCO2_1, 'cp':cp1, 'cpmol':cpmol1,
           'alpha_lt':alpha_lt,'alpha_gt':alpha_gt}
    dat2 = {'T': T2, 'p': p2, 'Gamma': Gamma2, 'rho': rho2, 'z': vres2, 'RH': RH2,
            'xN2': xN2_2, 'xO2': xO2_2, 'xH2O': xH2O_2, 'xCO2': xCO2_2, 'cp':cp2, 'cpmol':cpmol2,
           'alpha_lt':alpha_lt,'alpha_gt':alpha_gt}
    #############################################################
    # T(K), RH(unitless), p(Pa), Gamma(K/m), rho(kg/m3), x(molar mixing ratio)
    #############################################################
    return dat1,dat2
