# Load dependencies
import numpy as np
import matplotlib.pyplot as plt
from collections.abc import Sequence
from typing import Type,Dict,Tuple
import scipy.optimize as opt
            
# Helper functions

# thermodynamic corrective factors arising from RH(T)
def get_phi(par,T):
    alpha = np.where(T<par.Tmid,par.alpha_lt,par.alpha_gt)
    return -2*alpha*T*(T-par.Tmid)
def get_phipr(par,T):
    alpha = np.where(T<par.Tmid,par.alpha_lt,par.alpha_gt)
    return -2*alpha*T*(2*T-par.Tmid)

# Planck function
def get_B(nu:Sequence[float],T:int)->Dict[str,float]:
    # Planck emission ((W/m2)*cm) as a function of wavenumber (nu; array; cm-1) and temperature (T; scalar; K)
    nu = nu*100 # 1/cm -> 1/m
    h = 6.626e-34 # J*s
    c = 3.0e8     # m/s
    k = 1.38e-23  # J/K
    B = 2*h*c**2*nu**3/(np.exp((h*c*nu)/(k*T))-1) # kg*m/s^2 = W m^(-2) m^(+1)
    B = 100*B # W m^(-2) m^(+1) -> W*cm/m^2
    return {'B':B}

# Water vapor path (WVP; kg/m2)
def get_WVP(par:Dict[str,float],RH:float,T:float)->float:
    WVP0 = par.pinf/(par.Rv*par.Gamma) # reference water vapor path
    phi  = get_phi(par,T) # corrective RH factor
    return WVP0*RH*np.exp(-par.L/(par.Rv*T))/(par.L/(par.Rv*T)+phi)

# spectral functions for rotation band
def get_xirot_lo(nu,par):
    return par.Arot_lo*np.sin(2*np.pi*(nu-par.nurot)/par.Prot_lo)
def get_xirot_co(nu,par):
    return par.Arot_co*np.cos(2*np.pi*(nu-par.nurot)/par.Prot_co)
def get_xirot(nu,par):
    expr1 = get_xirot_lo(nu,par)
    expr2 = get_xirot_co(nu,par)
    return expr1+expr2
def get_dxirot_lo(nu,par):
    return par.Arot_lo*np.cos(2*np.pi*(nu-par.nurot)/par.Prot_lo)*(2*np.pi/par.Prot_lo)
def get_dxirot_co(nu,par):
    return -par.Arot_co*np.sin(2*np.pi*(nu-par.nurot)/par.Prot_co)*(2*np.pi/par.Prot_co)
def get_dxirot(nu,par):
    expr1 = get_dxirot_lo(nu,par)
    expr2 = get_dxirot_co(nu,par)
    return expr1+expr2
def get_lrot(nu,par):
    xirot  = get_xirot(nu,par)
    dxirot = get_dxirot(nu,par)
    lrot = 1/np.absolute(-par.lrot0**(-1) + dxirot/(1+xirot))
    return lrot

# spectral functions for vibration-rotation band
def get_xivib_lo(nu,par):
    return 0
def get_xivib_co(nu,par):
    return 0
def get_xivib(nu,par):
    return 0
def get_lvib(nu,par):
    return 0

# spectral representation of kref
def get_kref_rot(nu,par):
    return par.krot*np.exp(-(nu-par.nurot)/par.lrot0)*(1+get_xirot(nu,par))
def get_kref_vib(nu,par):
    return par.kvib*np.exp(-(par.nuvib-nu)/par.lvib0)*(1+get_xivib(nu,par))

# root solver for emitting wavenumber in rotation band
def nu1rot_from_rootsolve_new(target_blah:float,par:Dict[str,float]):
    import scipy.optimize as opt
    # solves for emitting wavenumber given [blah]=-ln(kj)
    # [blah] refers to terms lacking spectral dependence
    #######################################################
    def obj_func(x):
        # solve for nu where target_blah=-ln(kj)
        return (target_blah+np.log(get_kref_rot(x,par)))
    #######################################################
    # Return the roots of the (non-linear) equations defined by obj_func(x) = 0 given a starting estimate
    # There could be multiple roots, so we do a grid search with 25 cm-1 spacing
    nu1_solve = []
    nu_guesses = np.arange(par.nus[par.i0],par.nus[par.i1]+50,50)
    for i,nu_guess in enumerate(nu_guesses):
        tmp = opt.fsolve(obj_func,nu_guess,xtol=1e-8)[0]
        tmp = np.round(tmp,1) # round to nearest 0.1 cm-1
        # fsolve sometimes returns "garbage" roots: verify then add to list
        if abs(obj_func(tmp))<0.1: # threshold to be a root (conservative)
            if tmp not in nu1_solve:
                nu1_solve.append(tmp)
    if nu1_solve: # roots exist
        nu1_solve = np.array(nu1_solve)
        #print(nu1_solve)
        # keep the root closest to the nu1 in straight-line approximation
        nu1_sla = par.nurot + par.lrot0*(target_blah+np.log(par.krot))
        dnu1    = np.absolute(nu1_solve-nu1_sla)
        find    = np.where(dnu1==np.amin(dnu1))
        knu1    = int(find[0])
        # keep the closest root
        return nu1_solve[knu1]
    else: # no roots exist
        #print(f'error: no roots to be found')
        return np.nan

# root solver for emitting wavenumber in vibration-rotation band
def nu1vib_from_rootsolve_new(target_blah:float,par:Dict[str,float]):
    import scipy.optimize as opt
    # solves for emitting wavenumber given [blah]=-ln(kj)
    # [blah] refers to terms lacking spectral dependence
    #######################################################
    def obj_func(x):
        # solve for nu where target_blah=-ln(kj)
        return (target_blah+np.log(get_kref_vib(x,par)))
    #######################################################
    # Return the roots of the (non-linear) equations defined by obj_func(x) = 0 given a starting estimate
    # There could be multiple roots, so we do a grid search with 25 cm-1 spacing
    nu1_solve = []
    nu_guesses = np.arange(par.nus[par.i0],par.nus[par.i1]+50,50)
    for i,nu_guess in enumerate(nu_guesses):
        tmp = opt.fsolve(obj_func,nu_guess,xtol=1e-8)[0]
        tmp = np.round(tmp,1) # round to nearest 0.1 cm-1
        # fsolve sometimes returns "garbage" roots: verify then add to list
        if abs(obj_func(tmp))<0.1: # threshold to be a root (conservative)
            if tmp not in nu1_solve:
                nu1_solve.append(tmp)
    if nu1_solve: # roots exist
        nu1_solve = np.array(nu1_solve)
        #print(nu1_solve)
        # keep the root closest to the nu1 in straight-line approximation
        nu1_sla = par.nuvib - par.lvib0*(target_blah+np.log(par.kvib))
        dnu1    = np.absolute(nu1_solve-nu1_sla)
        find    = np.where(dnu1==np.amin(dnu1))
        knu1    = int(find[0])
        # keep the closest root
        return nu1_solve[knu1]
    else: # no roots exist
        #print(f'error: no roots to be found')
        return np.nan

def get_SSM1D(par,dataset: Dict[str, float]) -> Dict[str, float]:
    #############################################################
    # The parameters class defines the following attributes:
    #############################################################
    # kj: absorption coefficient at band-maximum of best-fit line (m2/kg)
    # lj0: exponential falloff coefficient of best-fit line (cm-1)
    # nuj: wavenumber at band-maximum (~150 or ~1500 cm-1)
    # Pj_lo: period of error oscillation (cm-1) of line-only component
    # Aj_lo: amplitude of error oscillation of line-only component
    # Pj_co: period of error oscillation (cm-1) of continuum-only component
    # Aj_co: amplitude of error oscillation of continuum-only component
    #############################################################
    if par.band == 'wv-rot-right':
        return _compute_ssm1d_rot_right(par,dataset)
    elif par.band == 'wv-vib-rot':
        return _compute_ssm1d_vib_rot(par,dataset)
    else:
        raise ValueError(f"Unrecognized band: {par.band}")

def _compute_ssm1d_rot_right(par,dataset):
    #############################################################
    # Unpack the inputs
    cp, z, T, p, rho, RH, Gamma = (dataset['RFM'][key] for key in ['cp', 'z', 'T', 'p', 'rho', 'RH', 'Gamma']) 
    #############################################################
    # use a root solver to find the emitting wavenumber
    # ln(D*p/pref*WVP) = -ln(kj)
    nu1rot = [] # cm-1 (nlev)
    for k in range(len(z)):
        target_blah = np.log(par.D*p[k]/par.pref*get_WVP(par,RH[k],T[k]))
        val = nu1rot_from_rootsolve_new(target_blah,par)
        nu1rot.append(val if val>=150 and val<=1000 else np.nan)
    nu1rot = np.array(nu1rot)

    fig, ax = plt.subplots()
    ax.plot(nu1rot,z/1e3)
    #############################################################
    # calculate blackbody emission as a function of (T,nu1rot)
    piB1rot = [] # (W*cm/m^2) (nlev)
    for k in range(len(z)):
        if not np.isnan(nu1rot[k]):
            val = get_B(nu1rot[k],T[k])
            piB1rot.append(np.pi*val['B'])
        else:
            piB1rot.append(np.nan)
    piB1rot = np.array(piB1rot)
    #############################################################
    # At each height, calculate Beta/z (nlev; do not take tropospheric mean)
    c1    = Gamma*par.Rd/par.ggr
    c2    = par.L/(par.Rv*T)
    phi   = get_phi(par,T)
    phipr = get_phipr(par,T)
    Boverz = (1 + c1*(c2+phi) + c1/(c2+phi)*(c2-phipr))*par.ggr/(par.Rd*T)
    #############################################################
    # Compute the emitting width at the emitting wavenumber
    lrot1 = []
    for k in range(len(z)):
        # wavenumber doing the emitting at this height
        nu1 = nu1rot[k]
        if not np.isnan(nu1):
            # use functional approximation of lrot
            ltmp = get_lrot(nu1,par)
            lrot1.append(ltmp)
        else:
            lrot1.append(np.nan)
    lrot1 = np.array(lrot1)
    #############################################################
    # Compute the SSM1D with varying emitting width
    Hrot = (1/(cp*rho))*(Boverz)*(piB1rot)*(lrot1) # (W/m3->K/s)
    #############################################################
    return {'nu1rot':nu1rot,'z1rot':z,'piB1rot':piB1rot,'Hrot':Hrot,'Boverz':Boverz,'lrot':lrot1}

def _compute_ssm1d_vib_rot(par,dataset):
    #############################################################
    # Unpack the inputs
    cp, z, T, p, rho, RH, Gamma = (dataset['RFM'][key] for key in ['cp', 'z', 'T', 'p', 'rho', 'RH', 'Gamma'])
    #############################################################
    # use a root solver to find the emitting wavenumber
    # ln(D*p/pref*WVP) = -ln(kj)
    nu1vib = [] # cm-1 (nlev)
    for k in range(len(z)):
        target_blah = np.log(par.D*p[k]/par.pref*get_WVP(par,RH[k],T[k]))
        val = nu1vib_from_rootsolve_new(target_blah,par)
        nu1vib.append(val if val>1000 and val<=1500 else np.nan)
    nu1vib = np.array(nu1vib)
    #############################################################
    # calculate blackbody emission as a function of (T,nu1vib)
    piB1vib = [] # (W*cm/m^2) (nlev)
    for k in range(len(z)):
        if not np.isnan(nu1vib[k]):
            val = get_B(nu1vib[k],T[k])
            piB1vib.append(np.pi*val['B'])
        else:
            piB1vib.append(np.nan)
    piB1vib = np.array(piB1vib)
    #############################################################
    # At each height, calculate Beta/z (nlev; do not take tropospheric mean)
    c1    = Gamma*par.Rd/par.ggr
    c2    = par.L/(par.Rv*T)
    phi   = get_phi(par,T)
    phipr = get_phipr(par,T)
    Boverz = (1 + c1*(c2+phi) + c1/(c2+phi)*(c2-phipr))*par.ggr/(par.Rd*T)
    #############################################################
    # Compute the emitting width at the emitting wavenumber
    lvib1 = []
    for k in range(len(z)):
        # wavenumber doing the emitting at this height
        nu1 = nu1vib[k]
        if not np.isnan(nu1):
            # use functional approximation of lvib
            #ltmp = get_lvib(nu1,nuvib,Pvib_lo,Avib_lo,Pvib_co,Avib_co,lvib0)
            lvib1.append(par.lvib0)
        else:
            lvib1.append(np.nan)
    lvib1 = np.array(lvib1)
    #############################################################
    # Compute the SSM1D with varying emitting width
    Hvib = (1/(cp*rho))*(Boverz)*(piB1vib)*(lvib1) # (W/m3->K/s)
    #############################################################
    return {'nu1vib':nu1vib,'z1vib':z,'piB1vib':piB1vib,'Hvib':Hvib,'Boverz':Boverz,'lvib':lvib1}


