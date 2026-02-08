
import numpy as np
from scipy.interpolate import interp1d
from dataclasses import dataclass
from typing import Optional, Callable, Union, Dict

R_UNIV = 8.314462618  # J/mol/K

@dataclass
class Gas:
    name: str
    M: float            # molar mass (kg/mol)
    cpmol: float        # molar heat capacity (J/mol/K) at your working T range
    # Dry-air mole fraction specification (relative to dry air, i.e., excluding H2O):
    # May be: float (scalar), callable f(z)->array, or array of len(z). None if fill_remainder=True.
    xdry: Optional[Union[float, np.ndarray, Callable[[np.ndarray], np.ndarray]]] = None
    # If True, this gas takes the remaining dry-air fraction at each level.
    fill_remainder: bool = False

    @property
    def cp_mass(self) -> float:
        return self.cpmol / self.M  # J/kg/K

    @property
    def R_specific(self) -> float:
        return R_UNIV / self.M      # J/kg/K


def get_custom_atm(par, vres1=np.arange(0, 3e4, 1e2), gases: Optional[Dict[str, Gas]] = None):
    """
    Flexible atmospheric composition model with dry-air bookkeeping.
    - H2O is computed from RH & saturation vapor pressure as in your original.
    - All other gases are specified as dry-air mole fractions and scaled by (1 - xH2O).
    - Supports adding gases like O3 via constant dry ppv, a profile function, or arrays.

    Returns:
        dat1: dict with T, p, rho, RH, cp, cpmol, Gamma, z, alpha_lt, alpha_gt,
              x (dict of per-gas mole fractions vs z), and backward-compatible xN2/xO2/xCO2/xH2O if present.
    """

    #############################################################
    # NEW: Build a default gas registry if none is provided
    #############################################################
    if gases is None:
        # Default set: H2O is computed separately; all others are "dry" gases.
        gases = {
            # This one will be computed from RH and saturation; do not set xdry for H2O.
            'H2O': Gas('H2O', M=par.MH2O, cpmol=par.cpmolH2O),

            # Dry-air gases (define their dry mole fractions). You can add more later.
            # Make one gas the "fill" so you don't need to re-normalize when adding others.
            'N2':  Gas('N2',  M=par.MN2,  cpmol=par.cpmolN2,  xdry=None,        fill_remainder=True),
            'O2':  Gas('O2',  M=par.MO2,  cpmol=par.cpmolO2,  xdry=par.xO2),
            'CO2': Gas('CO2', M=par.MCO2, cpmol=par.cpmolCO2, xdry=par.xCO2),
            # Example if you already have Argon in par: uncomment and adjust dry fraction
            # 'Ar':  Gas('Ar',  M=par.MAr,  cpmol=par.cpmolAr,  xdry=0.00934),
        }

    # Helper to evaluate a gas's dry fraction at the model levels z
    def _eval_xdry(g: Gas, z: np.ndarray) -> Optional[np.ndarray]:
        if g.xdry is None:
            return None
        if callable(g.xdry):
            arr = np.asarray(g.xdry(z), dtype=float)
        else:
            arr = np.asarray(g.xdry, dtype=float)
            if arr.ndim == 0:  # scalar
                arr = np.full_like(z, float(arr), dtype=float)
            elif arr.shape != z.shape:
                raise ValueError(f"xdry array for gas {g.name} must have shape {z.shape}, got {arr.shape}")
        # Bound to [0,1]
        return np.clip(arr, 0.0, 1.0)

    #############################################################
    # Temperature profile & tropopause
    #############################################################
    Gamma = np.ones([len(par.z)], dtype='float') * par.Gamma
    T = par.Ts - par.Gamma * par.z
    mask = np.where(T < par.Ttrp)
    T = np.where(T < par.Ttrp, par.Ttrp, T)
    Gamma[mask] = 0
    ktrp = int(np.amin(np.where(T == par.Ttrp)))  # index of tropopause
    ztrp = par.z[ktrp]

    #############################################################
    # Relative humidity profile (uniform or optimal) like before
    #############################################################
    if par.uniform:
        RH = np.ones([len(par.z)]) * par.RHmid
        alpha_lt = 0
        alpha_gt = 0
        print('initializing with uniform RH')
    else:
        RH = np.ones([len(par.z)]) * par.RHtrp
        foo = get_optimal_RH(T[0:(ktrp + 1)], par.Ts, par.Tmid, par.Ttrp, par.RHs, par.RHmid, par.RHtrp)
        RH[0:(ktrp + 1)] = foo['RH']
        alpha_lt = foo['alpha_lt']
        alpha_gt = foo['alpha_gt']
        print('initializing with non-uniform RH')

    #############################################################
    # Initialize state arrays
    #############################################################
    nlev = len(par.z)
    rho  = np.zeros_like(par.z, dtype=float)
    p    = np.zeros_like(par.z, dtype=float)
    p[0] = par.ps + RH[0] * get_esat_over_l(par, par.Ts)  # surface dry + vapor

    # Per-gas storage (mole fractions)
    x_all = {name: np.zeros(nlev, dtype=float) for name in gases.keys()}
    
    # Per-gas storage (mass fractions)
    w_all = {name: np.zeros(nlev, dtype=float) for name in gases.keys()}

    # Bulk properties
    Mave  = np.zeros(nlev, dtype=float)
    Rtot  = np.zeros(nlev, dtype=float)
    cp    = np.zeros(nlev, dtype=float)      # J/kg/K
    cpmol = np.zeros(nlev, dtype=float)      # J/mol/K

    # Pre-evaluate dry profiles for non-water gases (if callable/scalar/array)
    z = par.z
    dry_names = [g for g in gases.keys() if g != 'H2O']
    xdry_map: Dict[str, Optional[np.ndarray]] = {}
    fill_name: Optional[str] = None
    for name in dry_names:
        g = gases[name]
        xdry_map[name] = _eval_xdry(g, z)
        if g.fill_remainder:
            if fill_name is not None:
                raise ValueError("Only one gas can have fill_remainder=True.")
            fill_name = name
    if fill_name is None:
        raise ValueError("One dry gas must be designated fill_remainder=True (e.g., N2).")

    #############################################################
    # March upward
    #############################################################
    tol = 1e-8
    for k in range(nlev):
        dz = (par.z[k + 1] - par.z[k]) if (k < nlev - 1) else (par.z[1] - par.z[0])

        # --- Water vapor (H2O) from RH & CC, with stratospheric "freeze" at tropopause
        if k <= ktrp:  # troposphere
            pH2O = RH[k] * get_esat_over_l(par, T[k])
            x_all['H2O'][k] = max(pH2O / max(p[k], 1e-6), 1e-12)
        else:
            x_all['H2O'][k] = max(x_all['H2O'][ktrp], 1e-12)
            RH[k] = x_all['H2O'][k] * p[k] / max(get_esat_over_l(par, T[k]), 1e-6)

        # --- Dry-air bookkeeping at level k
        xH2O = x_all['H2O'][k]
        dry_total = max(1.0 - xH2O, 0.0)

        # Sum specified dry fractions at this level
        sum_spec = 0.0
        for name in dry_names:
            if name == fill_name:
                continue
            xi = xdry_map[name][k] if xdry_map[name] is not None else 0.0
            sum_spec += xi

        if sum_spec > 1.0 + 1e-10:
            raise ValueError(f"At z={z[k]:.1f} m the sum of specified dry fractions exceeds 1: {sum_spec:.6f}")

        # Assign dry fractions (scaled by dry_total)
        for name in dry_names:
            if name == fill_name:
                continue
            xi_dry = xdry_map[name][k] if xdry_map[name] is not None else 0.0
            x_all[name][k] = dry_total * xi_dry

        # Fill the remainder to the designated gas
        x_all[fill_name][k] = max(dry_total * (1.0 - sum_spec), 0.0)

        # --- Consistency checks
        x_sum = sum(x_all[name][k] for name in x_all.keys())
        if abs(x_sum - 1.0) > 1e-6:
            # Small numerical drift can happen; renormalize gently
            for name in x_all.keys():
                x_all[name][k] /= x_sum

        # --- Derived bulk properties at level k
        Mave[k]  = sum(x_all[name][k] * gases[name].M      for name in x_all.keys())
        cpmol[k] = sum(x_all[name][k] * gases[name].cpmol  for name in x_all.keys())

        # Mass fractions
        w = {name: x_all[name][k] * gases[name].M / Mave[k] for name in x_all.keys()}
        
        # save the per-gas mass fractions at this level
        for name in w_all.keys():
            w_all[name][k] = w[name]

        # cp and R (mass-based)
        cp[k]   = sum(w[name] * gases[name].cp_mass      for name in x_all.keys())
        Rtot[k] = sum(w[name] * gases[name].R_specific   for name in x_all.keys())

        # Density and pressure integration
        rho[k] = p[k] / (Rtot[k] * T[k])
        arg = -par.ggr / Rtot[k] * dz / T[k]
        if k < nlev - 1:
            p[k + 1] = p[k] * np.exp(arg)

    #############################################################
    # Interpolation to requested vertical resolution
    #############################################################
    def interpolate(var, vres):
        return interp1d(par.z, var, kind='cubic')(vres)

    T1      = interpolate(T, vres1)
    Gamma1  = interpolate(Gamma, vres1)
    p1      = interpolate(p, vres1)
    rho1    = interpolate(rho, vres1)
    RH1     = interpolate(RH, vres1)
    cp1     = interpolate(cp, vres1)
    cpmol1  = interpolate(cpmol, vres1)

    # Interpolate all gas mole-fraction profiles
    x1 = {}
    for name, arr in x_all.items():
        x1[name] = interpolate(arr, vres1)
        
    
    # NEW: Interpolate all gas mass-fraction profiles
    w1 = {}
    for name, arr in w_all.items():
        w1[name] = interpolate(arr, vres1)

    #############################################################
    # Backward-compatible keys for the common gases if present
    #############################################################
    dat1 = {
        'T': T1, 'p': p1, 'Gamma': Gamma1, 'rho': rho1, 'z': vres1, 'RH': RH1,
        'cp': cp1, 'cpmol': cpmol1,
        'alpha_lt': alpha_lt, 'alpha_gt': alpha_gt,
        'x': x1, 'w':w1,  # NEW: general access to all gases
        'Mave': interpolate(Mave, vres1),
    }

    # Populate legacy fields if those gases exist
    for legacy in ['N2', 'O2', 'H2O', 'CO2']:
        if legacy in x1:
            dat1[f'x{legacy}'] = x1[legacy]
        if legacy in w1:
            dat1[f'w{legacy}'] = w1[legacy]

    return dat1


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
    
def get_esat_over_l(par,T):
    import math
    # SVP over liquid from DAM (Pa)
    xi = 1 # factor to change the SVP
    pstar = xi*par.ptrip * (T/par.Ttrip)**((par.cpv-par.cvl)/par.rgasv) * math.exp( (par.E0v - (par.cvv-par.cvl)*par.Ttrip) / par.rgasv * (1/par.Ttrip - 1/T) )
    return pstar