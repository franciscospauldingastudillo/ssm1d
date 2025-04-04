# SSM1D: One-Dimensional Simple Spectral Model 

**SSM1D** is a lightweight, idealized model of radiative cooling based on spectral water vapor absorption features. It solves an analytic form of the cooling-to-space approximation of the full longwave radiative transfer equation under clear-sky conditions.

This package is designed to work in tandem with [`parameters-class`](https://github.com/franciscospauldingastudillo/parameters-class), a utility package for defining thermodynamic profiles and model parameters.

---

## 🔧 Installation

You can install SSM1D directly from GitHub:

```bash
pip install git+https://github.com/franciscospauldingastudillo/ssm1d.git
```

## 📦 Dependencies

You'll need to download the parameters-class package to use the SSM1D. You can install it directly from GitHub:
```bash
pip install git+ssh://git@github.com/franciscospauldingastudillo/ssm1d.git
```

## 🚀 Usage

```python
import ssm1d
import parameters
import rfmtools
```

### Example: Run the SSM1D
```python
# Step 0: Define the case and the resolution
# specify what gases and CIA-pairs to include (see generate_case)
gases     = ['H2O']
ciapairs  = []

# dynamically create an argument dictionary
def generate_args(exp, gases, ciapairs, **thermo):
    return {'exp': exp, **{f'gas{i+1}': gas for i, gas in enumerate(gases)}, 'valid_ciapairs': ciapairs, **thermo}
args = generate_args('earth',gases,ciapairs,RHs=0.75,RHmid=0.54,RHtrp=0.75,uniform=1,Ts=300,Tmid=250,Ttrp=200)
        
# create a class instance and generate an RFM case from argument dictionary (300K = 7.1 K/km; 315K = 5.0 K/km)
par = parameters.parameters(Gamma=6.5e-3)
par.generate_case(**args)

# vertical resolution
RFM      = np.arange(par.z[0],par.z[-1],1e2)

# default dataset
# z ~ m, p ~ Pa, T ~ K, Gamma ~ K/m, RH~unitless, hr~W/m3, srhr~cm*W/m3
# srtau~unitless, zsrtau1~m, crnus~cm-1, crsrtau~unitless, crzsrtau1~m
# piBsr~W*cm/m^2, piBbar~W/m^2
dataset  = ({'RFM':{
                    'z':{},'p':{},'T':{},'rho':{},'Gamma':{},'RH':{}, 'cp':{}, # THERMODYNAMICS
                     'alpha_lt':{}, 'alpha_gt':{}, # RELHUM PARAMS
                    'hr':{},'srhr':{}, # HEATING RATES
                    'srtau':{},'zsrtau1':{},'crnus':{},'crsrtau':{},'crzsrtau1':{}, # OPTICAL DEPTH
                   'piBsr':{},'piBbar':{}}, # BLACKBODY FLUX
             'SSM1D-LLA':{'z':RFM,'hr':{}},
             'SSM1D-QLLA':{'z':RFM,'hr':{}}}) 

# Step 1: Generate custom atmospheric profiles (p,T,z,x) at RFM and RFMi resolution
dat1 = ssm1d.atm.get_custom_atm(par, RFM)
keys = ['p', 'T', 'RH', 'rho', 'Gamma', 'cp', 'cpmol', 'alpha_lt', 'alpha_gt','z']
for key in keys:
    dataset['RFM'][key] = dat1[key]

# dynamically add molar mixing ratios to the dataset (signals to rfmtools how to build the RFM experiments)
for gas in gases:
    xgas_key = f'x{gas}'  # e.g., xN2, xCH4, xH2
    dataset['RFM'][xgas_key] = dat1[xgas_key]

# Step 2: Run the SSM1D for the rotation band
#############################################################
# NEW QLLA-rot
par.update_band('wv-rot-right')
dataset['par'] = par
# Run the SSM1D for the rotation band
dat1 = get_SSM1D_v2(dataset)

# Step 3: Run the SSM1D for the rotation-vibration band
#############################################################
# NEW LLA-VIB
par.update_band('wv-vib-rot')
dataset['par'] = par
dat2 = get_SSM1D_v2(dataset)
```


