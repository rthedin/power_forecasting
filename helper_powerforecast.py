import numpy as np
import xarray as xr

class color():
    RED = "\033[31m"
    GREEN = "\033[32m"
    BLUE = "\033[34m"
    RESET = "\033[0m"

def printred(*args, **kwargs):
    print( color.RED + " ".join(map(str,args))+ color.RESET, **kwargs)

def printgreen(*args, **kwargs):
    print( color.GREEN + " ".join(map(str,args))+ color.RESET, **kwargs)

def printblue(*args, **kwargs):
    print( color.BLUE + " ".join(map(str,args))+ color.RESET, **kwargs)

def inred(string):
    '''
        Transform text to red
        Example call:
            print(f"Black {inred('this should be red')}")
            a=1.3213
            print(f"This is black {inred(f'{a:.2f}')}, and the number is red")
    '''
    if not isinstance(string,str):
        string = str(string)
    return color.RED + string + color.RESET


def ingreen(string):
    if not isinstance(string,str):
        string = str(string)
    return color.GREEN + string + color.RESET


def change_warning(action):
    '''
    Silencing and re-enabling of custom warnings for convenience
    
    Input:
    ------
    action: str
        Following warnings.filterwarnings. Typical use: "ignore" or "default"
    '''
    
    warnings.filterwarnings(action, message="More than one time coordinate present for variable")
    warnings.filterwarnings(action, message="Calling float on a single element Series is deprecated and will raise a TypeError in the future")
    warnings.filterwarnings(action, message="Will not remove GRIB file because it previously existed.")



def calc_qois(ds):
    ''' Calculates some quantities of interest '''
    
    # ----------- Wind speeds at 10 and 80 m
    ds['_wspd at 80 m'] = (ds[':UGRD:80 m above ground']**2 + ds[':VGRD:80 m above ground']**2)**0.5
    ds['_wspd at 10 m'] = (ds[':UGRD:10 m above ground']**2 + ds[':VGRD:10 m above ground']**2)**0.5
    
    # ----------- Heat flux
    # Net radiation (positive into the surface, negative out of the surface). All the long/shortwaves are positive in value from the ds
    #         down shortwave  +  up shortwave        +   down longwave      +  up longwave        
    Rn = ds[':DSWRF:surface'] - ds[':USWRF:surface'] + ds[':DLWRF:surface'] - ds[':ULWRF:surface'] 
    ds['Rn'] = Rn
    # Sensible heat flux (positive when energy into the domain, out of the surface)
    H = ds[':SHTFL:surface']
    # Latent heat flux (positive when energy into the domain, out of the surface)
    LS = ds[':LHTFL:surface']
    # Soil heat flux (flipping it so positive when energy into the domain, out of the surface)
    G = ds[':GFLUX:surface']
    # Total budget (positive out of the surface)
    ds['_energybudget W/m2'] = H + LS - G 
    
    return ds

