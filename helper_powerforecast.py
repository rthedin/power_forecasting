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
