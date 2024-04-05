import numpy as np
import pandas as pd
import xarray as xr
import os, sys
import re
import multiprocessing 
import warnings
import itertools

from timezonefinder import TimezoneFinder
import datetime
import time

import matplotlib.pyplot as plt
from herbie import Herbie# , FastHerbie
from toolbox import EasyMap, pc, ccrs

from helper_powerforecast import *
from description_all_vars import desc_sfc, desc_subh


class powerForecast:

    def __init__ (self,
                  loc,
                  datetime_oi_local = None,
                  datetime_oi_utc = None,
                  fxx_period = '18h',
                  model = 'hrrr',
                  product = 'sfc',
                  avgAroundLoc = False,
                  verbose=1):

        '''
        Initialize the power forecast object

        Inputs:
        ------
        loc: tuple of scalars
            Location of interest given in (long, lat)
        datetime_oi_local, datetime_oi_utc: np.datetime64
            Datatime of interest, either in the local timezone or in
            UTC. Only one of the two should be provided.
        fxx_period: string or int
            Period of forecast, starting at the datetime. Given
            as a time string ('1d', '4h', '120min') or integer in hours
        model: str
            Model used for forecasting. Many available, only hrrr currently
            tested.
        product: str
            Product data. Surface, pressure, nativel and subhourly products
            are available.
        avgAroundLoc: bool (default False)
            Whether or not to average the resulting values around the
            location of interest with a 3x3 grid of points
        verbose: scalar
            Verbosity level

        Useful methods
        --------------
        show_all_variables()

        '''

        self.loc               = loc
        self.datetime_oi_local = datetime_oi_local
        self.datetime_oi_utc   = datetime_oi_utc
        self.fxx_period        = fxx_period
        self.model             = model
        self.product           = product
        self.avgAroundLoc      = avgAroundLoc
        self.verbose           = verbose


        self._check_inputs()

        self._determine_all_valid_fxx()




    def _determine_all_valid_fxx(self):
        '''
        Using the fxx_period requested by the user, create an array of desired forecast values that
        are also valid and exist.
        Forecast values are:
        fxx = 0 through 48, step=1: availble for available for runs initialized at 00z, 06z, 12z, 18z
        fxx = 0 through 18, step=1: available for runs initialized at all other hours.
        '''

        # Determine available forecast hours given start time   
        if self.datetime_oi_tz.hour in (0,6,12,18):
            available_fxx_hours = np.arange(0,48)+1
        else:
            available_fxx_hours = np.arange(0,18)+1
            
            
        # Parse the requested forecast hours
        if isinstance(self.fxx_period, str):
            try:
                self.fxx_period_delta = pd.Timedelta(self.fxx_period)
            except Exception as e:
                print(e)
                raise ValueError(f"Temporal unit should be given, like `h`, `d`, `min`. Error: {str(e)}")
            self.fxx_period = self.fxx_period_delta.total_seconds()/3600
            
        if isinstance(self.fxx_period, float):
            if self.fxx_period%1.0 != 0:
                raise ValueError(f"Forecast only available for integer number of hours. Received {self.fxx_period} hours.")
            self.fxx_period = int(self.fxx_period)
            
        desired_fxx_hours = np.arange(0,self.fxx_period)+1

        if not set(desired_fxx_hours).issubset(set(available_fxx_hours)):
            raise ValueError(f"Not all desired forecast hours are available.\nAvailable: {available_fxx_hours}\nRequested: {desired_fxx_hours}")

        # Add 0-forecast (analysis)
        self.desired_fxx_hours = np.insert(desired_fxx_hours, 0, 0)


    def _determine_all_valid_datetime(self):

        start = self.datetime_oi_tz - pd.Timedelta(self.fxx_period, 'h')
        end   = self.datetime_oi_tz + pd.Timedelta(self.fxx_period, 'h')

        self.desired_datetime_period_tz =  pd.date_range(start, end=end, freq='1h').tolist()


    def _check_inputs(self):
    
        # Operations on loc, checking if CONUS. Need to relax check if offshore
        top = 49.3457868 # north lat
        left = -124.7844079 # west long
        right = -66.9513812 # east long
        bottom =  24.7433195 # south lat
        if not bottom <= self.loc[1] <= top and not left <= self.loc[0] <= right:
            raise ValueError(f"Location given is not within the CONUS")
    
    
        # Operations on datetimes
        datetime_oi_local = pd.to_datetime(self.datetime_oi_local)
        datetime_oi_utc   = pd.to_datetime(self.datetime_oi_utc)
    
        if datetime_oi_local is None and datetime_oi_utc is not None:
            self.datetime_oi_tz = datetime_oi_utc.tz_localize('UTC')
        elif datetime_oi_local is not None and datetime_oi_utc is None:
            tz = TimezoneFinder().timezone_at(lng=self.loc[0], lat=self.loc[1])
            self.datetime_oi_tz = datetime_oi_local.tz_localize(tz)
        else:
            raise ValueError(f'One and only one of datetime_oi_local and datetime_oi_utc should be given.')
    
    
    
        # Operations on product
        if self.product == 'sfc':
            self.description = desc_sfc
        elif self.product == 'subh':
            self.description = desc_subh
        else:
            raise ValueError('No description table for {self.product} product.')
    
    
    
    
    def _get_desired_hrrr_variables(self):
    
        # Several relevant variables. User should add to list desired variables.
        # Allows regex (typical wildcard * should be .*, where the dot is any character).
        desired_vars = ['UGRD', 'VGRD', 'WIND', 'TMP', 'HPBL', 'GUST',
                        'HGT:.*mb', 'HGT:surf', 'MAXUW', 'MAXVW', 'SFCR', 'FRICV',
                        'SHTFL', 'LHTFL', 'GFLUX', 'WRF:surf', '.*CSH:0-1000','POT']
    
        # search the strings above in the dict and create a list
        desired_hrrr_vars = []
    
        for curr_var in desired_vars:
            matching_vars = [k for k in list(self.description.keys()) if re.search(curr_var, k)]
            desired_hrrr_vars.append(matching_vars)
        if self.verbose >1: print(f"{ingreen(curr_var)}: {matching_vars}")
    
        # Flatten the list
        self.desired_hrrr_vars = [i for sublist in desired_hrrr_vars for i in sublist]
    
    
    
    
    def get_HRRR(self):
        '''
        Get all grib2 files from HRRR given the start date and all fxx values requested
        '''
    
    
        # First get the current non-forecast time
        H = _get_single_HRRR(curr_fxx = 0)
    
    
        # Now, loop over the requested forecast, if any
        for curr_fxx in self.desired_fxx_hours:
            H = _get_single_HRRR(curr_fxx)
    
    
    
    def _get_single_HRRR(self, curr_fxx):
        '''
        For all the desired values, open the grib file and get the value of interest, adding to a dictionary.
        Several warnings are issues by herbie; let's silence those where we know it's okay and then re-enable them.
        '''
    
        change_warning('ignore')
    
        # Open grib file. Herbie cannot handle timezones, so we will convert it manually
        datetime_oi_hrrr = self.datetime_oi_tz.tz_convert(tz=None)  # converts to UTC and makes it timezone-naive
        H = Herbie(datetime_oi_hrrr, model=self.model, product=self.product, fxx=curr_fxx, verbose=False)
    
        H.download()
        
        # Create empty dataset
        all_values = xr.Dataset.from_dict( {
            "coords": {
                "datetime": {"dims": "datetime",
                             "data": [self.datetime_oi_tz]},
                "fxx":      {"dims": "fxx",
                             "data": [curr_fxx]},
            },
            "dims": {"datetime", "fxx"},
            "data_vars": {
            },
        } )
    
        for var in self.desired_hrrr_vars:
            ds = H.xarray(var, remove_grib=False)
            dsi = ds.herbie.nearest_points(self.loc)
            # In dsi, there are two data variables: the one we're interested in, and `gribfile_projects`. Let's get the relevant one.
            alldatavars = list(dsi.data_vars)
            var_oi =  [i for i in alldatavars if i != 'gribfile_projection']
            if len(var_oi) > 1:
                print(f"WARNING: Variable {var} has {len(var_oi)} data entries: {var_oi}. Getting {var_oi[0]}")
            var_oi = var_oi[0]
            value = dsi[var_oi].values[0]
            # Add to dataset
            all_values[var] = (('datetime', 'fxx'), [[value]])
    
        change_warning('default')
    
        return all_values
    
    
    
    
    def _get_multiple_HRRR_par(list_of_datetime_fxx):
        
        t0 = time.time()
    
        ncores = multiprocessing.cpu_count()
    
        print(f"Getting {len(list_of_datetime_fxx)} GRIB files using {ncores} cores.")
        if __name__ == '__main__':
            pool = multiprocessing.Pool(processes=ncores)
            allds = pool.starmap(self._get_single_HRRR, list_of_datetime_fxx)
            pool.close()
            pool.join()
    
        ds = xr.merge(allds)
    
        tf = time.time()
        if verbose>0:
            print(f'Done. Finished in {time.strftime("%Hh%Mm%Ss", time.gmtime(tf-t0))}.')
        
        ds = calc_qois(ds)
        return ds
    
    
    
    
    def _get_all_combinations_datetime_fxx(self):
        list_of_datetime_fxx_allcomb = list(itertools.product(self.desired_datetime_period_tz, self.desired_fxx_hours))
        ds = self._get_multiple_HRRR_par(list_of_datetime_fxx_allcomb)
    
    
    def _get_forecast_only(self):
        list_of_datetime_fxx_fxx = list(itertools.product([self.datetime_oi_tz], self.desired_fxx_hours))
        ds_fxx = self._get_multiple_HRRR_par(list_of_datetime_fxx_fxx)
    
    
    def _get_analysis_only(self, analysis_len_in_h=24*10):
        all_analysis =  pd.date_range(self.datetime_oi_tz, end=self.datetime_oi_tz + pd.Timedelta(analysis_len_in_h, 'h'), freq='1h').tolist()
        list_of_datetime_fxx_anl = list(itertools.product(all_analysis, [0]))
        ds_anl = self._get_multiple_HRRR_par(list_of_datetime_fxx_anl)
    
    
    
    
    
    
    
    
    
    
    
    def show_all_variables():
        
        # Get a hrrr file just to have the variables
        H = Herbie(self.datetime_oi_tz.tz_convert(tz=None), model=self.model, product=self.product, fxx=0, verbose=False)
    
        vars_all = list(H.inventory(searchString=".")['search_this'])
        # Variables represented as :SYMBOL:level:analysis_or_fxx. Let's drop the anl/fxx part
        vars_all = [':'.join(i.split(':')[:3]) for i in vars_all]
        param = [self.description[v]['parameter'] for v in vars_all]
        desc =  [self.description[v]['description'] for v in vars_all]
        level = [self.description[v]['level'] for v in vars_all]
    
        print("AVAILABLE VARIABLES")
        print("-------------------")
        for i, p in enumerate(param):
            print(f"- {ingreen(f'{desc[i]} at {level[i]}')}. Search for {vars_all[i]}:.*")
            




if __name__ == "__main__":

    # --------------- USER INPUT ----------------
    loc = (-105, 39)
    datetime = "2023-01-01 00:00:00"
    fxx_period = '24h'
    avgAroundLoc = False
    # ---------- END OF USER INPUT --------------

    pf = powerForecast(loc, datetime, fxx_period, avgAroundLoc)

    # pf.get_HRRR()

