import numpy as np
import pandas as pd
import xarray as xr
import datetime
import os, sys

from herbie import Herbie
from toolbox import EasyMap, pc, ccrs

from helper_powerforecast import *
from description_all_vars import desc_sfc, desc_sub h
import wgrid2_variables


class powerForecast:

    def __init__ (self,
                  loc,
                  datetime,
                  fxx_period = '24h',
                  avgAroundLoc = False,
                  verbose=1):

        '''
        Initialize the power forecast object

        Inputs:
        ------

        loc: tuple of scalars
            Location of interest given in (long, lat)
        datetime: np.datetime64
            Datatime of interest
        fxx_period: string or int
            Period of forecast, starting at the datetime. Given
            as a time string ('1d', '4h', '120min') or integer in hours
        avgAroundLoc: bool (default False)
            Whether or not to average the resulting values around the
            location of interest with a 3x3 grid of points
        verbose: scalar
            Verbosity level

        '''

        self.loc          = loc
        self.datetime     = datetime
        self.fxx_period   = fxx_period
        self.avgAroundLoc = avgAroundLoc
        self.verbose      = verbose

        # Hard-coded default values
        self.model = 'hrrr'
        self.product = 'sfc'


        if product == 'sfc':
            self.description = desc_sfc
        elif product == 'subh':
            self.description = desc_subh
        else:
            raise ValueError('No description table for {self.product} product.')


        self._determine_all_valid_fxx()




        def _determine_all_valid_fxx(self):
            '''
            Using the fxx_period requested by the user, create an array of desired forecast values that
            are also valid and exist.
            Forecast values are:
            fxx = 0 through 48, step=1: availble for available for runs initialized at 00z, 06z, 12z, 18z
            fxx = 0 through 18, step=1: available for runs initialized at all other hours.
            '''
            pass




        def check_inputs(self):
            pass




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
            ''' Open grib2 file using Herbie '''

            H = Herbie(self.datetime,
                       model = self.model,
                       product = self.product,
                       fxx = curr_fxx
                       )
            return H


        def _is_location_valid(self):
            pass


        def _get_all_variables_available_in_grib(self):
            '''
            Get a list of all the variables available in the current
            GRIB2 file by using a generic regex search string
            '''




        def show_all_variables(self):
	    vars_all = list(H.inventory(searchString=".")['search_this'])
	    param = [description[v]['parameter'] for v in vars_all]
	    desc =  [description[v]['description'] for v in vars_all]
	    level = [description[v]['level'] for v in vars_all]

	    print("AVAILABLE VARIABLES")
	    print("-------------------")
	    for i, p in enumerate(param):
		print(f"- {ingreen(f'{desc[i]} at {level[i]}')}. Search for {vars_all[i]}")

            




if __name__ == "__main__":

    # --------------- USER INPUT ----------------
    loc = (-105, 39)
    datetime = "2023-01-01 00:00:00"
    fxx_period = '24h'
    avgAroundLoc = False
    # ---------- END OF USER INPUT --------------

    pf = powerForecast(loc, datetime, fxx_period, avgAroundLoc)

    # pf.get_HRRR()

