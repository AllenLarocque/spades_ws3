# This setups up ws3 and loads various default ws3 parameters
# Defaults are overwritten as defined in spades_ws3.R. Indirectly some of these may be defined in global.R.

import os
from os import listdir
from os.path import isfile, join, dirname
import sys

if enable_debugpy:
    import debugpy
    try:
        debugpy.listen(("localhost", 5678))
    except RuntimeError:
        # Already listening
        pass

    if not debugpy.is_client_connected():
        print("Waiting for debugger attach...")
        debugpy.wait_for_client()
    
########################################################################################################
import ws3
from ws3.forest import ForestModel, Action
from ws3.spatial import ForestRaster
from ws3.common import clean_vector_data, reproject_vector_data, rasterize_stands, hash_dt, warp_raster
########################################################################################################

import numpy as np
import pandas as pd
from pandas import DataFrame as DF
import geopandas as gpd
import matplotlib.pyplot as plt
import rasterio
import rasterio.plot
import seaborn as sns
from functools import partial
import fiona
import folium
import pathlib
try:
   import cPickle as pickle
except:
   import pickle
import copy
from pathlib import Path
import rasterio
import shutil

from spadesws3 import clean_shapefiles, rasterize_inventory, read_basenames, compile_basecodes, bootstrap_forestmodel, bootstrap_areas, schedule_harvest_optimize, schedule_harvest_areacontrol, sda

# configure paths and global default variables
scenario_name = 'base'
sda_mode = 'randblk'       # 'randpxl'
obj_mode = 'min_harea'     # 'max_harea'
horizon = 10
period_length = 10
yields_period_length = 10
yields_x_unit = 'years'
time_step = 1
max_age = 1000
save_fm = False
target_path = join(dat_path, 'targets.csv')
yld_path = dat_path
tolerance = 10.
clean_inv = False
rasterize_inv = False
gdb_path = lambda bn: '%s/gis/gdb/%s.gdb' % (dat_path, bn)
shp_path = lambda bn: '%s/gis/shp/%s.shp' % (dat_path, bn)
tif_path = lambda bn: '%s/tif/%s' % (dat_path, bn)
shp_name = 'stands'
age_col = 'age'
theme_cols = ['theme0', 'theme1', 'theme2', 'theme3']
compress = 'lzw'
dtype = rasterio.uint8
base_year = 2015
prop_names = [u'THLB', u'AU', u'LdSpp', u'Age2015', u'Shape_Area']
prop_types = [(u'theme0', 'str:10'),
              (u'theme1', 'str:1'),
              (u'theme2', 'str:5'), 
              (u'theme3', 'str:50'), 
              (u'age', 'int:5'), 
              (u'area', 'float:10.1')]
tvy_name = 'totvol'
snk_epsg = 3005 # ESPG:3005 corresponds to NAD83/BC Albers
sns.set_style('dark')
hdt_path = '%s/hdt' % dat_path
coast_bn = ['tsa01', 'tsa02', 'tsa03'] # FIX ME: totally bogus (check GIS data for real coast TSA codes) 
oe_harvest = '_age >= 100 and _age <= 400'
action_params = {'harvest':{'oe':oe_harvest,
                            'mask':('?', '1', '?', '?'),
                            'is_harvest':True,
                            'targetage':0}}
util = 0.85
proc_mode = 'seq_tsa'
deg_mode = 'classic'
obj_mode = 'min_harea'
run_sda = 1 # set to True for production runs
fm_horizon = horizon
cap_age = 900
twopass_cacut_factor = 1.50
run_sda = True
raster_d = 250
hdt = {bn:pickle.load(open('%s/hdt_%s.pkl' % (hdt_path, bn), 'rb')) for bn in basenames}

def kwargs():
    basecodes = compile_basecodes(hdt, basenames, theme_cols)
    kwargs = {'basenames':basenames,
              'model_name':'foo',
              'model_path':dat_path,
              'base_year':int(base_year),
              'yld_path':yld_path,
              'tif_path':tif_path,
              'horizon':int(horizon),
              'period_length':int(period_length),
              'max_age':int(max_age),
              'basecodes':basecodes,
              'action_params':action_params,
              'hdt':hdt,
              'add_null_action':True,
              'tvy_name':tvy_name,
              'compile_actions':True,
              'yields_x_unit':yields_x_unit,
              'yields_period_length':int(yields_period_length),
              'verbose':False}
    return kwargs    


def bootstrap_forestmodel_kwargs():
    return bootstrap_forestmodel(**kwargs())


def _save_forestmodel(fm, year, verbose=False):
    """
    Save forest model object using pickle, clearing solver models first.
    Solver models (_model attributes in Problem objects) cannot be pickled as they appeal to C++ libraries etc.
    
    """
    # Get output_path from globals (set by R via py$output_path)
    base_path = globals().get('output_path', dat_path)
    
    # Save to fm_checkpoints subdirectory in output path
    save_path = join(base_path, 'fm_checkpoints')
    os.makedirs(save_path, exist_ok=True)
    filename = join(save_path, 'fm_year_%04d.pkl' % year)
    
    if verbose:
        print('Saving forest model to: %s' % filename)
    
    def clear_solver_models(obj, visited=None):
        '''Recursively clear _model attributes from Problem objects to enable pickling'''
        if visited is None:
            visited = set()
        obj_id = id(obj)
        if obj_id in visited:
            return
        visited.add(obj_id)
        
        # Clear _model attribute if present
        if hasattr(obj, '_model') and obj._model is not None:
            obj._model = None
        
        # Recursively process dictionaries, lists, and object attributes
        try:
            if isinstance(obj, dict):
                for v in obj.values():
                    clear_solver_models(v, visited)
            elif isinstance(obj, (list, tuple)):
                for v in obj:
                    clear_solver_models(v, visited)
            elif hasattr(obj, '__dict__'):
                for attr_value in obj.__dict__.values():
                    clear_solver_models(attr_value, visited)
        except (RuntimeError, RecursionError):
            pass  # Skip if recursion limit reached
    
    # Clear solver models before saving (they'll be recreated when needed)
    clear_solver_models(fm)
    
    # Save the forest model
    with open(filename, 'wb') as f:
        pickle.dump(fm, f)
    
    if verbose:
        print('Saved forest model to: %s' % filename)


def simulate_harvest(fm, basenames, year,
                     planning_period_freq = 10, # Allen added
                     mode='optimize',
                     target_scalefactors=None,
                     mask_area_thresh=0.,
                     verbose=False,
                     mgmt_unit_theme=None,
                     workers=1,
                     save_fm=None):
    # Check module-level save_fm (set by R via py$save_fm), then use parameter if provided
    try:
        import sys
        current_module = sys.modules[__name__]
        global_save_fm = getattr(current_module, 'save_fm', False)
    except:
        global_save_fm = globals().get('save_fm', False)
    
    # Use parameter if provided, otherwise use global
    if save_fm is None:
        save_fm = global_save_fm
    elif not save_fm:  # If parameter is False/0, still check global
        save_fm = global_save_fm
    
    # Ensure save_fm is a boolean (handle R logical values)
    if not isinstance(save_fm, bool):
        if save_fm is None:
            save_fm = False
        elif isinstance(save_fm, (int, float)):
            save_fm = bool(int(save_fm)) if save_fm != 0 else False
        else:
            save_fm = bool(save_fm)
    
    # Convert workers to int (R numeric values come through as float)
    workers = int(workers) if workers is not None else 1
    
    if verbose or save_fm:
        print('simulate_harvest: save_fm=%s, year=%d' % (save_fm, year))
    
    bootstrap_areas(fm, basenames, tif_path, hdt, year, new_dts=False)
    fm.reset()

    # Only schedule harvests at year 0 or multiples of period_length (Allen added)
    if year == 0 or year % planning_period_freq == 0:  # '%' is the modulo operator. It returns the remainer of a division between the two numbers

      if mode == 'optimize':
           schedule_harvest_optimize(fm, basenames, p_max_hv=target_scalefactors,
                                     mgmt_unit_theme=mgmt_unit_theme, workers=workers)
          #profile_schedule_harvest_optimize(fm, basenames, target_scalefactors,
          #                                  mgmt_unit_theme, workers)
      elif mode == 'areacontrol':
          schedule_harvest_areacontrol(fm,
                                       target_scalefactors=target_scalefactors,
                                       mask_area_thresh=mask_area_thresh,
                                       verbose=verbose)
      else: # bad mode value
          raise ValueError('Undefined optimizer mode. "optimize" and "areacontrol" should work')
      
      # Save forest model after optimizer runs if enabled
      if save_fm:
          try:
              _save_forestmodel(fm, year, verbose)
              if verbose:
                  print('Successfully saved forest model for year %d' % year)
          except Exception as e:
              print('ERROR saving forest model for year %d: %s' % (year, str(e)))
              if verbose:
                  import traceback
                  traceback.print_exc()
              # Don't raise - allow simulation to continue
    
    ## After harvests are scheduled, spatialize them by running SDA (Spatial Disturbance Allocator)
    # This connects the aspatial harvests to the raster map
    # SDA runs once at year 0 and then every multiple of planning_period_freq
    if year == 0 or year % planning_period_freq == 0:  # Allen added. Runs SDA at year 0 and also every multiple of the planning_period_freq
      sda(fm, basenames, 1, tif_path, hdt, sda_mode=sda_mode, verbose=verbose)
    
    # SDA runs every year
    #sda(fm, basenames, 1, tif_path, hdt, sda_mode=sda_mode, verbose=verbose)


# Define optimizer speed profiler function (UNTESTED):
def profile_schedule_harvest_optimize(fm, basenames, target_scalefactors, 
                                      mgmt_unit_theme, workers):
    """
    Profile schedule_harvest_optimize() and print the 20 slowest functions.
    """
    import cProfile
    import io
    import pstats
    pr = cProfile.Profile()
    pr.enable()

    # ---- Run the actual function ----
    result = schedule_harvest_optimize(fm, basenames, p_max_hv=target_scalefactors, 
                                       mgmt_unit_theme=mgmt_unit_theme, workers=workers)

    pr.disable()

    # ---- Print top 20 functions by cumulative time ----
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats('cumtime')
    ps.print_stats(20)
    print(s.getvalue())

    return result
