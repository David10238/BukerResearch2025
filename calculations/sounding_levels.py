import numpy as np
from metpy_dataframe_loader import Sounding
from metpy.units import units
import metpy.calc as mpcalc

class LCL:
  def makeLCL(sounding:Sounding):
    return LCL(mpcalc.lcl(sounding.pressure * units.mbar, sounding.temperature*units.degC, sounding.dewpoint*units.degC))

  def __init__(self, lcl):
    self.lcl = lcl

class LFC:
  def makeLFC(sounding:Sounding):
    # make sure data exists
    msk = np.isnan(sounding.pressure) | np.isnan(sounding.temperature) | np.isnan(sounding.dewpoint)
    if np.all(msk):
      return None
    

    # ensure pressure is decreasing
    cleaned_pressure = sounding.pressure[~np.isnan(sounding.pressure)]
    if np.any(np.diff(cleaned_pressure) > 0):
      return None
    

    lfc_data = mpcalc.lfc(sounding.pressure * units.mbar, sounding.temperature*units.degC, sounding.dewpoint*units.degC)
    
    # find closest pressure
    pressure = lfc_data[0]
    abs_vals = np.abs(sounding.pressure * units.mbar - pressure)
    if(np.all(np.isnan(abs_vals))):
      return None
    idx = np.nanargmin(abs_vals)
    lfc_height_pressure = sounding.height[idx]

    # find closest temperature
    temperature = lfc_data[1]
    abs_vals = np.abs(sounding.temperature * units.degC - temperature)
    if(np.all(np.isnan(abs_vals))):
      return None
    idx = np.nanargmin(abs_vals)
    lfc_height_temperature = sounding.height[idx]

    # todo maybe fix here
    #assert lfc_height_pressure == lfc_height_temperature, f'My assumptions about lfc method need to be checked. Heights were {lfc_height_pressure} and {lfc_height_temperature}'

    height = lfc_height_pressure
    return LFC(height, pressure, temperature)

  def __init__(self, height:float, pressure:float, temperature:float):
    self.height = height
    self.pressure = pressure
    self.temperature = temperature