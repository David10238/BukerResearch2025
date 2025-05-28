
import pandas as pd
import numpy as np
from metpy.units import units

import metpy.calc as mpcalc

class LFC:
  def makeLFC(sounding):
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

    # todo 
    #assert lfc_height_pressure == lfc_height_temperature, f'My assumptions about lfc method need to be checked. Heights were {lfc_height_pressure} and {lfc_height_temperature}'

    height = lfc_height_pressure
    return LFC(height, pressure, temperature)

  def __init__(self, height, pressure, temperature):
    self.height = height
    self.pressure = pressure
    self.temperature = temperature




# time (second)
# pressure (mbar)
# height (meter)
# temperature (degC)
# dewpoint (degC)
# speed (meter)
class Sounding:
  def find_lifted_condensation_level(self):
    return mpcalc.lcl(self.pressure * units.mbar, self.temperature*units.degC, self.dewpoint*units.degC)

  def __init__(self, df:pd.DataFrame):
    self.time = df['time'].to_numpy()
    self.pressure = df['pressure'].to_numpy()
    self.height = df['height'].to_numpy()
    self.temperature = df['temperature'].to_numpy()
    self.dewpoint = df['dewpoint'].to_numpy()
    self.wind_u = df['wind_u'].to_numpy()
    self.wind_v = df['wind_v'].to_numpy()
    self.speed = df['speed'].to_numpy()
    self.direction = df['direction'].to_numpy()
    self.longitude = df['longitude'].to_numpy()
    self.latitude = df['latitude'].to_numpy()

    self.as_df = df
    
    self.level_free_convection = LFC.makeLFC(self)

def __load_sounding_df(file_name:str)->pd.DataFrame:
  clean_directory = "Full-CP20-sounding-dataset-clean"
  file_path = f'{clean_directory}/{file_name}'

  df = pd.read_csv(file_path)
  df = df.drop(df.columns[0], axis=1)
  return df

def __convert_df_to_metpy(df:pd.DataFrame)->dict[str, np.array]:
  u = df['wind_u'].to_numpy() * (units.meter / units.second)
  v = df['wind_v'].to_numpy() * (units.meter / units.second)

  return pd.DataFrame({
    'time' : df['time'],
    'pressure' : df['pressure'],
    'height' : df['height'],
    'temperature' : df['temperature'],
    'dewpoint' : df['dewpoint'],
    'speed' : mpcalc.wind_speed(u, v),
    'direction' : mpcalc.wind_direction(u, v),
    'longitude' : df['longitude'],
    'latitude': df['latitude'],
    'wind_u' : df['wind_u'],
    'wind_v' : df['wind_v']
  })


def load_and_convert_sounding(file_name:str) -> Sounding:
  return Sounding(__convert_df_to_metpy(__load_sounding_df(file_name)))