import pandas as pd
import numpy as np
from metpy.units import units
import metpy.calc as mpcalc
from metpy.calc import cape_cin, parcel_profile, bulk_shear

class Shear:
  def __init__(self, u:float, v:float):
    self.u = u
    self.v = v

class StormEnergy:
  def __init__(self, cape:float, cin:float):
    self.cape = cape
    self.cin = cin

class LFC:
  def __init__(self, height:float, pressure:float, temperature:float):
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
  
  def calculate_energy(self, maxHeight:int = -1)-> StormEnergy | None:
    if(maxHeight < 0):
      maxHeight = self.height.max()    

    H = self.height
    msk = H <= maxHeight

    p = self.pressure
    p = p[msk]

    T = self.temperature
    T = T[msk]

    Td = self.dewpoint
    Td = Td[msk]

    # metpy is fussy if this mask isn't applied
    msk = np.isnan(p) | np.isnan(T) | np.isnan(Td)
    p = p[~msk] * units.mbar
    T = T[~msk] * units.degC
    Td = Td[~msk] * units.degC

    try:
      prof = parcel_profile(p, T[0], Td[0]).to('degC')
      energy = cape_cin(p, T, Td, prof)
      return StormEnergy(energy[0], energy[1])
    except:
      return None
  
  def calculate_bulk_shear(self, low:int = -1, high:int = -1):
    h = self.height
    low = np.min(h) if low < 0 else low
    high = np.max(h) if high < 0 else high

    msk = (h >= low) & (h <= high)

    p = self.pressure[msk]
    u = self.wind_u[msk]
    v = self.wind_v[msk]

    msk = np.isnan(p) | np.isnan(u) | np.isnan(v)

    p = p[~msk] * units.mbar
    u = u[~msk] * units.meter / units.second
    v = v[~msk] * units.meter / units.second

    try:
      shear = bulk_shear(p, u, v)
      return Shear(shear[0], shear[1])
    except:
      return None
    
  def calculate_lfc(self)->LFC | None:

    lfc_data = mpcalc.lfc(self.pressure * units.mbar, self.temperature*units.degC, self.dewpoint*units.degC)

    # find closest pressure
    pressure = lfc_data[0]
    abs_vals = np.abs(self.pressure * units.mbar - pressure)
    if(np.all(np.isnan(abs_vals))):
      return None
    
    idx = np.nanargmin(abs_vals)
    lfc_height = self.height[idx]

    temperature = lfc_data[1]

    height = lfc_height
    return LFC(height, pressure, temperature)

def __load_sounding_df(file_name:str)->pd.DataFrame:
  clean_directory = "Full-CP20-sounding-dataset-clean"
  file_path = f'{clean_directory}/{file_name}'

  df = pd.read_csv(file_path)
  df = df.drop(df.columns[0], axis=1)
  return df

def __convert_df_to_metpy(df:pd.DataFrame)->pd.DataFrame:
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

def load_and_convert_sounding(file_path:str) -> Sounding:
  return Sounding(__convert_df_to_metpy(__load_sounding_df(file_path)))