import pandas as pd
import numpy as np
from metpy.units import units
import metpy.calc as mpcalc
from metpy.calc import cape_cin, parcel_profile, bulk_shear, lifted_index, mixed_parcel
from numpy import concatenate

from typing import TypedDict
import pandas as pd

class StormMetadata:
  def __init__(self, latitude:float, longitude:float, storm_direction:float, storm_speed:float, tornado:bool, ef_rating:int | None):
    self.latitude = latitude
    self.longitude = longitude
    self.storm_direction = storm_direction
    self.storm_speed = storm_speed
    self.tornado = tornado
    self.ef_rating = ef_rating
    pass

class StormMetadataLoader:
  def __init__(self):
    df = pd.read_csv("snd-storm-attributes.csv")

    date = df['yyyymmdd'].tolist()
    time = df['hhmmss'].tolist()
    storm_number = df['storm#'].tolist()
    latitude = df["latstorm"].tolist()
    longitude = df["lonstorm"].tolist()
    storm_direction = df['dirstm'].tolist()
    storm_speed = df['spdstm'].tolist()
    tornado = df['tor'].tolist()

    self._metadata = dict()

    for index in range(len(df)):
      file_name = f'snd-{date[index]}-{str(time[index]).rjust(6, "0")}-{storm_number[index]}-10m.zagl'
      
      i_tornado = tornado[index]

      self._metadata[file_name] = StormMetadata(
        latitude=latitude[index],
        longitude=longitude[index],
        storm_direction=storm_direction[index],
        storm_speed=storm_speed[index],
        tornado = i_tornado >= 0,
        ef_rating = i_tornado if (i_tornado >= 0 and i_tornado <= 5) else None
      )

  def of_storm(self, fileName:str) -> StormMetadata:
    return self._metadata[fileName]

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

class BunkersMotion:
  def __init__(self, leftMotion:tuple[float, float], rightMotion:tuple[float, float], meanWind:tuple[float, float]):
    self.leftMotion = leftMotion
    self.rightMotion = rightMotion
    self.meanWind = meanWind

# time (second)
# pressure (mbar)
# height (meter)
# temperature (degC)
# dewpoint (degC)
# speed (meter)
class Sounding:
  def __init__(self, metadata:StormMetadata, df:pd.DataFrame):
    self.metadata = metadata
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
    
  def calculate_lifted_index(self)->float | None:
    p = self.pressure * units.mbar
    T = self.temperature * units.degC
    Td = self.dewpoint * units.degC
    h = self.height * units.m
    
    try:
      # Calculate 500m mixed parcel
      parcel_p, parcel_t, parcel_td = mixed_parcel(p, T, Td, depth=500 * units.m, height=h)

      # Replace sounding temp/pressure in lowest 500m with mixed values
      above = h > 500 * units.m
      press = concatenate([[parcel_p], p[above]])
      temp = concatenate([[parcel_t], T[above]])

      # Calculate parcel profile from our new mixed parcel
      mixed_prof = parcel_profile(press, parcel_t, parcel_td)

      # Calculate lifted index using our mixed profile
      return lifted_index(press, temp, mixed_prof)[0]
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
  
  def calculate_bunkers_storm_motion(self)->BunkersMotion | None:
    index_surface = 0
    index_5km = (5000 // 10) - 1

    surface_wind = (self.wind_u[index_surface], self.wind_v[index_surface])
    wind_5km = (self.wind_u[index_5km], self.wind_v[index_5km])

    if np.isnan(surface_wind[0]) | np.isnan(surface_wind[1]) | np.isnan(wind_5km[0]) | np.isnan(wind_5km[1]):
      return None

    # Mean wind
    mean_wind = ((surface_wind[0] + wind_5km[0]) / 2, (surface_wind[1] + wind_5km[1]) / 2)
    
    # Shear vector
    shear_vector = (wind_5km[0] - surface_wind[0], wind_5km[1] - surface_wind[1])
    
    # Perpendicular shear vector (rotate 90 degrees to the right)
    perp_shear_vector = (-shear_vector[1], shear_vector[0])
    
    # Scale the perpendicular shear vector (7.5 m/s is typical)
    scale_factor = 7.5 / np.sqrt(perp_shear_vector[0]**2 + perp_shear_vector[1]**2)
    scaled_perp_shear = (perp_shear_vector[0] * scale_factor, perp_shear_vector[1] * scale_factor)
    
    # Right-moving storm motion
    right_motion = (mean_wind[0] - scaled_perp_shear[0], mean_wind[1] - scaled_perp_shear[1])
    
    # Left-moving storm motion
    left_motion = (mean_wind[0] + scaled_perp_shear[0], mean_wind[1] + scaled_perp_shear[1])

    return BunkersMotion(left_motion, right_motion, mean_wind)
  
  def calculate_superhelicity(self, depth:int=1000)-> float | None:
    mask = self.height > depth | np.isnan(self.height) | np.isnan(self.wind_u) | np.isnan(self.wind_v)

    heights = self.height[~mask]
    u = self.wind_u[~mask]
    v = self.wind_v[~mask]

    if len(heights) < 10:
      return None

    # Calculate vertical shear (first derivatives)
    du_dz = np.diff(u) / np.diff(heights)
    dv_dz = np.diff(v) / np.diff(heights)
    
    # Calculate shear of shear (second derivatives)
    d2u_dz2 = np.diff(du_dz) / np.diff(heights[:-1])
    d2v_dz2 = np.diff(dv_dz) / np.diff(heights[:-1])
    
    # Calculate superhelicity
    superhelicity = 0
    for i in range(len(d2u_dz2)):
        term1 = dv_dz[i] * d2u_dz2[i]
        term2 = du_dz[i] * d2v_dz2[i]
        superhelicity += (term1 + term2) * (heights[i + 2] - heights[i + 1])  # Approximate integral
    
    return superhelicity
  
  def calculate_srh(self, storm_motion:tuple[float, float], depth:int=0)->float | None:
    mask = self.height > depth | np.isnan(self.height) | np.isnan(self.wind_u) | np.isnan(self.wind_v)

    heights = self.height[~mask]
    u = self.wind_u[~mask]
    v = self.wind_v[~mask]

    if len(heights) < 10:
      return None

    # Calculate storm-relative wind
    u_relative = u - storm_motion[0]
    v_relative = v - storm_motion[1]
    
    # Calculate vertical shear
    du_dz = np.diff(u) / np.diff(heights)
    dv_dz = np.diff(v) / np.diff(heights)
    
    # Calculate horizontal vorticity
    vorticity_x = -dv_dz  # ∂v/∂z
    vorticity_y = du_dz  # -∂u/∂z
    
    # Calculate the dot product of storm-relative wind and horizontal vorticity
    srh = 0
    for i in range(len(du_dz)):
        dot_product = (u_relative[i] * vorticity_x[i] + v_relative[i] * vorticity_y[i])
        srh += dot_product * (heights[i + 1] - heights[i])  # Approximate integral
    
    return srh
  
  def calculate_left_srh(self, depth:int=0)->float | None:
    storm_motion = self.calculate_bunkers_storm_motion()
    if not storm_motion:
      return None
    return self.calculate_srh(storm_motion.leftMotion, depth)
  
  def calculate_right_srh(self, depth:int=0)->float | None:
    storm_motion = self.calculate_bunkers_storm_motion()
    if not storm_motion:
      return None
    return self.calculate_srh(storm_motion.rightMotion, depth)


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

def load_and_convert_sounding(metadataLoader:StormMetadataLoader, file_path:str) -> Sounding:
  return Sounding(metadataLoader.of_storm(file_path.split('/')[-1]), __convert_df_to_metpy(__load_sounding_df(file_path)))