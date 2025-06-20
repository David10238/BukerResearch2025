
import pandas as pd
import numpy as np
from metpy.units import units
import metpy.calc as mpcalc

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