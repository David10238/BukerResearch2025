from metpy_dataframe_loader import Sounding

from metpy.calc import cape_cin, parcel_profile
from metpy.units import units
import numpy as np

class StormEnergy:
  def __init__(self, cape:float, cin:float):
    self.cape = cape
    self.cin = cin

def makeStormEnergy(sounding:Sounding, maxHeight:int = -1):
  if(maxHeight < 0):
    maxHeight = sounding.height.max()    

  H = sounding.height
  msk = H <= maxHeight

  p = sounding.pressure
  p = p[msk]

  T = sounding.temperature
  T = T[msk]

  Td = sounding.dewpoint
  Td = Td[msk]

  # make sure data exists
  msk = np.isnan(p) | np.isnan(T) | np.isnan(Td)

  if np.all(msk):
    return None
  
  p = p[~msk] * units.mbar
  T = T[~msk] * units.degC
  Td = Td[~msk] * units.degC

  # ensure pressure is decreasing
  cleaned_pressure = p[~np.isnan(p)]
  if np.any(np.diff(cleaned_pressure) > 0):
    return None
  
  prof = parcel_profile(p, T[0], Td[0]).to('degC')
  energy = cape_cin(p, T, Td, prof)

  return StormEnergy(energy[0], energy[1])