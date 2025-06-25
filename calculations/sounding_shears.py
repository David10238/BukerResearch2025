from metpy_dataframe_loader import Sounding

from metpy.calc import bulk_shear
from metpy.units import units
import numpy as np

class Shear:
  def __init__(self, u:float, v:float):
    self.u = u
    self.v = v

def makeShear(sounding: Sounding, low:int = -1, high:int = -1):
  h = sounding.height
  if low < 0:
    low = np.min(sounding.height)
  if high < 0:
    high = np.max(sounding.height)

  msk = (h >= low) & (h <= high)

  p = sounding.pressure[msk]
  u = sounding.wind_u[msk]
  v = sounding.wind_v[msk]

  msk = np.isnan(p) | np.isnan(u) | np.isnan(v)

  if np.all(msk):
    return None
  
  p = p[~msk]
  u = u[~msk]
  v = v[~msk]

  p = p * units.mbar
  u = u * units.meter / units.second
  v = v * units.meter / units.second

  try:
    shear = bulk_shear(p, u, v)
    return Shear(shear[0], shear[1])
  except:
    print("P: ", p)
    print("U:", u)
    print("V", v)
    print()
    return None