from metpy_dataframe_loader import *
from calculations.sounding_shears import makeShear
import numpy as np

sounding = load_and_convert_sounding("snd-19940506-232513-storm1-10m.zagl")

print(np.isnan(sounding.wind_u))