
import numpy as np
import pandas as pd
from metpy.units import units
from metpy_dataframe_loader import *

distance = np.arange(1, 5) * units.meters

time = units.Quantity(np.arange(2, 10, 2), 'sec')

print(distance / time)

name = "snd-19940507-000304-storm2-10m.zagl"

raw_df = load_sounding_df(name)
metpy_df = convert_df_to_metpy(raw_df)

print(raw_df.head(10))
print(metpy_df.head(10))