
import numpy as np
import pandas as pd
from metpy.units import units

distance = np.arange(1, 5) * units.meters

time = units.Quantity(np.arange(2, 10, 2), 'sec')

print(distance / time)

path = "Full-CP20-sounding-dataset-clean/snd-19940409-205903-storm1-10m.zagl"

df = pd.read_csv(path)
df = df.drop(df.columns[0], axis=1)

print(df.head(10))