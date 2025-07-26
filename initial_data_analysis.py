import os

from sounding import load_and_convert_sounding, StormMetadataLoader
import numpy as np
import pandas as pd

loader = StormMetadataLoader()

storm_speed = []
srh = []
srh_left = []
srh_right = []
bulk_shear = []
tornado = []
ef_rating = []
cape = []
cin = []


for file in os.listdir("Full-CP20-sounding-dataset-clean"):
  sounding = load_and_convert_sounding(loader, file)

  storm_speed.append(sounding.metadata.storm_speed)

  srh.append(sounding.calculate_superhelicity())
  srh_left.append(sounding.calculate_left_srh())
  srh_right.append(sounding.calculate_right_srh())

  _bulk_shear = sounding.calculate_bulk_shear()
  bulk_shear.append((_bulk_shear.u**2 + _bulk_shear.v**2)**0.5 if _bulk_shear else None)

  tornado.append(1 if sounding.metadata.tornado else 0)
  ef_rating.append(sounding.metadata.ef_rating)

  _storm_energy = sounding.calculate_energy()
  cape.append(_storm_energy.cape if _storm_energy else None)
  cin.append(_storm_energy.cin if _storm_energy else None)

df = pd.DataFrame({
  "storm_speed": storm_speed,
  "srh": srh,
  "srh_left": srh_left,
  "srh_right": srh_right,
  "bulk_shear": bulk_shear,
  "tornado": tornado,
  "ef_rating": ef_rating,
  "cape": cape,
  "cin": cin,
})

df.to_csv("initial_data.csv")