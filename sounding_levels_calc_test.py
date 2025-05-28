
import pandas as pd
import metpy.calc as mpcalc
import metpy.units as units
import os

from metpy_dataframe_loader import *

name = "snd-19940507-000304-storm2-10m.zagl"
sounding = load_and_convert_sounding(name)

#print(sounding.as_df.head())

print("LCF Height:", sounding.level_free_convection.height)
print("LCF Pressure", sounding.level_free_convection.pressure)
print("LCF Temperature", sounding.level_free_convection.temperature)

# verify my assertion works
clean_directory = "Full-CP20-sounding-dataset-clean"
# process all files
all_files = os.listdir(clean_directory)

found_lfc = 0
missing_lfc = 0
for file_name in all_files:
  print("Found ", found_lfc)
  print("Missing ", missing_lfc)
  print(file_name)
  sounding = load_and_convert_sounding(file_name)
  lfc = sounding.level_free_convection
  if not lfc:
    missing_lfc = missing_lfc + 1
  else:
    found_lfc = found_lfc + 1