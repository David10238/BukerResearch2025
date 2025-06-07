
import pandas as pd
import metpy.calc as mpcalc
import metpy.units as units
import os

from metpy_dataframe_loader import *
from calculations.sounding_levels import *

name = "snd-19940507-000304-storm2-10m.zagl"
sounding = load_and_convert_sounding(name)

#print(sounding.as_df.head())

lfc = LFC.makeLFC(sounding)
lcl = LCL.makeLCL(sounding)

print("LCF Height:", lfc.height)
print("LCF Pressure", lfc.pressure)
print("LCF Temperature", lfc.temperature)

print("LCL Height", lcl.height)

# verify my assertion works
clean_directory = "Full-CP20-sounding-dataset-clean"
# process all files
all_files = os.listdir(clean_directory)

found_lfc = 0
missing_lfc = 0
for file_name in all_files:
  sounding = load_and_convert_sounding(file_name)
  lfc = LFC.makeLFC(sounding)
  if not lfc:
    missing_lfc = missing_lfc + 1
  else:
    found_lfc = found_lfc + 1

print("Found ", found_lfc)
print("Missing ", missing_lfc)