import os

from metpy_dataframe_loader import *
from calculations.sounding_shears import makeShear

name = "snd-19940507-000304-storm2-10m.zagl"
sounding = load_and_convert_sounding(name)

#print(sounding.as_df.head())

shear = makeShear(sounding)

print("Shear U", shear.u)
print("Shear V", shear.v)

# verify my assertion works
clean_directory = "Full-CP20-sounding-dataset-clean"
# process all files
all_files = os.listdir(clean_directory)

found_shear = 0
missing_shear = 0
for file_name in all_files: # snd-19940506-232513-storm1-10m. is the error
  #print(file_name)
  sounding = load_and_convert_sounding(file_name)
  lfc = makeShear(sounding)
  if not lfc:
    missing_shear = missing_shear + 1
  else:
    found_shear = found_shear + 1

print("Found ", found_shear)
print("Missing ", missing_shear)