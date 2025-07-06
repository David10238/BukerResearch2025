import os
from sounding import load_and_convert_sounding

name = "snd-19940507-000304-storm2-10m.zagl"
sounding = load_and_convert_sounding(name)

shear = sounding.calculate_bulk_shear()
if shear is not None:
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
  shear = sounding.calculate_bulk_shear()
  if not shear:
    missing_shear = missing_shear + 1
  else:
    found_shear = found_shear + 1

print("Found ", found_shear)
print("Missing ", missing_shear)