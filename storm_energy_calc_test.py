from sounding import load_and_convert_sounding
import os

name = "snd-20090609-235406-storm1-10m.zagl"
sounding = load_and_convert_sounding(name)

storm_energy = sounding.calculate_energy()

if storm_energy is not None:
  print("Cape:", storm_energy.cape)
  print("Cin:", storm_energy.cin)

# process all files
clean_directory = "Full-CP20-sounding-dataset-clean"
all_files = os.listdir(clean_directory)

found_energy = 0
missing_energy = 0
for file_name in all_files:
  sounding = load_and_convert_sounding(file_name)
  energy = sounding.calculate_energy()
  if not energy:
    missing_energy = missing_energy + 1
  else:
    found_energy = found_energy + 1

print("Found: ", found_energy)
print("Missing: ", missing_energy)