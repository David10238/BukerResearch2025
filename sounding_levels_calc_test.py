import os

from sounding import load_and_convert_sounding

name = "snd-19940507-000304-storm2-10m.zagl"
sounding = load_and_convert_sounding(name)

#print(sounding.as_df.head())

lfc = sounding.calculate_lfc()

if lfc is not None:
  print("LCF Height:", lfc.height)
  print("LCF Pressure", lfc.pressure)
  print("LCF Temperature", lfc.temperature)

# verify my assertion works
clean_directory = "Full-CP20-sounding-dataset-clean"
# process all files
all_files = os.listdir(clean_directory)

found_lfc = 0
missing_lfc = 0
for file_name in all_files:
  sounding = load_and_convert_sounding(file_name)
  lfc = sounding.calculate_lfc
  if not lfc:
    missing_lfc = missing_lfc + 1
  else:
    found_lfc = found_lfc + 1

print("Found ", found_lfc)
print("Missing ", missing_lfc)