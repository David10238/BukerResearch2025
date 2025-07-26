import os

from sounding import load_and_convert_sounding, StormMetadataLoader
import numpy as np

loader = StormMetadataLoader()

name = "snd-19940507-000304-storm2-10m.zagl"
sounding = load_and_convert_sounding(loader, name)

storm_motion = sounding.calculate_bunkers_storm_motion()
super_helicity = sounding.calculate_superhelicity()
srh_left = sounding.calculate_left_srh()
srh_right = sounding.calculate_right_srh()
print(storm_motion)
print(super_helicity)
print(srh_left)
print(srh_right)


clean_directory = "Full-CP20-sounding-dataset-clean"
all_files = os.listdir(clean_directory)

found_motions = 0
missed_motions = 0

for file_name in all_files:
  sounding = load_and_convert_sounding(loader, file_name)
  storm_motion = sounding.calculate_bunkers_storm_motion()
  super_helicity = sounding.calculate_superhelicity()
  srh_left = sounding.calculate_left_srh()
  srh_right = sounding.calculate_right_srh()

  if not storm_motion or not srh_left or not srh_right or not super_helicity or np.isnan(srh_left) or np.isnan(srh_right) or np.isnan(super_helicity):
    missed_motions = missed_motions + 1
  else:
    found_motions = found_motions + 1

print("done without errors")
print("found", found_motions)
print("missing", missed_motions)