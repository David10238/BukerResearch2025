import pandas as pd
import numpy as np
import re
import os

raw_directory = "Full-CP20-sounding-dataset-raw"
clean_directory = "Full-CP20-sounding-dataset-clean"

def process_file(name:str):
  input_path = f'{raw_directory}/{name}'
  output_path = f'{clean_directory}/{name}'

  input_file = open(input_path, "r")
  output_file = open(output_path, 'w+')

  # remove header to file
  while not input_file.readline().strip().startswith('---'):
    pass

  # setup column names
  column_names = ['time', 'pressure', 'temperature', 'dewpoint', 'wind_u', 'wind_v', 'longitude', 'latitude', 'height']
  output_file.write(','.join(column_names) + '\n')

  # strip whitespace
  for line in input_file:
    csv_line = re.sub('\s+', ',', line.strip())
    output_file.write(csv_line + '\n')
  
  input_file.close()
  output_file.close()

  # remove missing values
  df = pd.read_csv(output_path)
  df = df.replace(-9999.0, np.nan)
  df.to_csv(output_path)

# process all files
all_files = os.listdir(raw_directory)
for file_name in all_files:
  process_file(file_name)