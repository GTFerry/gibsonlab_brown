#!/bin/bash

# Number of lines per file
n=$1
if [ -z "$n" ]; then
  echo "Usage: $0 <number of lines per file>"
  exit 1
fi

# Create a temporary root directory
root_dir="temp_root"
mkdir -p $root_dir

# Exclude the header and count the number of lines
total_lines=$(tail -n +2 input.csv | wc -l)

# Number of files needed
num_files=$(( (total_lines + n - 1) / n ))

# Generate the split CSV files
awk -v n="$n" -v header="$(head -n 1 input.csv)" 'NR==1 {next} {filename=sprintf("'${root_dir}'/input_part%02d.csv", int((NR-2)/n)+1); if (NR % n == 2) print header > filename; print >> filename}' input.csv

# Iterate over each part, create a directory, copy R script, and link to shared renv
for i in $(seq 1 $num_files); do
  temp_dir="${root_dir}/dir_${i}"
  mkdir -p $temp_dir
  cp "${root_dir}/input_part$(printf "%02d" $i).csv" "${temp_dir}/input.csv"
  cp "clustering-script.R" "${temp_dir}/clustering-script.R"

  # Create a symbolic link to the shared renv directory
  ln -s "$(pwd)/renv" "${temp_dir}/renv"

  # Run the R script in background
  (cd $temp_dir && Rscript clustering-script.R) &
done

# Wait for all background processes to complete
wait

# Optionally: Remove the temporary root directory after completion
# rm -rf $root_dir
