#!/bin/bash

# Path to your exported file
FILE_PATH="renv_exported_packages.txt"

# Activate conda environment
conda activate my_r_env

# Loop through the file and install packages via conda
while IFS=, read -r package version
do
    echo "Installing $package version $version..."
    
    # Attempt to install the R package via conda
    conda install "r-$package=$version" --yes
    
    # If the package isn't available via conda, use R to install it
    if [ $? -ne 0 ]; then
        echo "$package not found in conda, trying CRAN..."
        Rscript -e "install.packages('$package', version = '$version', repos = 'https://cloud.r-project.org/')"
    fi

done < "$FILE_PATH"
