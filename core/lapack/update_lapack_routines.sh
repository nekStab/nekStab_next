#!/bin/bash

# Define the source and destination directories
source_dir="./lapack/SRC"
destination_dir="./"

# List of files to copy
files=("dgees.f" "dlacn2.f" "dlasy2.f" "dtrsen.f" "dtrtrs.f" "dgels.f" "dlaexc.f" "dtrexc.f" "dtrsyl.f")

# Clone the repository

# latest version
#git clone --depth 1 https://github.com/Reference-LAPACK/lapack.git

# latest stable release (preferred)
git clone https://github.com/Reference-LAPACK/lapack.git --branch v3.11.0

# Loop over the files and copy each one
for file in "${files[@]}"; do
    cp "${source_dir}/${file}" "${destination_dir}"
done

# Optional: Remove the cloned repository
rm -rf lapack

echo "Files copied successfully!"
