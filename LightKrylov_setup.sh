#!/bin/bash

# Welcome message
echo "LightKrylov Setup Script"
echo "--------------------"
echo "This script will perform the following actions upon your confirmation:"
echo "1. Remove existing LightKrylov directory (if found)."
echo "3. Clone the LightKrylov repository."
echo "4. Build and export LightKrylov."
echo "--------------------"

# Detect the operating system
OS=$(uname)
should_clone="no"

# Install dependencies
read -p "Do you want to install/update necessary packages? (y/n): " confirm
if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
    echo "Installing dependencies..."
    if [ "$OS" == "Linux" ]; then
        conda config --add channels conda-forge
        conda install fpm
    elif [ "$OS" == "Darwin" ]; then
        brew tap fortran-lang/fortran
        brew update
        brew install fpm
    else
        echo "Unsupported operating system. Exiting."
        exit 1
    fi
else
    echo "Skipping dependencies installation."
fi

# Check for existing Nek5000 directory
if [ -d "LightKrylov" ]; then
    read -p "Directory LightKrylov exists. Do you want to remove it? (y/n): " confirm
    if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
        echo "Removing LightKrylov directory..."
        rm -rf LightKrylov
        should_clone="yes"
    else
        echo "LightKrylov directory will be retained."
    fi
else
    should_clone="yes"
fi

# Clone LightKrylov repository
if [ "$should_clone" == "yes" ]; then
    read -p "Do you want to clone the LightKrylov repository? (y/n): " confirm
    if [ "$confirm" == "y" ] || [ "$confirm" == "Y" ]; then
        echo "Cloning LightKrylov repository..."
        git clone --depth=1 https://github.com/nekStab/LightKrylov.git
        cd LightKrylov
        git checkout main
    else
        echo "Skipping cloning."
    fi
else
    cd LightKrylov
fi

fpm clean --all

# conda install -c conda-forge fortran-compiler openmpi libopenblas fpm

if command -v mpiifort >/dev/null 2>&1; then

    fpm install --compiler 'mpiifort' --flag "-Ofast -xHost -g -traceback"
    #fpm test --verbose  --compiler 'mpiifort' --flag "O0 -g -traceback -Wall"

else
    fpm install --compiler 'mpif90' --flag "-O3 -march=native -funroll-loops -ffast-math"
    fpm install --compiler 'mpif90' --flag "-O3 -march=native -funroll-loops -ffast-math -g -fbacktrace -fanalyzer -fbounds-check"


    fpm test --compiler 'mpif90' --flag "-O3 -march=native -funroll-loops -ffast-math -g -fbacktrace"
    fpm clean --all
    fpm install --compiler 'mpif90' --flag "-O3 -march=native -funroll-loops -ffast-math -g -fbacktrace -fanalyzer"
    
    #fpm test --verbose --flag "-Wall -Wextra -g -fbacktrace -fbounds-check"

fi

echo "LightKrylov setup complete."