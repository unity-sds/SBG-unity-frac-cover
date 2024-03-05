#!/bin/bash

# Install script to 1) install Julia via conda, and 2) install the Julia dependencies for this project

set -x

# Define some useful directories
pge_dir=$( cd "$(dirname "$0")" ; pwd -P )
app_dir=$(dirname ${pge_dir})

cd $app_dir

# Get SpectralUnmixing repo
git clone https://github.com/EnSpec/SpectralUnmixing.git -b v0.2.1-sister
specun_dir="$app_dir/SpectralUnmixing"

# Install Julia and then install Julia dependencies
conda create -n spectral-unmixing -y -c conda-forge julia=1.7 python=3.8 gdal=3.6.0 pandas=2.0.1 awscli
source activate spectral-unmixing
python -m pip install awscli
pip install Pillow
pip install pystac==1.8.4
pip install spectral==0.23.1

#needed for unity-sds
pip install poetry
git clone https://github.com/unity-sds/unity-py.git
cd unity-py
git switch develop
poetry install 

# Download snow climatology dataset
#aws s3 cp s3://sister-ops-registry/packages/LIN10A1_snow_climatology_13day.tif .

git clone https://github.com/EnSpec/hytools.git
cd hytools
pip install .
cd ..

echo 'Spectral Unmixing repo is here: ' $specun_dir
pushd $specun_dir
export JULIA_SSL_CA_ROOTS_PATH=""
julia -e 'using Pkg; Pkg.activate("."); Pkg.add(path="https://github.com/kmsquire/ArgParse2.jl"); Pkg.instantiate()'
julia -e 'using Pkg; Pkg.activate("."); Pkg.precompile()'
export JULIA_PROJECT=`pwd`
echo $JULIA_PROJECT
popd
cd $app_dir

