#!/bin/bash
set -e

ENV_NAME=IsoLamp

echo "Creating conda environment..."
conda env create -f IsoLamp_env.yml

echo "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

echo "Installing oarfish via cargo..."
cargo install --root $CONDA_PREFIX oarfish

echo "Installation complete."
