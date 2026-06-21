#!/bin/bash
set -e

ENV_NAME=IsoLamp

echo "Creating conda environment..."
conda env create -f IsoLamp_env.yml

echo "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate $ENV_NAME

echo "Installing oarfish..."
if cargo install --root "$CONDA_PREFIX" oarfish; then
    echo "Installed latest oarfish."
else
    echo "Latest version unavailable, installing oarfish 0.9.0."
    cargo install --root "$CONDA_PREFIX" oarfish --version 0.9.0
fi

if [ -x "$CONDA_PREFIX/bin/oarfish" ]; then
    echo "oarfish installed successfully:"
    "$CONDA_PREFIX/bin/oarfish" --version
else
    echo "ERROR: oarfish installation failed."
    exit 1
fi

echo "Installation complete."
