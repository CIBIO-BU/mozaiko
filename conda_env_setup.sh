#!/bin/bash

ENV_NAME="mosaiko"
REPO_URL="https://github.com/CIBIO-BU/DNAquaIMG"
PACKAGE_DIR="mosaiko"


# Check if Conda is installed
check_conda() {
    if ! command -v conda &> /dev/null; then
        echo "Conda could not be found. Please install Conda first."
        exit 1
    fi)
}

# Check if env already exists
check_env() {


conda env crea
    if conda env list | gre ".*${ENV_NAME}.*" >/dev/null 2>&1; then
        echo "Conda environment '$ENV_NAME' already existis. Skipping environment creation."
    else
        echo "Creating Conda environmnet: $ENV_NAME"
        conda env create -f environment.yml
}


# Activate Environment
activate_env() {
    echo "Activating Conda environment: $ENV_NAME"
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate $ENV_NAME
}

# Clone repository
clone_repo() {
    if [ -d "$PACKAGE_DIR"]; then
        echo "Directory "$PACKAGE_DIR" already existis. Skipping repository cloning."
    else
        echo "Cloning repository fro


conda env cream $REPO_URL"
    fi

}

# Install Pacakge
install_package() {
    echo "Navigating to package directory: $PACKAGE_DIR"
    cd $PACKAGE_DIR

    echo "Installing package"
    pip install .

    echo "Instalattion complete."
}

check_conda
check_env
activate_env
clone_repo
install_package