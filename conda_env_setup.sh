#!/bin/bash

ENV_NAME="dnaquaimg"
REPO_URL="https://github.com/CIBIO-BU/DNAquaIMG"
PACKAGE_DIR="DNAquaIMG"

# Check if Conda is installed
check_conda() {
    if ! command -v conda &> /dev/null; then
        echo "Conda could not be found. Please install Conda first."
        exit 1
    fi
}

# Check if env already exists
check_env() {
    if conda env list | grep ".*${ENV_NAME}.*" >/dev/null 2>&1; then
        echo "Conda environment '$ENV_NAME' already exists. Skipping environment creation."
    else
        echo "Creating Conda environment: $ENV_NAME"
        conda env create -f environment.yml
    fi
}

# Activate Environment
activate_env() {
    echo "Activating Conda environment: $ENV_NAME"
    source /opt/conda/etc/profile.d/conda.sh
    conda activate $ENV_NAME
}

# Clone repository
clone_repo() {
    if [ -d "$PACKAGE_DIR" ]; then
        echo "Directory $PACKAGE_DIR already exists. Skipping repository cloning."
    else
        echo "Cloning repository from $REPO_URL."
        git clone $REPO_URL
    fi
}

# Install Package
install_package() {
    echo "Navigating to package directory: $PACKAGE_DIR"
    cd $PACKAGE_DIR || exit

    echo "Installing package"
    pip install .

    echo "Installation complete."
}

main() {
    check_conda
    check_env
    activate_env
    clone_repo
    install_package
}

main