#!/bin/bash

ENV_NAME="dnaquaimg"
REPO_URL="git@github.com:CIBIO-BU/DNAquaIMG.git"
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
    if conda env list | grep -q "${ENV_NAME}"; then
        echo "Conda environment '$ENV_NAME' already exists. Skipping environment creation."
    else
        echo "Creating Conda environment: $ENV_NAME"
        conda env create -f environment.yml || { echo "Failed to create Conda environment"; exit 1; }
    fi
}

# Activate Environment
activate_env() {
    echo "Activating Conda environment: $ENV_NAME"

    # Ensure Conda is initialized
    # This step might be needed to properly configure the shell for Conda
    CUR_SHELL=shell.$(basename -- "${SHELL}")
    eval "$(conda "$CUR_SHELL" hook)"

    set -e
    conda activate "$ENV"
    echo "INFO: conda environment $ENV created and activated"
}

# Clone repository
clone_repo() {
    if [ -d "$PACKAGE_DIR" ]; then
        echo "Directory $PACKAGE_DIR already exists. Skipping repository cloning."
    else
        echo "Cloning repository from $REPO_URL."
        git clone "$REPO_URL" || { echo "Failed to clone repository"; exit 1; }
    fi
}

# Install Package
install_package() {
    echo "Navigating to package directory: $PACKAGE_DIR"
    cd "$PACKAGE_DIR" || { echo "Directory $PACKAGE_DIR does not exist"; exit 1; }

    echo "Installing package"
    pip install . || { echo "Failed to install package"; exit 1; }

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