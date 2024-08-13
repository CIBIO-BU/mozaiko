#!/bin/bash

ENV_NAME="dnaquaimg"
REPO_URL="git@github.com:CIBIO-BU/DNAquaIMG.git"
PACKAGE_DIR="DNAquaIMG"
CRABS_RELEASE="https://github.com/gjeunen/reference_database_creator/archive/refs/tags/v0.1.7.tar.gz"
CRABS_DIR="external_scripts/crabs"
CRABS_DIR=$(basename -s ".tar.gz" -a $CRABS_DIR)

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

# Install CRABS v0.1.7
install_crabs_release() {
    echo "This tool requires CRABS v0.1.7 for downstream analysis"
    echo "Checking if CRABS v0.1.7 is installed"

    crabs_output=$(crabs --version | tail -n 1)
    crabs_version=${crabs_output##* }

    if [ "$crabs_version" != '0.1.7' ]; then
        if [ -n "$crabs_version" ]; then
            echo "CRABS is installed with the wrong version. Please remove current CRABS installation an install 0.1.7"
        else
            echo "CRABS is not installed"
            echo "Downloading CRABS v0.1.7 from $CRABS_RELEASE"
            wget "$CRABS_RELEASE" -O "${CRABS_DIR}.tar.gz" || { echo "Failed to download CRABS"; exit 1; }

            echo "Unzipping CRABS"
            tar -xzf "${CRABS_DIR}.tar.gz" -C "$(dirname "$CRABS_DIR")" || { echo "Failed to unzip CRABS"; exit 1; }

            echo "Navigating to CRABS directory: $CRABS_DIR"
            cd "$CRABS_DIR" || { echo "Directory $CRABS_DIR does not exist"; exit 1; }

            echo "Installing CRABS"
            pip install . || { echo "Failed to install CRABS"; exit 1; }

            echo "CRABS Installation complete."
        fi
    else
        echo echo "CRABS is installed with 0.1.7 version."
    fi
}

main() {
    check_conda
    check_env
    activate_env
    clone_repo
    install_package
    install_crabs_release
}

main