#!/bin/bash

ENV_NAME="dnaquaimg"
REPO_URL="git@github.com:CIBIO-BU/DNAquaIMG.git"
PACKAGE_DIR="DNAquaIMG"
CRABS_RELEASE="https://github.com/gjeunen/reference_database_creator/archive/refs/tags/v0.1.7.tar.gz"
EXTERNAL_SCRIPTS_DIR="external_scripts"
CRABS_ARCHIVE="crabs.tar.gz"
CRABS_DIR="reference_database_creator-0.1.7"

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
    conda activate "$ENV_NAME"
    echo "INFO: conda environment $ENV_NAME created and activated"
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
    echo "mosaiko requires CRABS v0.1.7 for downstream analysis"
    echo "Checking if CRABS v0.1.7 is installed"

    crabs_output=$(crabs --version | tail -n 1)
    crabs_version=${crabs_output##* }

    if [ "$crabs_version" != '0.1.7' ]; then
        if [ -n "$crabs_version" ]; then
            echo "CRABS is installed with the wrong version. Please remove current CRABS installation an install 0.1.7"
        else
            echo "CRABS is not installed. Installing CRABS v0.1.7..."

            echo "Creating $EXTERNAL_SCRIPTS_DIR if it doesn't exist..."
            mkdir -p "$EXTERNAL_SCRIPTS_DIR" || { echo "Failed to create directory $EXTERNAL_SCRIPTS_DIR"; exit 1; }

            echo "Moving to $EXTERNAL_SCRIPTS_DIR..."
            cd "$EXTERNAL_SCRIPTS_DIR" || { echo "Directory $EXTERNAL_SCRIPTS_DIR does not exist"; exit 1; }

            echo "Downloading CRABS v0.1.7"
            wget "$CRABS_RELEASE" -O "$CRABS_ARCHIVE" || { echo "Failed to download CRABS"; exit 1; }

            echo "Unzipping CRABS"
            tar -xzf "$CRABS_ARCHIVE" || { echo "Failed to unzip CRABS"; exit 1; }

            echo "Navigating to CRABS directory"
            cd "$CRABS_DIR" || { echo "Directory reference_database_creator-0.1.7 does not exist"; exit 1; }

            echo "Installing CRABS"
            pip install . || { echo "Failed to install CRABS"; exit 1; }

            echo "Deleting CRABS archive"
            cd ..
            rm -f "$CRABS_ARCHIVE"

            echo "CRABS Installation complete."
        fi
    else
        echo "CRABS is installed with version 0.1.7."
    fi
}

install_cutadapt_package() {
    echo "mosaiko requires cutadapt for downstream analysis"

    echo "checking if cutadapt is already installed"

    if command -v cutadapt &> /dev/null; then
        current_version=$(cutadapt --version | cut -d ' ' -f2)
        required_version="4.9"

        if [ "$(printf '%s\n' "$required_version" "$current_version" | sort -V | head -n1)" = "$required_version" ]; then
            echo "Cutadapt version $current_version is already installed and meets the minimum requirement."
            return 0
        else
            echo "Cutadapt is installed but version $current_version is outdated. Minimum required version is $required_version."
        fi
    else
        echo "Cutadapt is not installed."
    fi

    echo "Installing/Upgrading cutadapt package..."

    if ! apt-get install -y pipx python3-venv; then
        echo "Failed to install pipx and python3-venv. Please check your system and try again."
        return 1
    fi

    if ! pipx install --force cutadapt; then
        echo "Failed to install cutadapt. Please check the error messages and try again."
        return 1
    fi

    if command -v cutadapt &> /dev/null; then
        echo "Cutadapt has been successfully installed."
        cutadapt --version
        return 0
    else
        echo "Cutadapt installation failed. Please check the error messages and try again."
        return 1
    fi


}

main() {
    check_conda
    check_env
    activate_env
    clone_repo
    install_package
    install_crabs_release
    install_cutadapt_package
}

main