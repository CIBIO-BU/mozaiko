#!/bin/bash

ENV_NAME="mozaiko"
REPO_URL="git@github.com:CIBIO-BU/mozaiko.git"
PACKAGE_DIR="mozaiko"
CRABS_RELEASE="https://github.com/gjeunen/reference_database_creator/archive/refs/tags/v0.1.7.tar.gz"
EXTERNAL_SCRIPTS_DIR="external_scripts"
CRABS_ARCHIVE="crabs.tar.gz"
CRABS_DIR="reference_database_creator-0.1.7"

# If a token argument is given, modify REPO_URL to use HTTPS with the token
if [ $# -gt 0 ]; then
    TOKEN="$1"
    REPO_URL="https://${TOKEN}@github.com/CIBIO-BU/mozaiko.git"
    echo "Using token-based repository URL."
else
    echo "Using default SSH repository URL."
fi

# Continue with the rest of your script
echo "Repository URL: $REPO_URL"

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
    if ! command -v git &> /dev/null; then
        echo "Git is not installed. Attempting to install git..."
        sudo apt-get update
        sudo apt-get install -y git
    fi

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

    echo "Installing mozaiko package"
    pip install . || { echo "Failed to install package"; exit 1; }

    echo "Installation complete."
}

# Install CRABS v0.1.7
install_crabs_release() {
    #echo "mozaiko requires CRABS v0.1.7 for downstream analysis"
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

            echo "Applying coverage calculation patch..."
            CRABS_PATH=$(which crabs)
            if [ -z "$CRABS_PATH" ]; then
                echo "Error: Could not locate crabs executable after installation"
                exit 1
            fi

            # Create backup
            cp "$CRABS_PATH" "${CRABS_PATH}.backup.$(date +%Y%m%d_%H%M%S)"

            # Apply the patch
            sed -i.bak '
                /print(f'\''filtering alignments based on parameter settings'\'')/ a\    whole_percent = float(COV)*100
                s/tcov >= float(COV) and/tcov >= whole_percent and/
                s/tcov >= float(COV)"/tcov >= whole_percent"/
            ' "$CRABS_PATH"

            # Restore executable permissions
            chmod +x "$CRABS_PATH"

            echo "CRABS patch applied successfully"

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
    #echo "mozaiko requires cutadapt for downstream analysis"
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

    if ! pip install cutadapt; then
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

install_vsearch() {
    #echo "mozaiko requires vsearch for downstream analysis"
    echo "checking if vsearch is already installed"

    if command -v vsearch &> /dev/null; then
        current_version=$(vsearch --version | cut -d ' ' -f2)
        required_version="2.13.3"

        if [ "$(printf '%s\n' "$required_version" "$current_version" | sort -V | head -n1)" = "$required_version" ]; then
            echo "Vsearch version $current_version is already installed and meets the minimum requirement."
            return 0
        else
            echo "Vsearch is installed but version $current_version is outdated. Minimum required version is $required_version."
        fi
    else
        echo "Vsearch is not installed."
    fi

    echo "Installing vsearch version 2.13.3 from source..."

    arch=$(uname -m)
    os=$(uname -s)

    case "${os}_${arch}" in
        "Linux_x86_64")
            url="https://github.com/torognes/vsearch/releases/download/v2.13.3/vsearch-2.13.3-linux-x86_64.tar.gz"
            ;;
        "Linux_aarch64")
            url="https://github.com/torognes/vsearch/releases/download/v2.13.3/vsearch-2.13.3-linux-aarch64.tar.gz"
            ;;
        "Linux_ppc64le")
            url="https://github.com/torognes/vsearch/releases/download/v2.13.3/vsearch-2.13.3-linux-ppc64le.tar.gz"
            ;;
        "Darwin_x86_64")
            url="https://github.com/torognes/vsearch/releases/download/v2.13.3/vsearch-2.13.3-macos-x86_64.tar.gz"
            ;;
        *)
            echo "Unsupported system architecture: ${os}_${arch}"
            return 1
            ;;
    esac

    temp_dir=$(mktemp -d)
    cd "$temp_dir" || return 1

    if ! wget "$url"; then
        echo "Failed to download vsearch binary. Please check your internet connection and try again."
        return 1
    fi

    tar_file=$(basename "$url")
    tar xzf "$tar_file"

    vsearch_binary=$(find . -name "vsearch" -type f -executable)

    if [ -z "$vsearch_binary" ]; then
        echo "Could not find the vsearch binary in the extracted files."
        return 1
    fi

    if sudo cp "$vsearch_binary" /usr/local/bin/vsearch; then
        echo "Vsearch has been successfully installed."
        vsearch --version
        return 0
    else
        echo "Failed to install vsearch. Please check the error messages and try again."
        return 1
    fi
}

install_catnip() {
    echo "Checking if catnit is installed."

    catnip_output=$(catnip --help | tail -n 1)

    if [ "$catnip_output" != 'catnip: command not found' ]; then
        echo "catnip is already installed."
    else
        echo "catnip is not installed. Installing catnip..."
        cd "$EXTERNAL_SCRIPTS_DIR" || { echo "Directory $EXTERNAL_SCRIPTS_DIR does not exist"; exit 1; }

        if [ $# -gt 0 ]; then
            TOKEN="$1"
            CATNIP_REPO_URL="https://${TOKEN}github.com/CIBIO-BU/catnip"
            echo "Using token-based repository UR for catnip."
            git clone "$CATNIP_REPO_URL" || { echo "Failed to clone catnip repository"; exit 1; }
        else
            echo "Using default SSH repository URL for catnip."
            git clone https://github.com/CIBIO-BU/catnip || { echo "Failed to clone catnip repository"; exit 1; }
        fi

        cd catmip || { echo "Directory catnip does not exist"; exit 1; }

        conda env create -f catnip-env.yml || { echo "Failed to create catnip Conda environment"; exit 1; }

        echo "Activating catnip Conda environment."
        conda activate catnip_env || { echo "Failed to activate catnip Conda environment"; exit 1; }

        echo "Installing catnip package."
        pip install -e . || { echo "Failed to install catnip package"; exit 1; }

        echo "catnip installation complete."
    fi
}

install_entry_points() {
    echo "Installing CLI commands."
    pip install -e .
}

main() {
    check_conda
    check_env
    activate_env
    clone_repo
    install_package
    install_catnip
    install_entry_points
    echo "mozaiko requires CRABS (v0.1.7), cutadapt and vsearch for downstream analysis."
    echo "Proceeding with installation..."
    install_crabs_release
    install_cutadapt_package
    echo "Instalation complete"
}

main