#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ENV_NAME="mozaiko"
REPO_URL="git@github.com:CIBIO-BU/mozaiko.git"
CATNIP_REPO_URL="https://github.com/CIBIO-BU/catnip.git"
PACKAGE_DIR="mozaiko"
CRABS_RELEASE="https://github.com/gjeunen/reference_database_creator/archive/refs/tags/v0.1.7.tar.gz"
EXTERNAL_SCRIPTS_DIR="${SCRIPT_DIR}/${PACKAGE_DIR}/external_scripts"
CRABS_ARCHIVE="crabs.tar.gz"
CRABS_DIR="reference_database_creator-0.1.7"

# Color codes for better output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Parse command line arguments
parse_args() {
    if [ $# -gt 0 ]; then
        TOKEN="$1"
        REPO_URL="https://${TOKEN}@github.com/CIBIO-BU/mozaiko.git"
        CATNIP_REPO_URL="https://${TOKEN}@github.com/CIBIO-BU/catnip.git"
        log_info "Using token-based repository URL"
    else
        log_info "Using default SSH repository URL"
    fi
    log_info "Repository URL: $REPO_URL"
}

# Check if Conda is installed
check_conda() {
    if ! command -v conda &> /dev/null; then
        log_error "Conda could not be found. Please install Conda first."
        echo "Visit: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
    log_info "Conda found: $(conda --version)"
}

# Initialize conda for bash
init_conda() {
    log_info "Initializing Conda for current shell"

    # Get conda base directory
    CONDA_BASE=$(conda info --base)

    # Source conda.sh to enable conda activate
    if [ -f "${CONDA_BASE}/etc/profile.d/conda.sh" ]; then
        source "${CONDA_BASE}/etc/profile.d/conda.sh"
    else
        log_error "Could not find conda.sh"
        exit 1
    fi
}

# Check if environment already exists
check_env() {
    if conda env list | grep -q "^${ENV_NAME} "; then
        log_warn "Conda environment '${ENV_NAME}' already exists"
        read -p "Do you want to remove and recreate it? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            log_info "Removing existing environment"
            conda env remove -n "${ENV_NAME}" -y
            return 1
        else
            log_info "Using existing environment"
            return 0
        fi
    fi
    return 1
}

# Create Conda environment with all dependencies
create_env() {
    log_info "Creating Conda environment: ${ENV_NAME}"

    conda create -y -n "${ENV_NAME}" \
        -c bioconda -c conda-forge \
        python=3.12 \
        pip \
        vsearch=2.13.3 \
        cutadapt=4.9 \
        biopython=1.84 \
        numpy=2.0.0 \
        pandas=2.2.2 \
        matplotlib=3.9.1 \
        tqdm \
        pyyaml \
        requests \
        wget \
        git || {
            log_error "Failed to create Conda environment"
            exit 1
        }

    log_info "Conda environment created successfully"
}

# Activate Environment
activate_env() {
    log_info "Activating Conda environment: ${ENV_NAME}"
    conda activate "${ENV_NAME}" || {
        log_error "Failed to activate environment"
        exit 1
    }
    log_info "Environment activated successfully"
}

# Clone repository
clone_repo() {
    if [ -d "$PACKAGE_DIR" ]; then
        log_warn "Directory $PACKAGE_DIR already exists"
        read -p "Do you want to update it? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            log_info "Updating repository"
            cd "$PACKAGE_DIR" || exit 1
            git pull || {
                log_error "Failed to update repository"
                exit 1
            }
            cd "$SCRIPT_DIR"
        else
            log_info "Skipping repository update"
        fi
    else
        log_info "Cloning repository from $REPO_URL"
        git clone "$REPO_URL" || {
            log_error "Failed to clone repository"
            exit 1
        }
    fi
}

# Install Package
install_package() {
    log_info "Installing mozaiko package"
    cd "$PACKAGE_DIR" || {
        log_error "Directory $PACKAGE_DIR does not exist"
        exit 1
    }

    pip install -e . || {
        log_error "Failed to install package"
        exit 1
    }

    log_info "Mozaiko package installation complete"
}

# Install CRABS v0.1.7
install_crabs_release() {
    log_info "Checking CRABS installation"

    if command -v crabs &> /dev/null; then
        crabs_output=$(crabs --version 2>&1 | tail -n 1)
        crabs_version=${crabs_output##* }

        if [ "$crabs_version" = "0.1.7" ]; then
            log_info "CRABS v0.1.7 is already installed"
            return 0
        else
            log_warn "CRABS version $crabs_version is installed, but v0.1.7 is required"
            log_warn "Please manually remove CRABS and rerun this script"
            return 1
        fi
    fi

    log_info "Installing CRABS v0.1.7"

    mkdir -p "$EXTERNAL_SCRIPTS_DIR" || {
        log_error "Failed to create directory $EXTERNAL_SCRIPTS_DIR"
        exit 1
    }

    cd "$EXTERNAL_SCRIPTS_DIR" || {
        log_error "Directory $EXTERNAL_SCRIPTS_DIR does not exist"
        exit 1
    }

    log_info "Downloading CRABS v0.1.7"
    wget "$CRABS_RELEASE" -O "$CRABS_ARCHIVE" || {
        log_error "Failed to download CRABS"
        exit 1
    }

    log_info "Extracting CRABS archive"
    tar -xzf "$CRABS_ARCHIVE" || {
        log_error "Failed to extract CRABS"
        exit 1
    }

    cd "$CRABS_DIR" || {
        log_error "Directory $CRABS_DIR does not exist"
        exit 1
    }

    log_info "Installing CRABS package"
    pip install . || {
        log_error "Failed to install CRABS"
        exit 1
    }

    cd "$EXTERNAL_SCRIPTS_DIR"
    rm -f "$CRABS_ARCHIVE"

    log_info "CRABS installation complete"
}

# Verify tools installation
verify_tools() {
    log_info "Verifying tool installations"

    local all_good=true

    if command -v cutadapt &> /dev/null; then
        log_info "cutadapt: $(cutadapt --version)"
    else
        log_error "cutadapt not found"
        all_good=false
    fi

    if command -v vsearch &> /dev/null; then
        log_info "vsearch: $(vsearch --version 2>&1 | head -n1)"
    else
        log_error "vsearch not found"
        all_good=false
    fi

    if command -v crabs &> /dev/null; then
        log_info "crabs: $(crabs --version 2>&1 | tail -n1)"
    else
        log_warn "crabs not found (optional)"
    fi

    if [ "$all_good" = false ]; then
        log_error "Some required tools are missing"
        return 1
    fi

    log_info "All required tools verified successfully"
    return 0
}

# Install catnip
install_catnip() {
    log_info "Checking catnip installation"

    if conda env list | grep -q "^catnip "; then
        log_info "catnip environment already exists"
        return 0
    fi

    log_info "Creating catnip environment"

    mkdir -p "$EXTERNAL_SCRIPTS_DIR"
    cd "$EXTERNAL_SCRIPTS_DIR" || {
        log_error "Directory $EXTERNAL_SCRIPTS_DIR does not exist"
        exit 1
    }

    if [ ! -d "catnip" ]; then
        log_info "Cloning catnip repository"
        git clone "$CATNIP_REPO_URL" || {
            log_error "Failed to clone catnip repository"
            exit 1
        }
    fi

    cd catnip || {
        log_error "Directory catnip does not exist"
        exit 1
    }

    conda env create -f catnip-env.yml || {
        log_error "Failed to create catnip Conda environment"
        exit 1
    }

    conda run -n catnip pip install -e . || {
        log_error "Failed to install catnip package"
        exit 1
    }

    log_info "catnip installation complete"
    log_info "To use catnip, run: conda activate catnip"
}

# Main installation function
main() {
    log_info "Starting mozaiko installation"

    parse_args "$@"
    check_conda
    init_conda

    local env_exists=false
    if check_env; then
        env_exists=true
    fi

    if [ "$env_exists" = false ]; then
        create_env
    fi

    activate_env
    clone_repo
    install_package

    log_info "Installing additional dependencies"
    install_crabs_release
    install_catnip

    verify_tools

    log_info "Installation complete!"
    log_info "To use mozaiko, run: conda activate ${ENV_NAME}"
}

main "$@"