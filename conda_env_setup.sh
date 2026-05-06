#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

ENV_NAME="mozaiko"
REPO_URL="git@github.com:CIBIO-BU/mozaiko.git"
PACKAGE_DIR="mozaiko"
CRABS_RELEASE="https://github.com/gjeunen/reference_database_creator/archive/refs/tags/v0.1.7.tar.gz"
EXTERNAL_SCRIPTS_DIR="${SCRIPT_DIR}/external_scripts"
CRABS_ARCHIVE="crabs.tar.gz"
CRABS_DIR="reference_database_creator-0.1.7"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# # Parse command line arguments
# parse_args() {
#     if [ $# -gt 0 ]; then
#         TOKEN="$1"
#         REPO_URL="https://${TOKEN}@github.com/CIBIO-BU/mozaiko.git"
#         log_info "Using token-based repository URL"
#     else
#         log_info "Using default SSH repository URL"
#     fi
#     log_info "Repository URL: $REPO_URL"
# }

check_conda() {
    command -v conda >/dev/null 2>&1 || {
        log_error "Conda not found. Install Miniconda first."
        exit 1
    }
}

# Check if env already exists
check_env() {
    if conda env list | grep -q "${ENV_NAME}"; then
        echo "Conda environment '$ENV_NAME' already exists. Skipping environment creation."
        return 0
    else
        return 1
    fi
}

init_conda() {
    source "$(conda info --base)/etc/profile.d/conda.sh"
}

# Create Conda environment with all dependencies
create_env() {
    log_info "Creating Conda environment: ${ENV_NAME}"

    conda env create -f environment.yml || {
        log_error "Failed to create Conda environment"
        exit 1
    }

    log_info "Conda environment created successfully"
}

activate_env() {
    log_info "Activating environment: $ENV_NAME"
    conda activate "$ENV_NAME" || true

    # Verify activation actually worked
    if [ "${CONDA_DEFAULT_ENV:-}" != "$ENV_NAME" ]; then
        log_error "Failed to activate conda environment '$ENV_NAME'"
        exit 1
    fi

    log_info "Environment '$ENV_NAME' activated successfully"
}

# # Clone repository
# clone_repo() {
#     if [ -d "$PACKAGE_DIR" ]; then
#         log_warn "Directory $PACKAGE_DIR already exists"

#         # Check if running interactively
#         if [ -t 0 ]; then
#             read -p "Do you want to update it? (y/N): " -n 1 -r
#             echo
#         else
#             log_warn "Non-interactive shell detected, skipping update"
#             REPLY="n"
#         fi

#         if [[ $REPLY =~ ^[Yy]$ ]]; then
#             log_info "Updating repository"
#             cd "$PACKAGE_DIR" || exit 1
#             git pull || {
#                 log_error "Failed to update repository"
#                 exit 1
#             }
#             cd "$SCRIPT_DIR"
#         else
#             log_info "Skipping repository update"
#         fi
#     else
#         log_info "Cloning repository from $REPO_URL"
#         git clone "$REPO_URL" || {
#             log_error "Failed to clone repository"
#             exit 1
#         }
#     fi
# }

# Install Package
install_package() {
    log_info "Installing mozaiko package"
    # cd "$PACKAGE_DIR" || {
    #     log_error "Directory $PACKAGE_DIR does not exist"
    #     exit 1
    # }

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

    log_info "Installing CRABS v0.1.7 from local archive"

    cd "$EXTERNAL_SCRIPTS_DIR" || {
        log_error "Directory $EXTERNAL_SCRIPTS_DIR does not exist"
        exit 1
    }

    CRABS_ARCHIVE="$EXTERNAL_SCRIPTS_DIR/crabs-0.1.7.zip"

    if [ ! -f "$CRABS_ARCHIVE" ]; then
        log_error "CRABS archive not found at $CRABS_ARCHIVE"
        exit 1
    fi

    log_info "Unzipping CRABS archive"
    unzip -o "$CRABS_ARCHIVE" || {
        log_error "Failed to unzip CRABS archive"
        exit 1
    }

    CRABS_DIR="$EXTERNAL_SCRIPTS_DIR/crabs-0.1.7"

    cd "$CRABS_DIR" || {
        log_error "Directory $CRABS_DIR does not exist after unzipping."
        exit 1
    }

    log_info "Installing CRABS package"
    pip install . || {
        log_error "Failed to install CRABS"
        exit 1
    }

    log_info "CRABS installation complete."
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

    if command -v catnip &> /dev/null; then
        log_info "catnip: $(catnip --version 2>&1 | head -n1)"
    else
        log_warn "catnip not found"
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

# Main installation function
main() {
    log_info "Starting mozaiko installation"

    # parse_args "$@"
    init_conda

    local env_exists=false
    if check_env; then
        env_exists=true
    fi

    if [ "$env_exists" = false ]; then
        create_env
    fi

    activate_env
    # clone_repo
    install_package

    log_info "Installing additional dependencies"
    install_crabs_release

    verify_tools

    log_info "Installation complete!"
    log_info "To use mozaiko, run: conda activate ${ENV_NAME}"
}

main