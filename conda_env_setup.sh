#!/usr/bin/env bash
# -*- coding: utf-8 -*-

#
# create a conda environment
ENV="$1"

if [ "$ENV-" == "-" ]; then
    CDIR=$(basename "$PWD")
    ENV="${CDIR}_env"
fi
set  -o pipefail
if conda env list | grep ".*${ENV}.*" >/dev/null 2>&1; then
    echo "conda $ENV already exists, skipping creation..."
else
    echo "env $ENV does not exist, creating it..."
    conda create -n "$ENV" -y
fi

# Next line is needed in otder to be able to activate the environment
CUR_SHELL=shell.$(basename -- "${SHELL}")
eval "$(conda "$CUR_SHELL" hook)"

#
set -e
conda activate "$ENV"
echo "INFO: conda environment $ENV created and activated"

REPOS=(-c r -c conda-forge)

# Packages to install
conda install -n "$ENV" -y "${REPOS[@]}" shellcheck
# complete here


echo "All done."
exit 0
