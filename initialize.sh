#!/bin/bash

find_in_conda_env(){
    conda env list | grep "${@}" >/dev/null 2>/dev/null
}

SOURCEDIR=$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )

X=${SOURCEDIR}/chisel/src/chisel
conda env config vars set CHISEDIR=$X -n SEACON

# Create helper chisel environment
help_env_name=SEACON_chisel
if find_in_conda_env "$help_env_name"; then
    echo "$help_env_name already exists"
else
    echo "setting up helper chisel environment"
    conda create -n SEACON_chisel
    conda install -n SEACON_chisel chisel
fi

