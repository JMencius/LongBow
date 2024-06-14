#!/bin/bash

source activate ont-longbow;
mkdir -p "$CONDA_PREFIX"/etc/conda/activate.d;
echo "alias longbow='python $PWD/longbow.py'" > "$CONDA_PREFIX"/etc/conda/activate.d/setup_alias.sh;
