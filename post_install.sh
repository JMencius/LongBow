#!/bin/bash

echo "Setting up post installtion"
BASE=$(conda info --base)
mkdir -p "$BASE"/envs/ont-longbow/etc/conda/activate.d;
echo "alias longbow='python $PWD/longbow.py'" > "$BASE"/envs/ont-longbow/etc/conda/activate.d/setup_alias.sh;
echo "Post installation finished";
echo "To test LongBow installation, use";
echo " "
echo "$ conda activate ont-longbow;longbow --version";
