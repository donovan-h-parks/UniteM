#!/bin/bash

# install python libraries
python -m pip install . --no-deps --ignore-installed --no-cache-dir -vvv

# copy main UniteM python script
mkdir -p ${PREFIX}/bin/

chmod +x bin/unitem
cp bin/unitem ${PREFIX}/bin/