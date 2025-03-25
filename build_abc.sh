#!/bin/bash

# Check if abc directory exists
if [ ! -d "abc" ]; then
    # Initialize and update abc submodule
    git submodule init
    git submodule update
fi

# Enter abc directory and build
cd abc
make -j$(nproc) 