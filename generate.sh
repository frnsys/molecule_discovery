#!/bin/bash

# Activate virtualenv
source ~/.pyenv/versions/data/bin/activate

while [ 1 ]; do
    # This one has to GPU unfortunately
    echo "Sampling..."
    CUDA_VISIBLE_DEVICES=0 python sample.py

    # Just use CPU here
    echo "Processing samples..."
    CUDA_VISIBLE_DEVICES='' python process_samples.py

    echo "Generating synthesis plans..."
    python plan.py

    echo "Generating synthesis plan images..."
    python make_plan_images.py

    echo "Updating site..."
    cd site
    python gendata.py
    bash update.sh
    cd ..

    echo "Done"
    sleep 600
done