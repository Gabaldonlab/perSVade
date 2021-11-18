#!/bin/bash
set -e

# this script gets a $FROM and a $TO versions of perSVade, and creates a $TO version based on the $FROM version, changing the content of ./scripts and ./installation/test_installation

# get the version to build
from_image=mikischikora/persvade:$1
to_image=mikischikora/persvade:$2

# check that the from image exists
docker run $from_image ls > /dev/null 2>&1 || (echo "FROM image $from_image shoud exist" && exit 1)

# check whether the image already exists
docker run $to_image ls > /dev/null 2>&1 && echo "TO image $from_image already exists. Skipping it's generation"

# check that the to image does not exist
echo "Creating image $to_image if it does not exist already..."
docker run $to_image ls > /dev/null 2>&1  ||  docker build -t $to_image --build-arg from_tag=$1 -f ./docker/Dockerfile_onlyChangeScripts . 

# check that the image works
echo 'Testing the image'
mkdir $PWD/docker/perSVade_testing_outputs_$2 || echo 'perSVade_testing_outputs already exists'
docker run -v "$PWD/docker/perSVade_testing_outputs_$2":/perSVade/installation/test_installation/testing_outputs $to_image python -u /perSVade/installation/test_installation/test_installation_modules.py

echo "SUCCESS: You could create the docker image of perSVade. You can publish this image with docker push $to_image"

# print the steps to reproduce the singularity image
singularity_image="./mikischikora_persvade_$2.sif"
echo "You can then create a singularity image with singularity build --docker-login $singularity_image docker://$to_image"
echo "You can test that the singularity image works correctly with mkdir perSVade_testing_outputs; singularity exec -B ./perSVade_testing_outputs:/perSVade/installation/test_installation/testing_outputs -e $singularity_image bash -c 'source /opt/conda/etc/profile.d/conda.sh && conda activate perSVade_env && /perSVade/installation/test_installation/test_installation_modules.py' "
