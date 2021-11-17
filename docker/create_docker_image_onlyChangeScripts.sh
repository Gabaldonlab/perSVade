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

dajhgdad

# print the steps to reproduce the singularity image
singularity_image="./mikischikora_persvade_$2.sif"
echo "You can then create a singularity image with singularity build --docker-login $singularity_image docker://$to_image"
echo "You can test that the singularity image works correctly with mkdir perSVade_testing_outputs; singularity exec -B ./perSVade_testing_outputs:/perSVade/installation/test_installation/testing_outputs -e $singularity_image bash -c 'source /opt/conda/etc/profile.d/conda.sh && conda activate perSVade_env && /perSVade/installation/test_installation/test_installation_modules.py' "


jgadjhgahgd 










# create a function that returns an error if the image does not exist
function check_image_exists {

   docker run $1 ls > /dev/null 2>&1 && echo "image $1 already exists"

}

# create base image with the environment
base_image=persvade_base:$tag
check_image_exists $base_image || {

    echo "$base_image does not exist. Generating"
    docker build -t $base_image -f ./docker/Dockerfile_base .

}

# add external_software dependencies
external_software_image=persvade_with_external_software:$tag
check_image_exists $external_software_image || {

    echo "$external_software_image does not exist. Generating"
    docker build -t $external_software_image --build-arg tag=$tag -f ./docker/Dockerfile_external_software . 

}

# add repeats dependencies
repeats_image=persvade_with_repeats:$tag
check_image_exists $repeats_image || {

    echo "$repeats_image does not exist. Generating"
    docker build -t $repeats_image --build-arg tag=$tag -f ./docker/Dockerfile_repeats . 

}

# add gridss_env
gridss_env_image=persvade_with_gridss_env:$tag
check_image_exists $gridss_env_image || {

    echo "$gridss_env_image does not exist. Generating"
    docker build -t $gridss_env_image --build-arg tag=$tag -f ./docker/Dockerfile_gridss_env . 

}

# add picard_env
picard_env_image=persvade_with_picard_env:$tag
check_image_exists $picard_env_image || {

    echo "$picard_env_image does not exist. Generating"
    docker build -t $picard_env_image --build-arg tag=$tag -f ./docker/Dockerfile_picard_env . 

}

# add all extra dependencies
final_image=mikischikora/persvade:$tag
check_image_exists $final_image || {

    echo "$final_image does not exist. Generating"
    docker build -t $final_image --build-arg tag=$tag -f ./docker/Dockerfile_extra_envs . 

}

# check that the image works
echo 'Testing the image'
mkdir $PWD/docker/perSVade_testing_outputs_$tag || echo 'perSVade_testing_outputs already exists'
docker run -v $PWD/docker/perSVade_testing_outputs_$tag:/perSVade/installation/test_installation/testing_outputs mikischikora/persvade:$tag python -u /perSVade/installation/test_installation/test_installation_modules.py

# success message
echo "SUCCESS: You could create the docker image of perSVade. You can publish this image with docker push mikischikora/persvade:$tag. You can then create a singularity image with singularity build --docker-login ./mikischikora_persvade_$tag.sif docker://mikischikora/persvade:$tag"