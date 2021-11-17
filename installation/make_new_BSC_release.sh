#!/bin/bash
set -e

# makes a new bsc release with the same base environment (only changing the scripts of perSVade)

# get the version
version=$1

# go to the releases folder
cd /home/mschikora/samba/scripts/perSVade/releases

# define dirs
code_gz="$version".tar.gz
code_dir=perSVade-"$version"
url="https://github.com/Gabaldonlab/perSVade/archive/refs/tags/$version.tar.gz"

# download the release
if [ ! -d $code_dir ]
then

	# download
	wget $url

	# extract
	tar -xvf $code_gz

	# remove
	rm $code_gz

fi

# get to the dir
cd $code_dir

# install extra dependencies
./installation/install_standalone_packages.sh

echo "perSVade version $version was created. You can now test it in the cluster with 'cd /gpfs/projects/bsc40/mschikora/scripts/perSVade/releases/perSVade-$version' and 'installation/test_installation/test_installation_modules.py &' "