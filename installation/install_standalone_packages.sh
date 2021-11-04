#!/bin/bash
set -e

# This script adds external software so that perSVade.py can run. It does not create any environmental changes

# define installation dir
installation_dir=$(readlink -f $0 | sed 's/\/install_standalone_packages.sh//g')

# go to the installation dir and make the external software one
cd $installation_dir
if [ -d external_software ] 
then 
	rm -rf external_software
fi
mkdir external_software
cd external_software

# download the gztool depending on the architecture
architechture=$(uname -i)

if [ $architechture = "x86_64" ]
then
	echo "installing gztool for architechture $architechture ..."
	wget https://github.com/circulosmeos/gztool/releases/download/v0.11.5/gztool-linux.x86_64
	mv gztool-linux.x86_64 gztool
	
elif [ $architechture = "aarch64" ]
then
	echo "installing gztool for architechture $architechture"
	wget https://github.com/circulosmeos/gztool/releases/download/v0.11.5/gztool-linux.aarch64
	mv gztool-linux.aarch64 gztool

else
	echo "WARNING: gztool can't be installed in architechture $architechture. perSVade will use more standard and slower tools instead. "
fi

# download extra dependencies (platform-independent)
wget https://github.com/PapenfussLab/gridss/releases/download/v2.9.2/gridss-2.9.2-gridss-jar-with-dependencies.jar
wget https://github.com/PapenfussLab/gridss/releases/download/v2.9.2/gridss.sh
wget https://github.com/PapenfussLab/clove/releases/download/v0.17/clove-0.17-jar-with-dependencies.jar

# give execution permission to all
chmod -R 755 $installation_dir/../*

# print success file
echo 'SUCCESS: You already installed all the dependencies that should be added to this github repo files'