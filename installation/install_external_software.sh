#!/bin/bash
set -e

# This script installs the external software

# define the path to the installation dir
installation_dir=$(readlink -f $0 | sed 's/\/install_external_software.sh//g')

# go to the installation dir
cd $installation_dir
if [ -d external_software ] 
then 
	rm -rf external_software
fi
mkdir external_software
cd external_software

# install the optimised lowess package
git clone https://github.com/histed/lowess-work-carljv
cd lowess-work-carljv
python setup.py install

# back to the external_software dir
cd $installation_dir/external_software

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
wget https://raw.githubusercontent.com/weiyuchung/CONY/master/CONY.R

# give execution permission to all
chmod 755 *;

echo 'SUCCESS: The non-conda external dependencies have been installed'