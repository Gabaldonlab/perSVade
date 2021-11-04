#!/bin/bash
set -e

# This script adds external software so that perSVade.py can run

# define installation dir
installation_dir=$(readlink -f $0 | sed 's/\/setup_environment.sh//g')

