#!/bin/bash -i

#################################################
# Date 			    : 01-09-2020 (v0_02)
# Description   : set up the tools; create primary and secondary databases; index all the databases
# Author	      : Sonal Rashmi
# Updated on    : 28/09/2020 Sanket
#################################################

abort()
{
    echo >&2 '
*****************************
*** Error in Installation ***
*****************************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

trap 'abort' 0

set -e

home_folder=$(pwd)

echo >&2 '
**************************************
*** Tools Installation 			******
**************************************
'
#fastp
if [[ ! -f $home_folder/external/fastp-0.20.1/fastp ]];
then
	cd $home_folder/external/fastp-0.20.1 && make clean && make && cd $home_folder
fi
#samtools

if [[ ! -f $home_folder/external/tabix/tabix  ]];
then
  cd $home_folder/external
  git clone --recursive https://github.com/samtools/tabix.git
  cd tabix && make && cd $home_folder
fi

if [[ ! -d $home_folder/external/htslib-1.10 ]];
then
	cd $home_folder/external/
	wget https://github.com/samtools/htslib/releases/download/1.10/htslib-1.10.tar.bz2
	tar -xf htslib-1.10.tar.bz2
	cd htslib-1.10
	make
	rm $home_folder/external/htslib-1.10.tar.bz2
	mv $home_folder/external/htslib-1.10 $home_folder/external/htslib
	cd $home_folder
else
	cd $home_folder/external/htslib && make && cd $home_folder
fi

if [[ ! -d $home_folder/external/samtools-1.10 ]];
then
	cd $home_folder/external/
	wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
	tar -xf samtools-1.10.tar.bz2
	cd samtools-1.10
	make
	rm $home_folder/external/samtools-1.10.tar.bz2
	cd $home_folder
else
	cd $home_folder/external/samtools-1.10/htslib-1.10 && make clean
	cd $home_folder/external/samtools-1.10 && make clean && ./configure && make && cd $home_folder
fi

if [[ ! -d "$home_folder/external/bcftools" ]];
then
	cd $home_folder/external/
	wget https://github.com/samtools/bcftools/releases/download/1.10/bcftools-1.10.tar.bz2
	tar -xf bcftools-1.10.tar.bz2
	cd bcftools-1.10
	make
	mv $home_folder/external/bcftools-1.10 $home_folder/external/bcftools
	rm $home_folder/external/bcftools-1.10.tar.bz2
	cd $home_folder
else
	cd $home_folder/external/bcftools/htslib-1.6 && make clean
	cd $home_folder/external/bcftools && make clean && ./configure && make && cd $home_folder

fi

# bamtools for freebayes

cd $home_folder/external/
git clone git://github.com/pezmaster31/bamtools.git
cd bamtools
mkdir build
cd build
cmake ../
make
cd $home_folder
#freebayes (freebayes requires g++, camke, the standard C and C++ development libraries, liblzma, pthread, and libbzip2.)

if [[ ! -d "$home_folder/external/freebayes" ]];
then
	cd $home_folder/external/
	wget http://clavius.bc.edu/~erik/freebayes/freebayes-5d5b8ac0.tar.gz
	tar xf freebayes-5d5b8ac0.tar.gz
  mv bamtools ./freebayes
	cd freebayes
	if [ ! -f $home_folder/external/freebayes/bin/freebayes ];
	then
    make -j4
		#echo "$SUDO_PASSWWORD" | sudo -S make && echo "$SUDO_PASSWWORD" | sudo -S make;
	else
		echo "Freebayes Installation Error !!!! Check the libraries."
		abort
	fi
	cd $home_folder/external/
	git clone --recursive https://github.com/vcflib/vcflib.git
	cd $home_folder/external/vcflib
	#echo "$SUDO_PASSWWORD" | sudo -S make -j
  make -j
	cd $home_folder/external/
	rm -rf $home_folder/external/freebayes/vcflib
	mv $home_folder/external/vcflib $home_folder/external/freebayes
else
	cd $home_folder/external/vcflib && make clean && make -j && cd $home_folder
	cd $home_folder/external/freebayes && make clean && make -j4 && cd $home_folder
fi

# minimap2

if [[ ! -f "$home_folder/external/minimap2" ]];
then
	cd $home_folder/external/minimap2 && make clean && make && cd $home_folder
fi

#medaka (conda)

conda_env=$(conda info | grep -i 'base environment')
medaka_env=$(conda info --envs | grep "medaka")
if [[ -z $conda_env ]];
then
	echo "Install Conda!!!!!!!!!!!!!"
	abort
elif [[ -z $medaka_env ]];
then
   source activate medaka
else
	conda create -n medaka -c conda-forge -c bioconda medaka
fi

#nanofile

pip_version=$(which pip3)
if [[ -z "$pip_version" ]]
then
	echo "pip3 not found in path"
	abort
else
	pip3 install nanofilt
fi

#lofreq

if [[ ! -f "$home_folder/external/lofreq_star-2.1.2/bin/lofreq" ]]
then
	cd $home_folder/external/lofreq_star-2.1.2 && ./bootstrap && ./configure && make clean && make && make install && cd $home_folder
fi

#bash databasesetup.sh
