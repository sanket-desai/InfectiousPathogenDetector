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
  rm -rf tabix
  git clone --recursive https://github.com/samtools/tabix.git
  cd tabix && make && cd $home_folder
fi

if [[ -d $home_folder/external/htslib-1.10 ]];
then
	cd $home_folder/external/
  rm -rf htslib
  git clone --recursive https://github.com/samtools/htslib.git
  cd htslib && ./configure && make
	cd $home_folder
else
  cd $home_folder/external/
  rm -rf htslib
  git clone --recursive https://github.com/samtools/htslib.git
  cd htslib && ./configure && make
	cd $home_folder
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

#cd $home_folder/external/
#rm -rf bamtools
#git clone git://github.com/pezmaster31/bamtools.git
#cd bamtools
#mkdir build
#cd build
#cmake ../
#make
#cd $home_folder
#freebayes (freebayes requires g++, camke, the standard C and C++ development libraries, liblzma, pthread, and libbzip2.)

rm -rf $home_folder/external/freebayes
git clone --recursive https://github.com/ekg/freebayes.git
cd freebayes && make

if [[ ! -f "$home_folder/external/freebayes/bin/freebayes" ]];
then
	#cd $home_folder/external/
	#git clone --recursive https://github.com/vcflib/vcflib.git
	#cd $home_folder/external/vcflib
	#echo "$SUDO_PASSWWORD" | sudo -S make -j
  #make -j
	#cd $home_folder/external/
	#rm -rf $home_folder/external/freebayes/vcflib
	#mv $home_folder/external/vcflib $home_folder/external/freebayes
  #cd freebayes && make
  echo "Freebayes installation failed!!"
  abort
fi

# minimap2

if [[ ! -f "$home_folder/external/minimap2/minimap2" ]];
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
	conda create -n medaka -c conda-forge -c bioconda medaka
fi

#nanofile

pip_version=$(which pip3)
if [[ -z "$pip_version" ]]
then
	echo "pip3 not found in path"
	abort
else
	sudo pip3 install nanofilt
fi

#lofreq

if [[ ! -f "$home_folder/external/lofreq_star-2.1.2/bin/lofreq" ]]
then
	cd $home_folder/external/lofreq_star-2.1.2 && ./bootstrap && ./configure && make clean && make && make install && cd $home_folder
fi

shell_name=$(echo $0)

$shell_name databasesetup.sh
