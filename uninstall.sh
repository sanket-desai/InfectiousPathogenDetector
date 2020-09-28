#!/bin/bash -i

#################################################
# Date 			    : 01-09-2020 (v0_02)
# Description   : Uninstall the tools and
# Author	      : Sonal Rashmi
# Updated on    : 28/09/2020 - Sanket 
#################################################

abort()
{
    echo >&2 '
*******************************
*** Error in Uninstallation ***
*******************************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

trap 'abort' 0

set -e

home_folder=$(pwd)

echo >&2 '
****************************************
*** Deleting External Tools 		****
****************************************
'


#samtools
if [[ -d $home_folder/external/samtools-1.10 ]];
then
	rm -rf $home_folder/external/samtools-1.10
fi

#bcftools
if [[ -d $home_folder/external/bcftools ]];
then
	rm -rf $home_folder/external/bcftools
fi

#htslib
if [[ -d $home_folder/external/htslib ]];
then
	rm -rf $home_folder/external/htslib
fi

#freebayes

if [[ -d $home_folder/external/freebayes ]];
then
	rm -rf $home_folder/external/freebayes*
fi

#vcflib

if [[ -d $home_folder/external/vcflib ]];
then
        rm -rf $home_folder/external/vcflib
fi


# minimap2

if [[ -d $home_folder/external/minimap2 ]];
then
	cd $home_folder/external/minimap2 && make clean && cd $home_folder
fi

# tabix

if [[ -d $home_folder/external/tabix ]];
then
  cd $home_folder/external/tabix && make clean && cd $home_folder
fi

#lofreq
if [[ -z $home_folder/external/lofreq_star-2.1.2/bin/lofreq ]]
then
	cd $home_folder/external/lofreq_star-2.1.2 && make clean cd $home_folder
fi


echo >&2 '
****************************************
*** Deleting the ipd database 		****
****************************************
'

if [[ -f $home_folder/data/primaryref/hspathoref/hspatho.fa ]]
then
	rm $home_folder/data/primaryref/hspathoref/*
fi

if [[ -f $home_folder/data/primaryref/pathoref/patho.fa ]]
then
	rm $home_folder/data/primaryref/pathoref/*
fi

if [[ -f $home_folder/data/secondaryref/secondary.fa ]]
then
	rm $home_folder/data/secondaryref/*
fi

if [[ -f $home_folder/data/annotation/hspatho.gff ]]
then
	rm $home_folder/data/annotation/*.gff
	rm $home_folder/data/annotation/*.gbk
	rm $home_folder/data/annotation/*.tsv

fi

if [[ -d $home_folder/data/annotation/gbktemp ]]
then
	rm -rf $home_folder/data/annotation/gbktemp/*
fi


if [[ -f $home_folder/version.info ]]
then
	rm $home_folder/version.info
fi

if [[ -f $home_folder/external/snpEff/data/ipd1060/genes.gbk ]]
then
	rm $home_folder/external/snpEff/data/ipd1060/genes.gbk
fi

cd $home_folder

trap : 0

echo >&2 '
****************************************
*** Uninstallation Successfully DONE ***
****************************************
'
