#!/bin/bash -i

#################################################
# Date 			    : 29-09-2020
# Description   : set up the tools; create primary and secondary databases; index all the databases
# Author	      : Sanket Desai
# Updated on    : 29/09/2020
#################################################

echo >&2 '
**************************************
*** Download the databases 		******
**************************************
'


mkdir -p $home_folder/data/{annotation,cov2moduleref,primaryref/{hspathoref,pathoref},secondaryref}

cd $home_folder/src
if [[ -f $home_folder/data/annotation/patho.ids ]]
then
	python3 ipdcreatedb.py
else
	echo "patho.ids Files not found !!!!!!!"
        abort
fi

cd $home_folder

LATEST=$(curl -s4 ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/ | grep -Eo 'release_[0-9]+' | sed 's/release_//g' | sort -nr | head -n1)

wget "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$LATEST/GRCh38.primary_assembly.genome.fa.gz" -P $home_folder/data/primaryref/hspathoref
wget "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_$LATEST/gencode.v$LATEST.annotation.gff3.gz" -P $home_folder/data/primaryref/hspathoref

gencode_fasta_file=$home_folder/data/primaryref/hspathoref/GRCh38.primary_assembly.genome.fa.gz
pathogen_fasta=$home_folder/data/primaryref/pathoref/patho.fa
if [ -f "$gencode_fasta_file" ] && [ -f "$pathogen_fasta" ];
then
	zcat $gencode_fasta_file > $home_folder/data/primaryref/hspathoref/hspatho.fa
	cat $pathogen_fasta >> $home_folder/data/primaryref/hspathoref/hspatho.fa
else
	echo "Fasta Files not found !!!!!!!"
	abort
fi

gencode_gff_file=$home_folder/data/primaryref/hspathoref/gencode.v$LATEST.annotation.gff3.gz
pathogen_gff=$home_folder/data/annotation/patho.gff

if [ -f "$gencode_gff_file" ] && [ -f "$pathogen_gff" ];
then
	zcat $gencode_gff_file > data/annotation/hspatho.gff
	cat $pathogen_gff >> data/annotation/hspatho.gff
else
	echo "Gff Files not found !!!!!!!"
	abort
fi

echo >&2 '
**************************************
*** Create Version document 	******
**************************************
'


date > $home_folder/version.info
echo -e "Hg38 Gencode version""\t"$LATEST >> $home_folder/version.info
echo -e "Number of Pathogens Included in Primary Database""\t"$(wc -l data/annotation/pathoannotation.tsv | awk '{print $1}') >> $home_folder/version.info
echo -e "Number of Pathogens Included in Secondary Database""\t"$(wc -l data/annotation/secondaryannotation.tsv | awk '{print $1}') >> $home_folder/version.info

rm $gencode_fasta_file $gencode_gff_file
#if [[ -f  $home_folder/data/annotation/patho.ids.gbk ]];
#then #
#	cat $home_folder/data/annotation/patho.ids.gbk >> $home_folder/external/snpEff/data/ipd1060/genes.gbk
#else
#	echo "patho.ids.gbk file not found."
#	abort
#fi
cd $home_folder/src
if [[ -f $home_folder/data/annotation/patho.ids.gbk ]]
then
	python3 snpeffinputgbkcreator.py
else
	echo "patho.ids.gbk Files not found !!!!!!!"
        abort
fi

cd $home_folder


echo >&2 '
**************************************
*** Indexing the ipd database 	******
**************************************
'

if [[ -f $home_folder/data/primaryref/pathoref/patho.fa ]];
then

	#hisat index
	$home_folder/external/hisat2-2.1.0/hisat2-build -p 4 $home_folder/data/primaryref/pathoref/patho.fa $home_folder/data/primaryref/pathoref/patho
	#samtools faidx
	$home_folder/external/samtools-1.10/samtools faidx $home_folder/data/primaryref/pathoref/patho.fa
	#minimap2 index
	$home_folder/external/minimap2/minimap2 -d $home_folder/data/primaryref/pathoref/patho.minimap2.ref $home_folder/data/primaryref/pathoref/patho.fa
	cnt=$(ls $home_folder/data/primaryref/pathoref/patho.* | wc -l)
	if [[ $cnt == 0 ]]
	then
		echo "Patho Hisat Index not generated"
		abort
	fi
else
	echo "Patho database not created"
	abort
fi

#indexing hspatho.fa (primary db + human)
if [[ -f $home_folder/data/primaryref/hspathoref/hspatho.fa ]];
then

	#hisat index
	$home_folder/external/hisat2-2.1.0/hisat2-build -p 4 $home_folder/data/primaryref/hspathoref/hspatho.fa $home_folder/data/primaryref/hspathoref/hspatho
	cnt=$(ls $home_folder/data/primaryref/hspathoref/hspatho.* | wc -l)
	if [[ $cnt == 0 ]]
	then
		echo "Hspatho Hisat Index not generated"
		abort
	fi
else
	echo "HsPatho database not created"
	abort
fi

#indexing secondary.fa (secondary db)
if [[ -f $home_folder/data/secondaryref/secondary.fa ]];
then
	#blast index
	$home_folder/external/ncbi-blast-2.10.0+/bin/makeblastdb -in $home_folder/data/secondaryref/secondary.fa -input_type fasta -dbtype nucl
	$home_folder/external/ncbi-blast-2.10.0+/bin/makembindex -input $home_folder/data/secondaryref/secondary.fa -output $home_folder/data/secondaryref/secondary.fa
	cnt=$(ls $home_folder/data/secondaryref/secondary.fa.* | wc -l)
	if [[ $cnt == 0 ]]
	then
		echo "Secondary Index not generated"
		abort
	fi
else
	echo "secondary database not created"
	abort
fi

#building snpEff file
if [[ -f $home_folder/external/snpEff/data/ipd1060/genes.gbk ]];
then
	cd $home_folder/external/snpEff
	java -Xmx40g -jar snpEff.jar build -genbank -v ipd1060
else
	echo "$home_folder/external/snpEff/data/ipd1060/genes.gbk doesn't exist"
	abort
fi

cd $home_folder

trap : 0

echo >&2 '
**************************************
*** Installation Successfully DONE ***
**************************************
'
