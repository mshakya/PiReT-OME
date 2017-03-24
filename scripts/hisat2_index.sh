#!/bin/env bash


# if there is no eukarya_fasta variable assigned
if [ -z ${eukarya_fasta+x} ]
	then
	echo "creating hisat2 index for prokaryotic reference sequences"
	echo "hisat2-build --large-index -p $numCPU $prokaryote_fasta $ref_index"
	hisat2-build -p $numCPU --large-index $prokaryote_fasta  $ref_index
	touch $ref_index".done"
elif [ -z ${prokaryote_fasta+x} ] 
	then
	echo "creating hisat2 index for eukaryotic reference sequences" 
	hisat2-build -p $numCPU --large-index $eukarya_fasta $ref_index
	touch $ref_index".done"
else
	echo "creating hisat2 indices for both eukaryotic and prokaryotic references sequences" 
	echo "hisat2-build --large-index -p $numCPU $eukarya_fasta,$prokaryote_fasta  $ref_index" 
	hisat2-build --large-index -p $numCPU $eukarya_fasta,$prokaryote_fasta  $ref_index 
	touch $ref_index".done"
fi

