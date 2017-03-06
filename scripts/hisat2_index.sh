#!/bin/sh

if [ -z ${eukarya_fasta+x} ] && [ -z ${prokaryote_fasta+x} ];
	then
	hisat2-build --large-index -p $numCPU $eukarya_fasta,$prokaryote_fasta  $ref_index;
	touch "$ref_index'.done'";
elif [ -z ${eukarya_fasta+x} ];
	then	
	hisat2-build -p $numCPU --large-index $prokaryote_fasta  $ref_index;
	touch "$ref_index'.done'";
elif [ -z ${prokaryote_fasta+x} ];
	then
	hisat2-build -p $numCPU --large-index $eukarya_fasta  $ref_index;
	touch "$ref_index'.done'";
fi

