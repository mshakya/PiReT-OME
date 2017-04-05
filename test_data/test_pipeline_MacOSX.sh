#!/usr/bin/env bash

find fastqs -iname "*test*R2.fastq" -exec echo `pwd`/{} \; > testR2.txt
find fastqs -iname "*test*R1.fastq" -exec echo `pwd`/{} \; | sed 's/R1\.fastq/R1\.fastq:/g' > testR1.txt
printf "samp1\nsamp2\nsamp3\nsamp4\nsamp5\nsamp6\n" > test_id.txt
printf "spleen\nspleen\nliver\nliver\nliver\nspleen\n" > test_gr.txt
paste -d "\0" testR1.txt testR2.txt > testR1R2.txt
paste test_id.txt testR1R2.txt test_gr.txt > test_ed.txt
rm test_id.txt testR1.txt testR2.txt testR1R2.txt test_gr.txt
echo -e "ID\tRawreads_files\tgroup" | cat - test_ed.txt > test_experimental_design.txt
rm test_ed.txt



printf "running pipeline now......\n"

printf "runPiReT-OME -test_kingdom eukarya \
-significant_pvalue 0.5 -exp test_experimental_design.txt \
-d pipeline_test_euk \
-eukarya_fasta data/eukarya_test.fa -index_ref_bt2 test_index \
-gff_eukarya data/eukarya_test.gff3"

runPiReT-OME -test_kingdom eukarya \
-significant_pvalue 0.5 -exp test_experimental_design.txt \
-d pipeline_test_euk \
-eukarya_fasta data/eukarya_test.fa -index_ref_bt2 test_index \
-gff_eukarya data/eukarya_test.gff3

printf "running pipeline finished"