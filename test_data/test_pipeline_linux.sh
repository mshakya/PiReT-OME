#!/usr/bin/env bash

# create a test experimental design file

find fastqs -iname "*test*R2.fastq" -exec readlink -f {} \; | sort > testR2.txt
find fastqs -iname "*test*R1.fastq" -exec readlink -f {} \; | sort | sed 's/R1\.fastq/R1\.fastq:/g' > testR1.txt
printf "samp1\nsamp2\nsamp3\nsamp4\nsamp5\nsamp6\n" > test_id.txt
printf "liver\nspleen\nspleen\nliver\nliver\nspleen\n" > test_gr.txt
paste test_id.txt testR1.txt testR2.txt test_gr.txt > test_ed.txt
rm test_id.txt testR1.txt testR2.txt test_gr.txt
sed -i 's/:\t/:/g' test_ed.txt
echo -e "ID\tRawreads_files\tgroup" | cat - test_ed.txt > test_experimental_design.txt
rm test_ed.txt


#
printf "running pipeline now"

printf "runPiReT-OME \
-exp test_experimental_design.txt \
-d pipeline_test_euk \
-eukarya_fasta data/eukarya_test.fa -index_ref_bt2 test_index \
-gff_eukarya data/eukarya_test.gff3"

runPiReT-OME \
-exp test_experimental_design.txt \
-d pipeline_test_euk \
-eukarya_fasta data/eukarya_test.fa -index_ref_ht2 test_index \
-gff_eukarya data/eukarya_test.gff3

printf "running pipeline finished"