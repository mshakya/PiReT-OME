# create a test experimental design file

find fastqs -iname "*test*R2.fastq" -exec readlink -f {} \; > testR2.txt
find fastqs -iname "*test*R1.fastq" -exec readlink -f {} \; | sed 's/R1\.fastq/R1\.fastq:/g' > testR1.txt
printf "samp1\nsamp2\nsamp3\n" > test_id.txt
printf "gr1\ngr2\ngr1\n" > test_gr.txt
paste test_id.txt testR1.txt testR2.txt test_gr.txt > test_ed.txt
rm test_id.txt testR1.txt testR2.txt test_gr.txt
sed -i 's/:\t/:/g' test_ed.txt
echo -e "ID\tRawreads_files\tgroup" | cat - test_ed.txt > test_experimental_design.txt 
rm test_ed.txt


#
printf "running pipeline now"
rm -r pipeline_test
perl ../scripts/runPipeline_rRNA_noqsub_commandline.pl -test_kingdom prokaryote \
-significant_pvalue 0.001 -exp test_experimental_design.txt \
-d pipeline_test \
-prokaryote_fasta data/Bacillus_anthracis__Ames_Ancestor_uid58083.fa \
-eukarya_fasta data/cavPor3.fa -index_ref_bt2 test_index \
-gff_prokaryote data/Bacillus_anthracis__Ames_Ancestor_uid58083.gff \
-test_method EdgeR -gene_coverage_fasta data/Bacillus_anthracis__Ames_Ancestor_uid58083.fa



printf "running pipeline finished"