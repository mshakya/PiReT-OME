
printf "running pipeline now"

printf "perl ../scripts/runPipeline_rRNA_noqsub_commandline.pl -test_kingdom both \
-significant_pvalue 0.001 -exp test_experimental_design.txt \
-d pipeline_test_both \
-prokaryote_fasta data/test_prok.fa \
-eukarya_fasta data/test_euk.fa -index_ref_bt2 test_index \
-gff_eukarya data/eukarya_test.gff3 -gff_prokaryote data/test_prok.gff \
-test_method both -gene_coverage_fasta data/test_prok.fa"

perl ../scripts/runPipeline_rRNA_noqsub_commandline.pl -test_kingdom both \
-significant_pvalue 0.001 -exp test_experimental_design.txt \
-d pipeline_test_both \
-prokaryote_fasta data/test_prok.fa \
-eukarya_fasta data/test_euk.fa -index_ref_bt2 test_index \
-gff_eukarya data/eukarya_test.gff3 -gff_prokaryote data/test_prok.gff \
-test_method both -gene_coverage_fasta data/test_prok.fa

printf "running pipeline finished"