bio fasta genomes.gb --type gene --id N --end 10 > fasta_ids.fa
bio fasta genomes.gb -end 10 > fasta_all1.fa
cat genomes.gb | bio fasta --end  10 > fasta_all2.fa
diff fasta_all1.fa fasta_all2.fa > nodiff.txt
bio fasta genomes.gb --end 100 --source --rename "{isolate}" > fasta_rename.fa
bio fasta genomes.gb --end 100 --source --rename alias.txt > fasta_alias.fa
cat genomes.gb | bio fasta --olap 29514 -e 10 --type CDS > fasta_olap1.fa
bio fasta genomes.gb --end 10 --type CDS > fasta_cds.fa
bio fasta genomes.gb --type CDS --translate > fasta_translate.fa
bio fasta GATTACA --frame -3 --translate > fasta_frame.fa
bio fasta genomes.gb --protein > fasta_protein.fa
cat fasta_cds.fa | bio fasta -e -3 > fasta_start.fa
cat fasta_cds.fa | bio fasta -s -3 > fasta_stop.fa
bio fasta GATTACA GTTAACA GTTTATA GTTT > fasta_multi.fa
bio fasta --gene S --protein  genomes.gb > fasta_s.fa
bio align fasta_s.fa > align_s.txt
bio align fasta_s.fa --table > align_s.tsv
bio align fasta_s.fa --vcf > align_s.vcf
