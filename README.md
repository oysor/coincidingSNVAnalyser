# coincidingSNVAnalyzer

Computational analysis of coinciding single nucleotide variants (SNVs) in germline and somatic mutation spectra.

https://sigven78.shinyapps.io/coincidingSNVAnalyser/

In short this application is a tool for the analysis and comparison of the mutational properties of coinciding DNA variation in the germline and the soma (i.e. cancer). The application employs processed data from five international sequencing consortias: Exome Aggregation Consortium (ExAC), 1000 Genomes Project (1000Genomes), database of short genetic variants (dbSNP), the International Cancer Genomics Consortium (ICGC) and Catalogue of Somatic Mutations in Cancer (COSMIC).

The application displays six different data comparisons, COSMIC versus 1000Genomes, COSMIC versus ExAC, COSMIC versus dbSNP, ICGC versus 1000Genomes, ICGC versus ExAC and ICGC versus dbSNP. Furthermore, variants within each comparison set are categorized into groups, first those that only occur as germline, second, those that occur only as somatic and third, those that occur both as germline and somatic (i.e. shared or coinciding). These three groups of variant datasets are hereby referred to as unique somatic variants, unique germline variants and coinciding variants respectively. 


/application

1. Put coinciding variant tables ([CSV files](https://drive.google.com/drive/folders/0B6GfJ6vekOM9QnJSRFVDVmZyODA?usp=sharing)) in folder application/data/tables/

/preprocessing

To quantify and build CSV files out of VCF files:
1. Go to folder /preprocessing | $ cd preprocessing
2. Download intersected VCF ([12 files](https://drive.google.com/drive/folders/0B6GfJ6vekOM9SVU4TlJvbzRQYms?usp=sharing))  
3. Download reference sequence (GRCh37, .fa format) and extract file in folder /chromFa  
4. Make folders (cosmic_exac, cosmic_oneKG, cosmic_dbsnp, icgc_exac, icgc_oneKG, icgc_dbsnp) | $ bash folderSetup.sh
5. Run python script for each intersected VCF | $ python quantify.py cosmic_exac.vcf.py 
6. Script output CSV files in respective folders.





