# coincidingSNVAnalyser

Computational analysis of coinciding single nucleotide variants (SNVs) in germline and somatic mutation spectra.

https://sigven78.shinyapps.io/thesisApp/

In short this application is a tool for analysing and comparing the mutation spectra in normal DNA variation and DNA variation related to cancer. In order to do this, the application use processed data from five international sequencing consortias: Exome Aggregation Consortium (ExAC), 1000 Genomes Project (1000Genomes), database of short genetic variants (dbSNP), the International Cancer Genomics Consortium (ICGC) and Catalogue of Somatic Mutations in Cancer (COSMIC).

The application displays up to six different comparisons, COSMIC & 1000Genomes, COSMIC & ExAC, COSMIC & dbSNP, ICGC & 1000Genomes, ICGC & ExAC and ICGC & dbSNP. Furthermore, variants of each comparison are categorized into groups, first those that only occur in germline, second, those that occur in somatic and third, those that occur in both germline and somatic (shared). These three sets of variant data are referred to as unique somatic variants, unique germline variants and coinciding variants respectively.


/application

1. Put coinciding variant tables ([CSV files](https://drive.google.com/drive/folders/0B6GfJ6vekOM9QnJSRFVDVmZyODA?usp=sharing)) in folder application/data/tables/

/preprocessing

To quantify and build CSV files out of VCF files:
1. Go to folder /preprocessing | cd preprocessing
2. Download intersected VCF ([12 files](https://drive.google.com/drive/folders/0B6GfJ6vekOM9SVU4TlJvbzRQYms?usp=sharing))  
3. Download reference sequence (GRCh37, .fa format) and place in folder /chromFa  
4. Make folders (cosmic_exac, cosmic_oneKG, cosmic_dbsnp, icgc_exac, icgc_oneKG, icgc_dbsnp) | $ bash setup.sh
5. Run python script for each intersected VCF file | $ python quantify.py cosmic_exac.vcf.py 






