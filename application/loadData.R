
## cosmic_1kg
cosmic_1kg = list(
  "coinciding" = fread("data/coinciding/COSMIC_oneKG_coinciding_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "germline" = fread("data/germline/COSMIC_oneKG_unique_germline_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "somatic"  = fread("data/somatic/COSMIC_oneKG_unique_somatic_genome.csv", header=T, sep=",", stringsAsFactors=F)
)
## COSMIC_ExAC
cosmic_ExAC = list(
  "coinciding" = fread("data/coinciding/COSMIC_ExAC_coinciding_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "germline" = fread("data/germline/COSMIC_ExAC_unique_germline_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "somatic" = fread("data/somatic/COSMIC_ExAC_unique_somatic_genome.csv", header=T, sep=",", stringsAsFactors=F)
)
## COSMIC_dbSNP
cosmic_dbSNP = list(
  "coinciding" = fread("data/coinciding/COSMIC_dbSNP_coinciding_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "germline" = fread("data/germline/COSMIC_dbSNP_unique_germline_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "somatic" = fread("data/somatic/COSMIC_dbSNP_unique_somatic_genome.csv", header=T, sep=",", stringsAsFactors=F)
)
## icgc_1kg
icgc_1kg = list(
  "coinciding"  = fread("data/coinciding/ICGC_oneKG_coinciding_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "germline"  = fread("data/germline/ICGC_oneKG_unique_germline_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "somatic" = fread("data/somatic/ICGC_oneKG_unique_somatic_genome.csv", header=T, sep=",", stringsAsFactors=F)
)
## icgc_ExAC
icgc_ExAC = list(
  "coinciding" = fread("data/coinciding/ICGC_ExAC_coinciding_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "germline" = fread("data/germline/ICGC_ExAC_unique_germline_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "somatic" = fread("data/somatic/ICGC_ExAC_unique_somatic_genome.csv", header=T, sep=",", stringsAsFactors=F)
)
## icgc_dbSNP
icgc_dbSNP = list(
  "coinciding" = fread("data/coinciding/ICGC_dbSNP_coinciding_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "germline" = fread("data/germline/ICGC_dbSNP_unique_germline_genome.csv", header=T, sep=",", stringsAsFactors=F),
  "somatic" = fread("data/somatic/ICGC_dbSNP_unique_somatic_genome.csv", header=T, sep=",", stringsAsFactors=F)
)
### coinciding variant tables
cosmic_1kg_table <- fread("data/tables/cosmic_oneKG_coinciding_variants_updated.csv", header=T, sep=",", stringsAsFactors=F)
cosmic_ExAC_table <- fread("data/tables/cosmic_exac_coinciding_variants_updated.csv", header=T, sep=",", stringsAsFactors=F)
cosmic_dbSNP_table <- fread("data/tables/cosmic_dbSNP_coinciding_variants_updated.csv", header=T, sep=",", stringsAsFactors=F)

cosmic_1kg_table_DL <- fread("data/tables/cosmic_oneKG_coinciding_variants.csv", header=T, sep=",", stringsAsFactors=F)
cosmic_ExAC_table_DL <- fread("data/tables/cosmic_exac_coinciding_variants.csv", header=T, sep=",", stringsAsFactors=F)
cosmic_dbSNP_table_DL <- fread("data/tables/cosmic_dbSNP_coinciding_variants.csv", header=T, sep=",", stringsAsFactors=F)

variantTables <- list(
  "cosmic_oneKG" = cosmic_1kg_table ,
  "cosmic_ExAC" = cosmic_ExAC_table,
  "cosmic_dbSNP" = cosmic_dbSNP_table,
  "cosmic_oneKG_DL" = cosmic_1kg_table_DL,
  "cosmic_ExAC_DL" = cosmic_ExAC_table_DL,
  "cosmic_dbSNP_DL" = cosmic_dbSNP_table_DL
)

databases = list("cosmic_ExAC"=cosmic_ExAC,
                 "cosmic_oneKG"=cosmic_1kg,
                 "cosmic_dbSNP"=cosmic_dbSNP,
                 "icgc_ExAC"=icgc_ExAC,
                 "icgc_oneKG"=icgc_1kg,
                 "icgc_dbSNP"=icgc_dbSNP
)

signatureTable <- fread("data/signatures_aetiologies.tsv", header=T, sep="\t", stringsAsFactors=F)



