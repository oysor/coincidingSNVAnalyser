import sys
import gzip
import re
import copy
from collections import OrderedDict
import cyvcf
import collections
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqUtils import nt_search
from Bio.Alphabet import generic_dna
import buildCSV as run

file1 = sys.argv[1]
database = sys.argv[1].split('.')[0]
chrom = '1'

oneKG_keys = ['EAS_AF_1KG',
              'EUR_AF_1KG',
              'AMR_AF_1KG',
              'SAS_AF_1KG',
              'AFR_AF_1KG']

exac_keys = ['AF_EAS_EXAC',
             'AF_NFE_EXAC',
             'AF_AMR_EXAC',
             'AF_SAS_EXAC',
             'AF_AFR_EXAC',
             'AF_Adj_EXAC',
             'AF_FIN_EXAC',
             'AF_OTH_EXAC']

dbSNP_keys = ['DBSNP_VALIDATION',
              'DBSNP_SUBMISSIONS']

cosmic_keys = ['COSMIC_CANCER_TYPE_GW',
               'COSMIC_CANCER_TYPE_ALL',
               'COSMIC_SITE_HISTOLOGY']

icgc_keys = ['ICGC_DONORS',
             'ICGC_PROJECTS']

## COSMIC ExAC
if database == "cosmic_exac":    
    scan_coinciding = run.COSMIC_ExAC()
    scan_unique = run.COSMIC()
    germline_keys = exac_keys
    somatic_keys = cosmic_keys
    unique_path = "unique_somatic/"
    coinciding_path = "coinciding_somatic/"
    path = "cosmic_exac/"
elif database == "exac_cosmic":
    scan_coinciding = run.COSMIC_ExAC()
    scan_unique = run.ExAC()
    germline_keys = exac_keys
    somatic_keys = cosmic_keys
    unique_path = "unique_germline/"
    coinciding_path = "coinciding_germline/"
    path = "cosmic_exac/"
## COSMIC oneKG
elif database == "cosmic_oneKG":
    scan_coinciding = run.COSMIC_oneKG()
    scan_unique = run.COSMIC()
    germline_keys = oneKG_keys
    somatic_keys = cosmic_keys
    unique_path = "unique_somatic/"
    coinciding_path = "coinciding_somatic/"
    path = "cosmic_oneKG/"
elif database == "oneKG_cosmic": 
    scan_coinciding = run.COSMIC_oneKG()
    scan_unique = run.oneKG()
    germline_keys = oneKG_keys
    somatic_keys = cosmic_keys
    unique_path = "unique_germline/"
    coinciding_path = "coinciding_germline/"
    path = "cosmic_oneKG/"
## COSMIC dbSNP
elif database == "cosmic_dbsnp":
    scan_coinciding = run.COSMIC_dbSNP()
    scan_unique = run.COSMIC()
    germline_keys = dbSNP_keys
    somatic_keys = cosmic_keys
    unique_path = "unique_somatic/"
    coinciding_path = "coinciding_somatic/"
    path = "cosmic_dbsnp/"
elif database == "dbsnp_cosmic":
    scan_coinciding = run.COSMIC_dbSNP()
    scan_unique = run.dbSNP()
    germline_keys = dbSNP_keys
    somatic_keys = cosmic_keys
    unique_path = "unique_germline/"
    coinciding_path = "coinciding_germline/"
    path = "cosmic_dbsnp/"
## ICGC ExAC   
elif database == "icgc_exac":
    scan_coinciding = run.ExAC()
    scan_unique = run.ICGC()
    germline_keys = exac_keys
    somatic_keys = icgc_keys
    unique_path = "unique_somatic/"
    coinciding_path = "coinciding_somatic/"
    path = "icgc_exac/"
elif database == "exac_icgc":
    scan_coinciding = run.ExAC()
    scan_unique = run.ExAC()
    germline_keys = exac_keys
    somatic_keys = icgc_keys
    unique_path = "unique_germline/"
    coinciding_path = "coinciding_germline/"
    path = "icgc_exac/"
## ICGC oneKG   
elif database == "icgc_oneKG":
    scan_coinciding = run.oneKG()
    scan_unique = run.ICGC()
    germline_keys = oneKG_keys
    somatic_keys = icgc_keys
    unique_path = "unique_somatic/"
    coinciding_path = "coinciding_somatic/"
    path = "icgc_oneKG/"
elif database == "oneKG_icgc":
    scan_coinciding = run.oneKG()
    scan_unique = run.oneKG()
    germline_keys = oneKG_keys
    somatic_keys = icgc_keys
    unique_path = "unique_germline/"
    coinciding_path = "coinciding_germline/"
    path = "icgc_oneKG/"
## ICGC dbSNP  
elif database == "icgc_dbsnp":
    scan_coinciding = run.dbSNP()
    scan_unique = run.ICGC()
    germline_keys = dbSNP_keys
    somatic_keys = icgc_keys
    unique_path = "unique_somatic/"
    coinciding_path = "coinciding_somatic/"
    path = "icgc_dbsnp/"
elif database == "dbsnp_icgc":
    scan_coinciding = run.dbSNP()
    scan_unique = run.dbSNP()
    germline_keys = dbSNP_keys
    somatic_keys = icgc_keys
    unique_path = "unique_germline/"
    coinciding_path = "coinciding_germline/"
    path = "icgc_dbsnp/"
else:
    print "wrong input file. Must either <database1>_<database2>.vcf.gz or <database2>_<database1>.vcf.gz"
    sys.exit(0)


## Reads the VCF header and stores the information in INFO_dict and VEP_dict
VEP_dict = OrderedDict()
vcf_header = gzip.open(file1, 'r')
for line in vcf_header:
    # Gathers all the info-CSQ VEP tags and stores them as keys in a dictionary  
    strings = line.split(',')
    if strings[0] == '##INFO=<ID=CSQ':
        for key in strings[3].split(':')[1].replace("\">\n","").strip().split('|'):
            VEP_dict[key] = None

    if line[0:2] != '##':
        break;        
vcf_header.close()

def is_snv(record):
    
    acgt = re.compile(r"[ACGT]{1}", re.IGNORECASE)
    
    if len(record.REF) == 1 and len(record.ALT[0]) == 1:
        if re.match(acgt, str(record.REF)) and re.match(acgt, record.ALT[0]):
            return True
    return False

def getReferencSequence(chrom):
    return SeqIO.read('chromFA/chr'+chrom+'.fa', "fasta")


def writeSummary(chrom,i,snv_count,coinciding_count,unique_count,not_snv,coinciding_skipped,unique_skipped,skipped):
    summary = open((path+"summary_"+database+"_"+unique_path[:-1]+".txt"), 'a+')
    summary.write(("####  chromosome "+chrom+" ####\n"))
    summary.write("parsed records: " + str(i)+"\n")
    summary.write("SNVs:" + str(snv_count)+"\n")
    summary.write("coinciding records: " +str(coinciding_count)+ "\n")
    summary.write("unique records: "+ str(unique_count)+"\n")
    summary.write("skipped records: "+ str((skipped + unique_skipped + coinciding_skipped + not_snv))+" ##\n")
    summary.write("-not SNVs:" + str(not_snv)+"\n")
    summary.write("-coinciding records skipped: " +str(coinciding_skipped)+ "\n")
    summary.write("-unique records skipped: "+ str(unique_skipped)+"\n")
    summary.close()

    
vcf_reader = cyvcf.Reader(open(file1,'r'))
ch_ref = getReferencSequence(chrom)

i = 0
snv_count = 0
not_snv = 0
coinciding_count = 0
coinciding_skipped = 0
unique_count = 0
unique_skipped = 0
skipped = 0



for record in vcf_reader:

    if record.CHROM != chrom:
        print "#### ",chrom," ####"
        print "variants parsed", i
        print "coinciding records:", coinciding_count
        print "unique records:", unique_count
        print "skipped records:", (skipped + unique_skipped + coinciding_skipped + not_snv)
        writeSummary(chrom,i,snv_count,coinciding_count,unique_count,not_snv,coinciding_skipped,unique_skipped,skipped)
        chrom = record.CHROM
        ch_ref = getReferencSequence(chrom)

    if 'CSQ' not in record.INFO:
        print record

    if is_snv(record) and 'CSQ' in record.INFO:
        snv_count += 1
        record_VEP = OrderedDict(zip(VEP_dict.keys(), record.INFO['CSQ'].split(',')))
        somatic_variant = any(key in record.INFO for key in somatic_keys)
        germline_variant = any(key in record.INFO for key in germline_keys)
        coinciding_variant = (somatic_variant and germline_variant)

        if unique_path == "unique_somatic/":
            unique_variant = (somatic_variant and not germline_variant)
        else:
            unique_variant = (germline_variant and not somatic_variant)

        if coinciding_variant:  
            if scan_coinciding.regVariant(record, record_VEP, ch_ref):
                coinciding_count += 1
            else:
                coinciding_skipped += 1

        elif unique_variant:
            if scan_unique.regVariant(record, record_VEP, ch_ref):
                unique_count += 1
            else:
                unique_skipped += 1         
        else:
            skipped += 1            
    else:
        not_snv += 1
    

print "variants parsed", i
print "skipped records:", (skipped + unique_skipped + coinciding_skipped + not_snv)
print "coinciding records", coinciding_count
print "unique records", unique_count
writeSummary("Summary",i,snv_count,coinciding_count,unique_count,not_snv,coinciding_skipped,unique_skipped,skipped)

scan_coinciding.writeToFile(scan_coinciding.getVariantPlot(), path+coinciding_path+database+"_variantPlot.csv", "variantPlot","coinciding",scan_coinciding.getHeader())
scan_coinciding.writeToFile(scan_coinciding.getContextPlot(), path+coinciding_path+database+"_contextPlot.csv", "contextPlot","coinciding",scan_coinciding.getHeader())
scan_coinciding.writeToFile(scan_coinciding.getConsequencePlot(), path+coinciding_path+database+"_consequencePlot.csv", "consequencePlot","coinciding",scan_coinciding.getHeader())
scan_unique.writeToFile(scan_unique.getVariantPlot(), path+unique_path+database+"_variantPlot.csv", "variantPlot",unique_path[:-1],scan_unique.getHeader())
scan_unique.writeToFile(scan_unique.getContextPlot(), path+unique_path+database+"_contextPlot.csv", "contextPlot",unique_path[:-1],scan_unique.getHeader())
scan_unique.writeToFile(scan_unique.getConsequencePlot(), path+unique_path+database+"_consequencePlot.csv", "consequencePlot",unique_path[:-1],scan_unique.getHeader())

if (database == "cosmic_exac" or database == "cosmic_oneKG" or database == "cosmic_dbSNP"):
    print len(scan_coinciding.variantTable)
    scan_coinciding.writeTable(scan_coinciding.variantTable, (path+"table/variantTable.csv"))
