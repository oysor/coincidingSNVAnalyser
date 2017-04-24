from collections import OrderedDict
import copy
import re
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from Bio.SeqUtils import nt_search
from Bio.Alphabet import generic_dna
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

class MutationType(object):

    def __init__(self):
        self.mType = OrderedDict()
        # Transitions
        self.mType["A>G:T>C"] = 0
        self.mType["C>T:G>A"] = 0
        # Tranversions
        self.mType["C>A:G>T"] = 0
        self.mType["C>G:G>C"] = 0
        self.mType["A>T:T>A"] = 0
        self.mType["A>C:T>G"] = 0
    
    def getDict(self):
        return self.mType
    
    def getType(self, record):
        REF = record.REF.upper()
        ALT = record.ALT[0].upper()
        
        mutation = REF+">"+ALT
        complement = Seq(mutation).complement()

        if (mutation+":"+complement) in self.mType:
            return (mutation+":"+complement)
        else:
            return (complement+":"+mutation)       

class MutationContext(object):

    def __init__(self):
        self.mutContext = OrderedDict()
        # including the context
        for mut in ["C>A","C>G", "C>T", "T>A","T>C","T>G"]:
            for b1 in "ACGT":
                for b2 in "ACGT":
                    self.mutContext[(b1+"["+mut+"]"+b2+":"+Seq((b1+"["+mut+"]"+b2)).complement())] = 0

    def getDict(self):
        return self.mutContext

    def getMutContext(self, seq, record):
        POS = record.POS
        REF = record.REF
        ALT = record.ALT[0]

        if seq[POS-2].upper() == 'N'or seq[POS].upper() == 'N':
            return None    
        
        context = seq[POS-2].upper()+"["+REF+">"+ALT+"]"+seq[POS].upper()
        complement =  Seq(context).complement()

        if (context+":"+complement) in self.mutContext:
            return (context+":"+complement)
        else:
            return (complement+":"+context)

#### data found in 1000Genomes and ExAC ####
class Population(object):

    def __init__(self, population_db):
        
        self.popFreq = OrderedDict()
        self.popFreq['Common'] = 0
        self.popFreq['LowFreq'] = 0
        self.popFreq['Rare'] = 0
        self.popFreq['VeryRare'] = 0

        self.population = OrderedDict()
        pop_1kg = ['EAS_AF_1KG',
                   'EUR_AF_1KG',
                   'AMR_AF_1KG',
                   'SAS_AF_1KG',
                   'AFR_AF_1KG']
        pop_ExAC = ['AF_EAS_EXAC',
                    'AF_NFE_EXAC',
                    'AF_AMR_EXAC',
                    'AF_SAS_EXAC',
                    'AF_AFR_EXAC',
                    'AF_FIN_EXAC',
                    'AF_OTH_EXAC',
                    'AF_Adj_EXAC',]
        
        if population_db == '1kg':
            for p in pop_1kg:
                 self.population[p] = copy.deepcopy(self.popFreq)
        elif population_db == 'ExAC':
             for p in pop_ExAC:
                 self.population[p] = copy.deepcopy(self.popFreq)

    def getDict(self):
        return self.population

    # returns a dictionairy of one or more populations found in this sample.
    # The frequency is based on minor allele frequency
    def getPopulation(self, record):   
        pop = {}       
        for k in record.INFO.keys():
            if k in self.population:
                MAF = record.INFO[k]
                if MAF > 0.0:
                    if MAF >= 0.05: 
                        pop[k] = 'Common'
                    elif (MAF < 0.05 and MAF >= 0.01):
                        pop[k] = 'LowFreq'
                    elif (MAF < 0.01 and MAF >= 0.001):
                        pop[k] = 'Rare'
                    elif (MAF < 0.001):
                        pop[k] = 'VeryRare'                
        return pop

#### data found in dbSNP  ####    
class Validate(object):

    def __init__(self):
        self.valDict = OrderedDict()
        self.valDict["byCluster"] = 0
        self.valDict["byFrequency"] = 0
        self.valDict["by1000G"] = 0
        self.valDict["byOtherPop"] = 0
        self.valDict["suspect"] = 0
        self.valDict["byHapMap"] = 0
        self.valDict["by2Hit2Allele"] = 0
        self.valDict["multiple"] = 0
        self.valDict["unknown"] = 0
     #   self.valDict["missing"] = 0
   
    def getDict(self):
        return self.valDict

     # return a list with all the dbSNP validations registered for this sample
    def getValidation(self, record):
        validation = []
        if 'DBSNP_VALIDATION' in record.INFO.keys():
            dbSNP_validation = record.INFO['DBSNP_VALIDATION']
        else:
            return []
        
        for v in dbSNP_validation.split(','):
            validation.append(v)

        if len(validation) > 1:
            return "multiple"
        else:
            return validation[0] 


####  data found in dbSNP ####
class dbSNPsubmission(object):
    
    def __init__(self):
        self.submDict = OrderedDict()
        self.submDict['1'] = 0
        self.submDict['2:5'] = 0
        self.submDict['6:10'] = 0
        self.submDict['>10'] = 0
 
    def getDict(self):
        return self.submDict
    
    def getSubmission(self, record):

        subm = 0
        if 'DBSNP_SUBMISSIONS' not in record.INFO.keys():
            return None
        else:
            subm = record.INFO['DBSNP_SUBMISSIONS'][0]

        if subm == 1:  return '1' 
        elif subm > 1 and subm <= 5:  return '2:5'
        elif subm > 5 and subm <= 10:  return '6:10'
        elif subm > 10:  return '>10'


#### data found in COSMIC ####        
class CancerType(object):

    def __init__(self):
        self.cancerDict = OrderedDict()

        ####  cancer types that we choose to register  ####
        self.cancerTypes = [
                      #"leukemia",
                      "chronic_lymphocytic_leukemia",
                      "acute_myeloid_leukemia",
                      "acute_lymphoblastic_b-cell_leukemia",
                      "acute_lymphoblastic_leukemia",
                      "lung_cancer", 
                      "breast_cancer", 
                      "malignant_melanoma",
                      "sarcoma", 
                      "stomach_cancer",
                      "prostate_cancer", 
                      "diffuse_large_B_cell_lymphoma", 
                      "glioma", 
                      "colorectal_cancer", 
                      "oesophageal_cancer",
                      "ovarian_cancer",
                      "pancreatic_cancer", 
                      "liver_cancer", 
                      #"bladder_cancer", 
                      "cervical_cancer", 
                      "cholangiocarcinoma", 
                      "urothelical_cancer",
                      "kidney_cancer", 
                      "pancancer",
                     ]
        
        for t in self.cancerTypes:
            cancerFreq = OrderedDict()
            cancerFreq['recurrent'] = 0
            cancerFreq['nonrecurrent'] = 0
            self.cancerDict[t] = cancerFreq

    def getDict(self):
        return self.cancerDict

    ### returns a dictionairy of one or more cancers from this sample
    def getCancer(self, record):

        cancer = {}
        if 'COSMIC_CANCER_TYPE_GW' not in record.INFO:
            return {}  

        for cancerType in record.INFO['COSMIC_CANCER_TYPE_GW'].split('&'):
            c = cancerType.split(':')
          
            # c[0] = cancer type, c[1] = cancer frequency
            if c[0] in self.cancerTypes:  
                if int(c[1]) > 1:
                    cancer[c[0]] = 'recurrent'  
                else:
                    cancer[c[0]] = 'nonrecurrent'


        if len(cancer.keys()) == 1 and  "pancancer" in cancer.keys(): 
            return {}
        else:
            return cancer

#### Information found in the CSQ tag added by Variant Effect Predictor (VEP).
class Consequence(object):

    def __init__(self):
        
        self.consequenceDict = OrderedDict()
        ###  http://www.ensembl.org/info/genome/variation/predicted_data.html#consequences
        self.consequenceList = ['splice_acceptor_variant',
                                'splice_donor_variant',
                                'stop_gained',
                                #'frameshift_variant',
                                #'inframe_insertion',
                                #'inframe_deletion',
                                'stop_lost',
                                'stop_retained_variant',
                                'start_lost',
                                'missense_variant',
                                'synonymous_variant',
                                #'transcript_ablation',
                                #'transcript_amplification',
                                #'protein_altering_variant',
                                'splice_region_variant',
                                #'incomplete_terminal_codon_variant',
                                #'coding_sequence_variant',
                                #'mature_miRNA_variant',
                                '5_prime_UTR_variant',
                                '3_prime_UTR_variant',
                                'non_coding_transcript_exon_variant',
                                'intron_variant',
                                'NMD_transcript_variant',
                                'non_coding_transcript_variant',
                                'upstream_gene_variant',
                                'downstream_gene_variant',
                                #'TFBS_ablation',
                                #'TFBS_amplification',
                                #'TF_binding_site_variant',
                                #'regulatory_region_ablation',
                                #'regulatory_region_amplification',
                                #'feature_elongation',
                                #'regulatory_region_variant',
                                #'feature_truncation',
                                'intergenic_variant']

        self.codingList = [
            'splice_acceptor_variant',
            'splice_donor_variant',
            'stop_gained',
            'stop_retained_variant',
            'stop_lost',
            'start_lost',
            'missense_variant',
            'synonymous_variant'
        ]

        self.nonCodingList = [
            'splice_region_variant',
            '5_prime_UTR_variant',
            '3_prime_UTR_variant',
            'non_coding_transcript_exon_variant',
            'intron_variant',
            'NMD_transcript_variant',
            'non_coding_transcript_variant',
            'upstream_gene_variant',
            'downstream_gene_variant',
            'intergenic_variant'
        ]

        for c in self.consequenceList:
            consCoding = OrderedDict()
            if c in self.codingList:
                consCoding['coding'] = 0
            else:
                consCoding['noncoding'] = 0  
            self.consequenceDict[c] = consCoding
        
        self.codingDict = OrderedDict()
        self.codingDict['coding'] = 0
        self.codingDict['noncoding'] = 0
              
    def getDict(self):
        return self.consequenceDict

    def getCodingDict(self):
        return self.codingDict
           
    ### Returns only the first of potentially several consequences for each sample
    def getConsequence(self, vep):     
        consequence = vep['Consequence'].split('&')[0] 
        if consequence in self.consequenceList:
            return consequence
        else:
            None

    ### Returns only the coding of the first consequence type   
    def getCodingType(self, vep):

        consequence = vep['Consequence'].split('&')[0]
        if consequence in self.consequenceList:
            return  self.consequenceDict[consequence].keys()[0]
        else:
            return None        

    def getGene(self, vep):
        return vep['Gene']
        
    def getBiotype(self, vep):
        return vep['BIOTYPE']
