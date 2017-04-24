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

import parseRecord as an

class WriteCSV(object):

    def writeToFile(self, fileDict, theFile, plot, dset, csv_header):
        with open(theFile, 'w') as f:
            f.write(csv_header)
            for k,n in fileDict.items():
                line = "%s,%s,%s,%d" % (plot,dset,k,n)
                f.write(line+"\n")
        f.close()
        print "Done: ", theFile

    def writeTable(self, table, outputFile):
    
        csv_header = "assembly,gdna_pos,consequence,gene,symbol,biotype,cosmic_mutation_id,cancer,info\n"
        with open(outputFile, 'w+') as f:
            f.write(csv_header)
            for line in table:
                f.write(line+"\n")
        f.close()
        print "Done: ", outputFile


        
class COSMIC(WriteCSV):
    
    def __init__(self):
        self.mt = an.MutationType()
        self.mc = an.MutationContext()
        self.cq = an.Consequence()
        self.ct = an.CancerType()
           
        self.csv_header = "plot,set,cancer,cancerFreq,coding,type,n\n"
        self.variantPlot = OrderedDict()
        self.contextPlot = OrderedDict()
        self.consequencePlot = OrderedDict()
        
        for cancer, cancer_freq_list in self.ct.getDict().items():
            for cancerFreq in cancer_freq_list.keys():
                for coding in self.cq.getCodingDict():
                        
                    for mutation in self.mt.getDict().keys():
                        self.variantPlot[cancer+","+cancerFreq+","+coding+","+mutation] = 0
                                 
                    for context in self.mc.getDict().keys():
                        self.contextPlot[cancer+","+cancerFreq+","+coding+","+context] = 0
                                 
                for consequence, coding  in self.cq.getDict().items():
                    self.consequencePlot[cancer+","+cancerFreq+","+coding.keys()[0]+","+consequence] = 0
                                 
    def regVariant(self, record, recordVEP, referenceSequence):
        
        mutation = self.mt.getType(record)
        context  = self.mc.getMutContext(referenceSequence, record)
        consequence = self.cq.getConsequence(recordVEP)
        coding = self.cq.getCodingType(recordVEP)
        cancer = self.ct.getCancer(record)
        
        if consequence: 
            if context: 
                if cancer: 
                    for cancerType, cancerFreq in  cancer.items():
                        self.variantPlot[cancerType +","+ cancerFreq  +","+ coding +","+mutation] += 1
                        self.contextPlot[cancerType +","+ cancerFreq  +","+ coding +","+ context] += 1
                        self.consequencePlot[cancerType +","+ cancerFreq  +","+ coding +","+ consequence] += 1
                        return True

    
    def getVariantPlot(self):
        return self.variantPlot
        
    def getContextPlot(self):
        return self.contextPlot
        
    def getConsequencePlot(self):
        return self.consequencePlot
    
    def getHeader(self):
        return self.csv_header
    

class ICGC(WriteCSV):

    def __init__(self):
        self.mt = an.MutationType()
        self.mc = an.MutationContext()
        self.cq = an.Consequence()
        self.csv_header = "plot,set,coding,type,n\n"
        self.variantPlot = OrderedDict()
        self.contextPlot = OrderedDict()
        self.consequencePlot = OrderedDict()
      
        for coding in self.cq.getCodingDict():
            
            for mutation in self.mt.getDict().keys():
                self.variantPlot[coding+","+mutation] = 0
                    
            for context in self.mc.getDict().keys():
                self.contextPlot[coding+","+context] = 0
                                 
        for consequence, coding in self.cq.getDict().items():
            self.consequencePlot[coding.keys()[0]+","+consequence] = 0
                                 
    def regVariant(self, record, recordVEP, referenceSequence):

        mutation = self.mt.getType(record)
        context  = self.mc.getMutContext(referenceSequence, record)
        consequence = self.cq.getConsequence(recordVEP)
        coding = self.cq.getCodingType(recordVEP)
        
        if consequence: 
            if context:
                self.variantPlot[coding +","+ mutation] += 1
                self.contextPlot[coding +","+ context] += 1
                self.consequencePlot[coding +","+ consequence] += 1
                return True
        return False
    
    def getVariantPlot(self):
        return self.variantPlot
        
    def getContextPlot(self):
        return self.contextPlot
        
    def getConsequencePlot(self):
        return self.consequencePlot
    
    def getHeader(self):
        return self.csv_header
    
        
class ExAC(WriteCSV):

    def __init__(self):
        self.mt = an.MutationType()
        self.mc = an.MutationContext()
        self.cq = an.Consequence()
        self.pc = an.Population('ExAC')
        self.csv_header = "plot,set,region,alleleFreq,coding,type,n\n"
        self.variantPlot = OrderedDict()
        self.contextPlot = OrderedDict()
        self.consequencePlot = OrderedDict()
      
        for pop, pop_freq_list in self.pc.getDict().items():
            for pop_freq in pop_freq_list.keys():
                for coding in self.cq.getCodingDict():
                        
                    for mutation in self.mt.getDict().keys():
                        self.variantPlot[pop+","+pop_freq+","+coding+","+mutation] = 0
                            
                    for context in self.mc.getDict().keys():
                        self.contextPlot[pop+","+pop_freq+","+coding+","+context] = 0
                                 
                for consequence, coding  in self.cq.getDict().items():
                    self.consequencePlot[pop+","+pop_freq+","+coding.keys()[0]+","+consequence] = 0
                                 
    def regVariant(self, record, recordVEP, referenceSequence):
        
        mutation = self.mt.getType(record)
        context  = self.mc.getMutContext(referenceSequence, record)
        consequence = self.cq.getConsequence(recordVEP)
        coding = self.cq.getCodingType(recordVEP)
        population = self.pc.getPopulation(record)   
        
        if consequence:
            if context:
                if population:
                    for pop, pop_freq in  population.items():
                        self.variantPlot[pop +","+ pop_freq  +","+ coding +","+ mutation] += 1
                        self.contextPlot[pop +","+ pop_freq  +","+ coding +","+ context] += 1
                        self.consequencePlot[pop +","+ pop_freq  +","+ coding +","+ consequence] += 1
                    return True       
        return False
    
    def getVariantPlot(self):
        return self.variantPlot
        
    def getContextPlot(self):
        return self.contextPlot
        
    def getConsequencePlot(self):
        return self.consequencePlot
    
    def getHeader(self):
        return self.csv_header

class oneKG(WriteCSV):

    def __init__(self):
        self.mt = an.MutationType()
        self.mc = an.MutationContext()
        self.cq = an.Consequence()
        self.pc = an.Population('1kg')
        self.csv_header = "plot,set,region,alleleFreq,coding,type,n\n"
        self.variantPlot = OrderedDict()
        self.contextPlot = OrderedDict()
        self.consequencePlot = OrderedDict()
       
        for pop, pop_freq_list in self.pc.getDict().items():
            for pop_freq in pop_freq_list.keys():
                for coding in self.cq.getCodingDict():

                    for mutation in self.mt.getDict().keys():
                        self.variantPlot[pop+","+pop_freq+","+coding+","+mutation] = 0
                                 
                    for context in self.mc.getDict().keys():
                        self.contextPlot[pop+","+pop_freq+","+coding+","+context] = 0
                                 
                for consequence, coding  in self.cq.getDict().items():
                    self.consequencePlot[pop+","+pop_freq+","+coding.keys()[0]+","+consequence] = 0
                                 
    def regVariant(self, record, recordVEP, referenceSequence):
        
        mutation = self.mt.getType(record)
        context  = self.mc.getMutContext(referenceSequence, record)
        consequence = self.cq.getConsequence(recordVEP)
        coding = self.cq.getCodingType(recordVEP)
        population = self.pc.getPopulation(record)

        if consequence: 
            if context:
                if population:
                    for pop, pop_freq in  population.items():
                        self.variantPlot[pop +","+ pop_freq  +","+ coding +","+ mutation] += 1
                        self.contextPlot[pop +","+ pop_freq  +","+ coding +","+ context] += 1
                        self.consequencePlot[pop +","+ pop_freq  +","+ coding +","+ consequence] += 1
                    return True
        return False

                    
    def getVariantPlot(self):
        return self.variantPlot
        
    def getContextPlot(self):
        return self.contextPlot
        
    def getConsequencePlot(self):
        return self.consequencePlot
    
    def getHeader(self):
        return self.csv_header
    
class dbSNP(WriteCSV):

    def __init__(self):
        self.mt = an.MutationType()
        self.mc = an.MutationContext()
        self.cq = an.Consequence()
        self.sm = an.dbSNPsubmission()
        self.va = an.Validate() 
        self.csv_header = "plot,set,validation,submission,coding,type,n\n"
        self.variantPlot = OrderedDict()
        self.contextPlot = OrderedDict()
        self.consequencePlot = OrderedDict()
        
        for val in self.va.getDict():
            for sub in self.sm.getDict():
                for coding in self.cq.getCodingDict():

                    for mutation in self.mt.getDict().keys():
                        self.variantPlot[val+","+sub+","+coding+","+mutation] = 0
                            
                    for context in self.mc.getDict().keys():
                        self.contextPlot[val+","+sub+","+coding+","+context] = 0
                                 
                for consequence, coding  in self.cq.getDict().items():
                    self.consequencePlot[val+","+sub+","+coding.keys()[0]+","+consequence] = 0
                                 
       
    def regVariant(self, record, recordVEP, referenceSequence):
        
        mutation = self.mt.getType(record)
        context  = self.mc.getMutContext(referenceSequence, record)
        consequence = self.cq.getConsequence(recordVEP)
        coding = self.cq.getCodingType(recordVEP)
        validation = self.va.getValidation(record)
        submission = self.sm.getSubmission(record)

        if consequence:
            if context: 
                if submission:
                    self.variantPlot[validation +","+ submission  +","+ coding +","+ mutation] += 1
                    self.contextPlot[validation +","+ submission  +","+ coding +","+ context] += 1
                    self.consequencePlot[validation +","+ submission  +","+ coding +","+ consequence] += 1
                    return True
        return False

                    
    def getVariantPlot(self):
        return self.variantPlot
        
    def getContextPlot(self):
        return self.contextPlot
        
    def getConsequencePlot(self):
        return self.consequencePlot
    
    def getHeader(self):
        return self.csv_header
        
class COSMIC_dbSNP(WriteCSV):

    def __init__(self):
        self.mt = an.MutationType()
        self.mc = an.MutationContext()
        self.cq = an.Consequence()
        self.ct = an.CancerType()
        self.vd = an.Validate()
        self.sm = an.dbSNPsubmission()
        self.csv_header = "plot,set,cancer,cancerFreq,validation,submission,coding,type,n\n"
        self.variantPlot = OrderedDict()
        self.contextPlot = OrderedDict()
        self.consequencePlot = OrderedDict()
        self.variantTable = []
        
        for cancer, cancer_freq_list in self.ct.getDict().items():
            for cancer_freq in cancer_freq_list.keys():
                for val in self.vd.getDict():
                    for sub in self.sm.getDict():
                        for coding in self.cq.getCodingDict():
                            
                            for mutation in self.mt.getDict().keys():
                                self.variantPlot[cancer+","+cancer_freq+","+val+","+sub+","+coding+","+mutation] = 0
                                 
                            for context in self.mc.getDict().keys():
                                self.contextPlot[cancer+","+cancer_freq+","+val+","+sub+","+coding+","+context] = 0
                                 
                        for consequence, coding  in self.cq.getDict().items():   
                            self.consequencePlot[cancer+","+cancer_freq+","+val+","+sub+","+coding.keys()[0]+","+consequence] = 0

    def regVariant(self, record, recordVEP, referenceSequence):
        
        mutation = self.mt.getType(record)
        context  = self.mc.getMutContext(referenceSequence, record)
        consequence = self.cq.getConsequence(recordVEP)
        coding = self.cq.getCodingType(recordVEP)
        
        cancer = self.ct.getCancer(record)
        coding = self.cq.getCodingType(recordVEP)
        validation = self.vd.getValidation(record)
        submission = self.sm.getSubmission(record)

        
        if  consequence: 
            if context: 
                if cancer: 
                    if submission:
                        if cancer:
                            for cancerType, cancerFreq in  cancer.items():
                                self.variantPlot[cancerType+","+ cancerFreq  +","+validation+","+ submission+ ","+coding+","+mutation] += 1
                                self.contextPlot[cancerType+","+ cancerFreq  +","+validation+","+ submission +","+coding+","+context] += 1
                                self.consequencePlot[cancerType+","+cancerFreq+","+validation+","+ submission +","+coding+","+consequence] += 1
                            self.putInVariantTable(record, recordVEP, context, validation, submission, cancer, coding)
                            return True
        return False


    def putInVariantTable(self, record, record_VEP, context, validation, submission, cancer, coding):
        gdnaPos = record.CHROM+":"+str(record.POS)+":"+record.REF+">"+record.ALT[0]
        variant = "GRCh37,"
        variant += gdnaPos+","
        variant += record_VEP['Consequence']+","
        variant += record_VEP['Gene']+","
        variant += record_VEP['SYMBOL']+","
        variant += record_VEP["BIOTYPE"]+","
        variant += ''.join(record.INFO["COSMIC_MUTATION_ID"].split(','))+","
        cancerList = []
        for cancerType in record.INFO['COSMIC_CANCER_TYPE_GW'].split('&'):
            c = cancerType.split(':')
            if int(c[1]) > 1:
                cancerFreq = "(recurrent)"
            else:
                cancerFreq = "(nonrecurrent)"
            if c[0] in cancer:   
                cancerList.append(cancerType+cancerFreq)

        variant += "&".join(cancerList)
        variant += ",variant consequence="+coding
        if len(record.INFO['DBSNP_VALIDATION'].split(',')) > 1:
            variant += ";validation method="+'&'.join(record.INFO['DBSNP_VALIDATION'].split(','))+"("+validation+")"
        else:
            variant += ";validation method="+record.INFO['DBSNP_VALIDATION'][0]
        variant += ";submissions="+str(record.INFO['DBSNP_SUBMISSIONS'][0])+"("+submission+")"

        self.variantTable.append(variant)
    
    def getVariantPlot(self):
        return self.variantPlot
        
    def getContextPlot(self):
        return self.contextPlot
        
    def getConsequencePlot(self):
        return self.consequencePlot
    
    def getHeader(self):
        return self.csv_header

    
class COSMIC_ExAC(WriteCSV):

    def __init__(self):
        self.mt = an.MutationType()
        self.mc = an.MutationContext()
        self.cq = an.Consequence()
        self.ct = an.CancerType()
        self.pc = an.Population('ExAC')
        self.csv_header = "plot,set,region,alleleFreq,cancer,cancerFreq,coding,type,n\n"
        self.variantPlot = OrderedDict()
        self.contextPlot = OrderedDict()
        self.consequencePlot = OrderedDict()
        self.variantTable = []
        
        for region, pop_freq_list in self.pc.getDict().items():
            for regionFreq in pop_freq_list.keys():
                for cancer, cancer_freq_list in self.ct.getDict().items():
                    for cancerFreq in cancer_freq_list.keys():
                        for coding in self.cq.getCodingDict():
                             
                             for mutation in self.mt.getDict().keys():
                                 self.variantPlot[region+","+regionFreq+","+cancer+","+cancerFreq+","+coding+","+mutation] = 0
                             for context in self.mc.getDict().keys():
                                 self.contextPlot[region+","+regionFreq+","+cancer+","+cancerFreq+","+coding+","+context] = 0
                                 
                        for consequence, coding  in self.cq.getDict().items(): 
                            self.consequencePlot[region+","+regionFreq+","+cancer+","+cancerFreq+","+coding.keys()[0]+","+consequence] = 0                            

    def regVariant(self, record, recordVEP, referenceSequence):
         
        mutation = self.mt.getType(record)
        context  = self.mc.getMutContext(referenceSequence, record)
        consequence = self.cq.getConsequence(recordVEP)
        coding = self.cq.getCodingType(recordVEP)
        cancer = self.ct.getCancer(record)
        coding = self.cq.getCodingType(recordVEP)
        population = self.pc.getPopulation(record)
    
        if consequence: 
            if context:
                if population:
                    if cancer: 
                        for region, regionFreq in population.items():
                            for cancerType, cancerFreq in cancer.items():
                                self.variantPlot[region +","+ regionFreq+","+cancerType+","+ cancerFreq  +","+ coding+","+mutation] += 1
                                self.contextPlot[region +","+ regionFreq+","+cancerType+","+ cancerFreq  +","+ coding+","+context] += 1
                                self.consequencePlot[region+","+regionFreq+","+cancerType+","+cancerFreq+","+coding+","+consequence] += 1
                                
                        self.putInVariantTable(record, recordVEP, context, consequence, population, cancer, coding)
                        return True
        return False

    def putInVariantTable(self, record, record_VEP, context, consequence, population, cancer, coding):
        gdnaPos = record.CHROM+":"+str(record.POS)+":"+record.REF+">"+record.ALT[0]
        variant = "GRCh37,"
        variant += gdnaPos+","
        variant += record_VEP['Consequence']+","
        variant += record_VEP['Gene']+","
        variant += record_VEP['SYMBOL']+","
        variant += record_VEP["BIOTYPE"]+","
        variant += ''.join(record.INFO["COSMIC_MUTATION_ID"].split(','))+","
        cancerList = []
        for cancerType in record.INFO['COSMIC_CANCER_TYPE_GW'].split('&'):
            c = cancerType.split(':')
            if int(c[1]) > 1:
                cancerFreq = "(recurrent)"
            else:
                cancerFreq = "(nonrecurrent)"
            if c[0] in cancer:   
                cancerList.append(cancerType+cancerFreq)

        variant += "&".join(cancerList)
        variant += ",variant consequece="+coding
        region = []
        for pop, freq in population.items():
            region.append(pop+":"+str(record.INFO[pop])+"("+freq+")")
        variant +=  ";population and allele frequency="+'&'.join(region)

    
        self.variantTable.append(variant)
    

        
    def getVariantPlot(self):
        return self.variantPlot
        
    def getContextPlot(self):
        return self.contextPlot
        
    def getConsequencePlot(self):
        return self.consequencePlot
    
    def getHeader(self):
        return self.csv_header

class COSMIC_oneKG(WriteCSV):

    def __init__(self):
        self.mt = an.MutationType()
        self.mc = an.MutationContext()
        self.cq = an.Consequence()
        self.ct = an.CancerType()
        self.pc = an.Population('1kg')
        self.csv_header = "plot,set,region,alleleFreq,cancer,cancerFreq,coding,type,n\n"
        self.variantPlot = OrderedDict()
        self.contextPlot = OrderedDict()
        self.consequencePlot = OrderedDict()
        self.variantTable = []
             
        for region, pop_freq_list in self.pc.getDict().items():
            for regionFreq in pop_freq_list.keys():
                for cancer, cancer_freq_list in self.ct.getDict().items():
                    for cancerFreq in cancer_freq_list.keys():
                        for coding in self.cq.getCodingDict():
                             
                             for mutation in self.mt.getDict().keys():
                                 self.variantPlot[region+","+regionFreq+","+cancer+","+cancerFreq+","+coding+","+mutation] = 0
                             for context in self.mc.getDict().keys():
                                 self.contextPlot[region+","+regionFreq+","+cancer+","+cancerFreq+","+coding+","+context] = 0
                        for consequence, coding  in self.cq.getDict().items():    
                            self.consequencePlot[region+","+regionFreq+","+cancer+","+cancerFreq+","+coding.keys()[0]+","+consequence] = 0                            

    def regVariant(self, record, recordVEP, referenceSequence):
        
        mutation = self.mt.getType(record)
        context  = self.mc.getMutContext(referenceSequence, record)
        consequence = self.cq.getConsequence(recordVEP)
        coding = self.cq.getCodingType(recordVEP)
        cancer = self.ct.getCancer(record)
        coding = self.cq.getCodingType(recordVEP)
        population = self.pc.getPopulation(record)

        if consequence: 
            if context:
                if population:
                    if cancer: 
                        for region, regionFreq in population.items():
                            for cancerType, cancerFreq in cancer.items():
                                self.variantPlot[region +","+ regionFreq+","+cancerType+","+ cancerFreq  +","+ coding+","+mutation] += 1
                                self.contextPlot[region +","+ regionFreq+","+cancerType+","+ cancerFreq  +","+ coding+","+context] += 1
                                self.consequencePlot[region+","+regionFreq+","+cancerType+","+cancerFreq+","+coding+","+consequence] += 1
                        self.putInVariantTable(record, recordVEP, context, consequence, population, cancer, coding)
                        return True
        return False

    def putInVariantTable(self, record, record_VEP, context, consequence, population, cancer, coding):
        gdnaPos = record.CHROM+":"+str(record.POS)+":"+record.REF+">"+record.ALT[0]
        variant = "GRCh37,"
        variant += gdnaPos+","
        variant += record_VEP['Consequence']+","
        variant += record_VEP['Gene']+","
        variant += record_VEP['SYMBOL']+","
        variant += record_VEP["BIOTYPE"]+","
        variant += ''.join(record.INFO["COSMIC_MUTATION_ID"].split(','))+","
        cancerList = []
        for cancerType in record.INFO['COSMIC_CANCER_TYPE_GW'].split('&'):
            c = cancerType.split(':')
            if int(c[1]) > 1:
                cancerFreq = "(recurrent)"
            else:
                cancerFreq = "(nonrecurrent)"
            if c[0] in cancer:   
                cancerList.append(cancerType+cancerFreq)

        variant += "&".join(cancerList)
        variant += ",variant consequece="+coding
        region = []
        for pop, freq in population.items():
            region.append(pop+":"+str(record.INFO[pop])+"("+freq+")")
        variant +=  ";population and allele frequency="+'&'.join(region)

        self.variantTable.append(variant)

    
    def getVariantPlot(self):
        return self.variantPlot
        
    def getContextPlot(self):
        return self.contextPlot
        
    def getConsequencePlot(self):
        return self.consequencePlot
    
    def getHeader(self):
        return self.csv_header


        

