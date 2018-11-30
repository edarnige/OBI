'''
Projet Alone in the Dark
myProject
2018.11.28
Juan Garcia - Eden Darnige - Alexandre Lambard
Python 2
'''

import myBio as bio

codon2aa = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

#### Examples from chapter 3 ####
# if bio.isDNA(seq):
#     result = bio.translate_frame(seq,2,codon2aa)
# else:
#     print('Sorry, this is not a nucleic sequence')

# bio.writeFASTA(result,'./Haemophilus_frames.fasta')

#####


#### Chapter 4 ####


def getGeneticCode(NCBI_ID):


myTable = getGeneticCode(1)


#Dictionary of DNA code is saved in variable 'code'
code = bio.readFASTA('./Haemophilus_influenzae_s723.fasta')
#Sequence only is saved in variable 'seq'
seq = code.values()[0]
print 'The length of the sequence is' len(seq)


def findORF(seq, threshold, codeTable):
    #problem 4.6.1.2 gives algorithm
    isGene3(seq)
    translate_frame(seq,complement,frame,codonTable)
   


listORFs =findORF(DNAseq,3*90,myTable)
print len(listORFs) #Number of ORFs

# for each ORF:
print(len(listORFs) ) # Number of ORFs
print(listORFs[0]['start'])
print(listORFs[0]['stop'])
print(listORFs[0]['frame']) # 1, 2, 3, -1, -2, or -3
print(listORFs[0]['length']) # as a tuple expressed in bp
print(listORFs[0]['protein']) # M.....*
print(listORFs[1]['start'])


def getLengths(orf_list):
    for i in orf_list:
        print "The length of ORF",i,"is", len(i)
    
lengthlist = getLengths(listORFs)

def getLongestORF(orflist):
    max = orflist[0]
    for i in orflist:
        if i > max:
            max = i
    return max

orf = getLongestORF(listORFs)

def getTopLongestORF(orflist,value):
    size=len(orflist)
    limit = value * size # = amount of elements in this percentile
    for i in range((size-limit),size) #??
        return i


topOrfs = getTopLongestORF(listORFs,0.25) # 25% percentile




dicto={}


def readCSV(filename, separator,data):


def writeCSV(filename, separator):

error  = proj.writeCSV('/my_home/my_project/my_data.csv',';',dicto)
result = readCSV('/my_home/my_project/my_data.csv',';')

if (dicto == result):
    print('OK')
else:
    print('KO')

def compare(orflist1,orflist2):


##### Genebank parser #####









