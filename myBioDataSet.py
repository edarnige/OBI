'''
Projet Alone in the Dark
myProject
2018.11.28
Juan Garcia - Eden Darnige - Alexandre Lambard
Python 2
'''


import myBio as mybio
import myProject as proj
import csv


seq = mybio.readFASTA('/home/parallels/Desktop/Omics/OBI/Haemophilus_influenzae_s723.fasta')
codeTable = proj.getGeneticCode(1)
listORFs= proj.findORF(seq,3*420,codeTable)

print "Threshold 420:"
# Number of ORFs for this threshold
print "how many ORFs:",len(listORFs)  

# #Example of ORF data
print "Example ORF 0"
print(listORFs[0]['start'])
print(listORFs[0]['stop'])
print(listORFs[0]['frame']) # 1, 2, 3, -1, -2, or -3
print(listORFs[0]['length']) #  in bp
print(listORFs[0]['protein']) # M.....*


#List of lengths of all ORFs in this threshold
lengthlist = proj.getLengths(listORFs)
#print "lengthlist",lengthlist

#Get info about the longest ORF
orf = proj.getLongestORF(listORFs)
#print "The longest ORF:",orf

#returns dictionary of ORFs in the top x percentile 
orfs = proj.TopLongestORF(listORFs,0.25)




dicto=listORFs[:]
# proj.writeCSV('ORF_420.csv',dicto)


#Threshold 300
print "Threshold 300:"
listORFs= proj.findORF(seq,3*300,codeTable)
print "how many ORFs:",len(listORFs)
dicto=listORFs[:]
#proj.writeCSV('ORF_300.csv',dicto)


#Threshold 210
print "Threshold 210:"
listORFs= proj.findORF(seq,3*210,codeTable)
print "how many ORFs:",len(listORFs)
dicto=listORFs[:]
#proj.writeCSV('ORF_210.csv',dicto)


#Threshold 90
print "Threshold 90:"
listORFs= proj.findORF(seq,3*90,codeTable)
print "how many ORFs:",len(listORFs)
dicto=listORFs[:]
#proj.writeCSV('ORF_90.csv',dicto)



#Threshold 0
print "Threshold 0"
listORFs= proj.findORF(seq,3*0,codeTable)
print "how many ORFs:",len(listORFs)
dicto=listORFs[:]
#proj.writeCSV('ORF_0.csv',dicto)
