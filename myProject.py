#from td1 import *

import td1 as mybio

def getGeneticCode(trans_table) : 
    tableau={
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
    if trans_table==1: 
        return tableau
    elif trans_table==2 : 
        tableau["AGA"]="_"
        tableau["AGG"]="_"
        tableau["ATA"]="M"
        tableau["TGA"]="W"
        return tableau
    elif trans_table==3 : 
        tableau["ATA"]="M"
        tableau["CTT"]="T"
        tableau["CTC"]="T"
        tableau["CTA"]="T"
        tableau["CTG"]="T"
        tableau["TGA"]="W"
        return tableau
    elif trans_table==4: 
        tableau["TGA"]="W"
        return tableau
    elif trans_table==5 : 
        tableau["AGA"]="S"
        tableau["AGG"]="S"
        tableau["ATA"]="M"
        tableau["TGA"]="W"
        return tableau
    elif trans_table==6 : 
        tableau["TAA"]="Q"
        tableau["TAG"]="Q"
        return tableau
    else:
        print "NCBI ID incorrect"

seq=mybio.readFASTA("/media/juagarcia/JUANMA/Python/BioPython/tes.fas")


def findORF(seq,threshold,codeTable) : 
    #seq=seq.values()[0]
    ORF=[] 

    threshold=threshold//3
    count=0
    seq_frames=mybio.translate_frame(seq,12,codeTable)
    
    seq=seq.values()[0]
    for i in range(0,len(seq_frames.values())) : 
        if (mybio.countWord(seq_frames.values()[i],codeTable["ATG"])>0 and mybio.countWord(seq_frames.values()[i],"_")>0) : 
            index=0

            while index<len(seq_frames.values()[i]) : 
                AA_end=seq_frames.values()[i].find("_",index)
                AA_start=seq_frames.values()[i].find(codeTable["ATG"],index)
                print("e ",AA_end)
                print("s ",AA_start)
                if (AA_end==-1 or AA_start==-1) : 
                    break
                if AA_end-AA_start>= threshold : 
                    frame=int(seq_frames.keys()[i].replace("frame",""))
                    if frame<0 : 
                        start_index=len(seq)-(abs(frame)-1)-AA_end*3-3
                        end_index=len(seq)-(abs(frame)-1)-AA_start*3-3
                    else : 
                        start_index=(frame-1)+AA_start*3
                        end_index=(frame-1)+AA_end*3
                    
                    ORF.append({"start": start_index , "stop": end_index , "frame": frame , "protein":seq_frames.values()[i][AA_start:AA_end] , "length": ((AA_end-AA_start)*3) , "length1" : end_index-start_index})
                index=AA_end+1

    return ORF
            #  sub_seq=seq_frames.values()[i].split("_")

            #  for j in range(0,len(sub_seq)) :
            #      if "M" in sub_seq[j] : 
            #          index=sub_seq[j].index("M")
            #          for k in index : 
            #              if len(sub_seq[j])-k>= threshold : 
            #                  ORF[count]={}
            #                  ORF["start"]=k*(j+1)
            #                  ORF["stop"]=len
            #                  ORF["frame"]=seq_frames.keys()[i].replace("frame","")
            #                  ORF["protein"]=sub_seq[j][index[k]:]
            #                  ORF["length"]=len(sub_seq[j][index[k]:])*3
            #                  count+=1

                                 

             
def getLengths(orfList) : 
    length=[]
    for i in range(0,len(orfList)) : 
        value=orfList[i]["length"]
        length.append({"ORF": i , "length": value})

    return length

def getLongestORF(orfList) : 
    length=getLengths(orfList)
    max_val=0
    for i in range(0,len(length)) : 
        if length[i]["length"]>max_val : 
            index=i
            max_val=length[i]["length"]
    return orfList[index]
listofORF=findORF(seq,3,getGeneticCode(1))
print listofORF
for i in range(0,len(listofORF)) : 
    check=1
    if listofORF[i]["length"]!=listofORF[i]["length1"] : 
        check=0

print check
def TopLongestORF(orfList,value) : 
    import math

    length=len(orfList)

    nombre=math.floor(length * value)

    return orfList[-nombre:]


##### READ AND WRITE CSV FILES 


def readCSV(filename) : 
    csv_file=open(filename,'r')
    reader=csv.reader(csv_file)
    mydict=dict((rows[0],rows[1]) for rows in reader)
    #mydict={k:int(v) for k, v in mydict.iteritems()}  #string to int
    print mydict

dicto={"a": 1 , "b": 2 , "c": 3, "d": 4}

def writeCSV(filename) : 
    csv_file=open(filename,'w')
    for key in dicto.keys() : 
        csv_file.write("%s,%s\n"%(key.dicto[key]))








