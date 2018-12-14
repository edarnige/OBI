'''
Projet Alone in the Dark
myProject
2018.11.28
Juan Garcia - Eden Darnige - Alexandre Lambard
Python 2
'''

import MyBio as mybio

def getGeneticCode(trans_table) : 
    """Returns a genetic code table
    
    This function returns a genetic code table, where DNA codons and aminoacids are related. This 
    function is based on the information extracted from the NCBI website https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. 
    Based on the standard code table, modifications are added following the guidelines given on the webpage.
    In order to avoid time-consuming processes, only a limited number of genetic tables will be 
    loaded
    This functions has been written by Eden Darnige, Juan Manuel Garcia and Alexandre Lambard
    
    Args: 
        trans_table: NCBI genetic table ID. This ID is a non-zero number
        
    Returns: 
        A dictionnary where codons are keys and aminoacids are values. This is the format: 
        {"Codon_1": "AA_1" , "Codon_2": "AA_2" , ...}
     
    Raises: 
        NumericValue: NCBI ID is not correct  
    
    """
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
        raise NumericValue("NCBI ID incorrect")


def findORF(seq,threshold,codeTable) : 
    """Returns a list of Open Reading Frames (ORF)
    
    This function retrieves all the ORF present in a DNA sequence, extracting several data about each 
    one of them. Since it works with protein sequences to find ORFs, a conversion on the indexes and on the
    threshold is neccesary to be made. A threshold is used to filter out small ORFs.
    This function has been written by Eden Darnige, Juan Manuel Garcia and Alexandre Lambard
    
    Args: 
        seq: A dictionnary whose key is a description of the DNA sequence (i.e., header of a FASTA file) 
        and whose value is a string of a DNA sequence
        threshold: A number, expressed in nucleotide base-pairs. 
        codeTable: A genetic code table, like the one that can be obtained using getGeneticCode. 
        
    Returns: 
        A list of dictionnaries, where information about an ORF is stored in each one of them. The format 
        of the dictionnaries is as follows: 
            ...
            listofORF[n]["start"]: starting index in the DNA sequence of nth ORF 
            listofORF[n]["stop"]: ending index in the DNA sequence of the nth ORF
            listofORF[n]["frame"]: reading frame 
            listofORF[n]["protein"]: protein sequence that corresponds to the ORF
            listofORF[n]["length"]: length, in base-pairs, of the ORF
            ...
            
    Raises: 
        TypeError: seq is not a dictionnary
        
    """
    ORF=[] 
    if type(seq) not dict : 
        raise TypeError
        
    threshold=threshold//3 # Since protein sequences are taken, threshold must be modified
    count=0
    seq_frames=mybio.translate_frame(seq,12,codeTable)
    
    seq=seq.values()[0]
    for i in range(0,len(seq_frames.values())) : 
        if (mybio.countWord(seq_frames.values()[i],codeTable["ATG"])>0 and mybio.countWord(seq_frames.values()[i],"_")>0) : 
            index=0

            while index<len(seq_frames.values()[i]) : 
                AA_end=seq_frames.values()[i].find("_",index)  # Start protein index
                AA_start=seq_frames.values()[i].find(codeTable["ATG"],index) # Stop protein index 
                if (AA_end==-1 or AA_start==-1) : 
                    break
                if AA_end-AA_start>= threshold : 
                    frame=int(seq_frames.keys()[i].replace("frame",""))
                    '''
                    In the following step, the conversion of start and stop protein index into a valid nucleic 
                    sequence index is done. Since the frame value given might be negative, two different situations
                    are considered.
                    '''
                    if frame<0 : 
                        start_index=len(seq)-(abs(frame)-1)-AA_end*3-3 
                        end_index=len(seq)-(abs(frame)-1)-AA_start*3-3
                    else : 
                        start_index=(frame-1)+AA_start*3
                        end_index=(frame-1)+AA_end*3
                    
                    ORF.append({"start": start_index , "stop": end_index , "frame": frame , "protein":seq_frames.values()[i][AA_start:AA_end] , "length": ((AA_end-AA_start)*3) , "length1" : end_index-start_index})
                index=AA_end+1

    return ORF

             
def getLengths(orfList) : 
    """ Returns the length of ORFs
    
    This function retrieves the length of all the ORF of a DNA sequence. This will be stored in a list with 
    dictionnaries, where the index and the length is added for each ORF. 
    This functions has been written by Eden Darnige, Juan Manuel Garcia and Alexandre Lambard
    
    Args: 
        orfList: A list of ORF, obtained from the findORF function.
        
    Returns: 
        A list with dictionnaries with index and length value for each ORF. 
        
    Raises: 
                
    """
    length=[]
    for i in range(0,len(orfList)) : 
        value=orfList[i]["length"]
        length.append({"ORF": i , "length": value})

    return length

def getLongestORF(orfList) : 
    """Returns the longest ORF 
    
    This function returns an ORF whose length is the highest one of a list of ORF. 
    It has been written by Eden Darnige, Juan Manuel Garcia and Alexandre Lambard

    Args: 
        orfList: A list of ORF, obtained from the findORF function.
        
    Returns: 
        A list with information about the ORF. This will be extracted from the ORF list obtained from 
        findORF. Therefore, it will have the same format as the ORF list (same keys within the dictionnary)
        
    Raises: 
            
    """
    length=getLengths(orfList)
    max_val=0
    for i in range(0,len(length)) : 
        if length[i]["length"]>max_val : 
            index=i
            max_val=length[i]["length"]
    return orfList[index]

def TopLongestORF(orfList,value) : 
    """Extracts the n% longest ORFs 
    
    This function retrieves the longest ORFs from a ORF list. To do so, a value which will be used as a threshold 
    should be given. An index is calculated using the value given by the user.  The ORF list is sorted and the 
    longest ORFs are extracted. 
    It has been written by Eden Darnige, Juan Manuel Garcia and Alexandre Lambard
    
    Args: 
        orfList: A list of ORF, obtained from the findORF function.
        value: A number between 0 and 1, which represents the nth longest ORFs
        
    Returns: 
        A list of ORFs with information about each of them. This list is extracted from the main ORF list. 
        
    Raises: 
           
    """
    import math
    from operator import itemgetter
    length=len(orfList)

    nombre=math.floor(length * value)
    nombre=int(nombre) # Avoid using floating number
    
    orfList=sorted(orfList,key=itemgetter('length')) # Sort list according to values in ["length"] key 
    
    return orfList[-nombre:]


##### READ AND WRITE CSV FILES 

def readCSV(filename,separator) : 
    """Prints a cvs file 
    This function prints a csv file by storing each row in a dictionary. All the dictionnaries are grouped in
    a list. 
    This function has been written by Eden Darnige 
    
    Args: 
        filename: filename of the csv file, as a string
        separator: row separator, as a string
    
    Returns: 
        A list with dictionnaries, where the keys are the headers and the values are the data stored 
        in the csv file 
        
    Raises: 
        
    
    """
    import csv
    ORF=[]
    with open(filename,'r') as csv_file: 
        reader=csv.DictReader(csv_file,delimiter=separator)
        line_count=0
        for row in reader :
            dicc={} 
            if line_count==0: 
                keys=row
                dicc=rewrite(dicc,keys,"key")
                line_count+=1
            dicc=rewrite(dicc,row,"value") 
            ORF.append(dicc)
            line_count+=1
    return ORF

def writeCSV(filename,separator,dictionary) : 
    '''Writes a ORF dictionnary into a CSV file
    This function loads a Python dictionnary into a CSV file. Keys of the dictionnary will be put in the 
    first line of this file as header, whereas the values will be in their corresponding columns. 
    This function has been written by Eden Darnige. 
    
    Args: 
        filename: filename of the csv file, as a string
        separator: row separator, as a string
        dictionnary: ORF dictionnary where the data is stored
    
    Returns: 
        A CSV file with the content of the dictionnary 
        
    Raises: 
        
    '''
    import csv
    with open(filename,'w') as file_name: 
        file_temp_writer=csv.writer(file_name,delimiter=separator)
        keys=dictionary[0].keys()
        file_temp_writer.writerow(keys)
        for d in dictionary:
            temp_list=[]
            for value in d.values():
                temp_list.append(value)
            file_temp_writer.writerow(temp_list)


         









