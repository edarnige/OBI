
'''
Projet Alone in the Dark
myBio
2018.11.28
Juan Garcia - Eden Darnige - Alexandre Lambard
Python 2
'''

def mystere(prot):
    nb_cys = 0
    for aa in prot:
        if aa == "C" :
            nb_cys = nb_cys + 1
    return nb_cys



def isDNA(seq) : 
    for i in seq : 
        if (i!="c" and i!="t" and i!="a" and i!="g" ) : 
            return False 
    return True



def countPro(seq) : 
    cont=0
    for i in seq : 
        if (i.isupper() and i=="P") : 
            cont=cont+1
            
    print "There are %d Proline aa in the protein sequence" %cont
    return cont


def countAll(seq,aa) : 
    cont=0
    for i in seq : 
        if (i==aa and i.isupper()) : 
            cont=cont+1
    print "There are %d %s aa in the protein sequence" %(cont,aa)
    return cont

def oneWord(seq,start,wlen) : 
    if start>=0 and start+wlen<len(seq): 
        
        word=""
        it=0
        while it<wlen : 
            
            word=word+seq[start+it]
            it=it+1
        return word    
    else : 
        print "Starting index is not correct"
    


def countWord(seq,word) : 
    if word in seq : 
        cont=0
        for i in range(0,len(seq)-len(word)) : 
            w=oneWord(seq,i,len(word))
            if w==word : 
                cont=cont+1
        return cont
    else : 
        print "The word is not part of the sequence"


def isCodonStart(seq,pos) : 
    codon=oneWord(seq,pos,3)
    if codon=="ATG" : 
        return True
    else : 
        return False

def isCodonStop(seq,pos) : 
    codon=oneWord(seq,pos,3)

    if (codon=="TAA" or codon=="TAG" or codon=="TGA") : 
        return True
    else : 
        return False


def isGene(seq) : 
    rem=len(seq)%3 
    cod_start=0
    cod_fin=0
    for i in range(0,len(seq)-rem,3) : 
        if isCodonStart(seq,i) : 
            cod_start=1
        elif isCodonStop(seq,i) : 
            cod_fin=1
    if cod_start==1 and cod_fin==1 : 
        return True
    else : 
        return False

def isGene3(seq) : 
    rem=len(seq)%3
    for i in range(0,len(seq)-rem) :
        if isCodonStart(seq,i) :
            frame_start=i%3 
            if isCodonStop(seq,i) : 
                frame_stop=i%3
                if frame_start==frame_stop : 
                    return True
                else : 
                    return False 

## 3

def readFASTA(filename) :   
    fasta = {}
    file_one=open(filename,'r')
    for line in file_one:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            active_sequence_name = line[1:]
            if active_sequence_name not in fasta:
                fasta[active_sequence_name] = []
            continue
        sequence = line.upper()
        fasta[active_sequence_name].append(sequence)
    
    fich=""
    for i in range(len(fasta[active_sequence_name])) : 
        fich=fich+fasta[active_sequence_name][i]
    
    fasta[active_sequence_name]=fich

    return fasta

def writeFASTA(seq,filename) : 
    if (type(seq)==str and type(filename)==str and isDNA(seq)) : 
        fichier=open(filename,'w')
        header=raw_input("Enter header of FASTA file: ")
        fichier.write(">"+header)
        fichier.write("\n")
        fichier.write(seq.upper())
        fichier.close()
    else : 
        print "Not working..."


def translate(nucl_seq,codonTable):
    protein_seq=''
    for n in range(0, len(nucl_seq), 3):
        if nucl_seq[n:n+3] in codonTable:
            protein_seq += codonTable[nucl_seq[n:n+3]]
    return protein_seq

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


def complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    bases=seq
    bases = [complement[base] for base in bases] 
    return ''.join(bases)

def translate_frame(seq,frame,codonTable) : 
    seq=seq.values()[0]
    sequences={}
    rev_seq=complement(seq)
    rev_seq=rev_seq[::-1]
    for j in range(3) :
        protein_seq='' 
        rev_protein_seq=''
        for i in range(j,len(seq),3) : 
            if seq[i:i+3] in codonTable:
                protein_seq += codonTable[seq[i:i+3]]
        sequences["frame{}".format(j+1)]=protein_seq

        for i in range(j,len(rev_seq),3) : 
            if seq[i:i+3] in codonTable:
                rev_protein_seq += codonTable[rev_seq[i:i+3]]
        sequences["frame{}".format(-(j+1))]=rev_protein_seq


        ### frame-1 est frame-3, et vice versa
    if (-3<=frame<=3 and frame!=0) : 
        string="frame"+str(frame)
        return {'string': sequences[string]}
    elif frame==-6 : 
        return {"frame-1": sequences["frame-1"],"frame-2": sequences["frame-2"],"frame-3": sequences["frame-3"]}
    elif frame==6 : 
        return  {"frame1": sequences["frame1"],"frame2": sequences["frame2"],"frame3": sequences["frame3"]}
    elif frame==12 : 
        return sequences
    else : 
        print "Frame error"






    

