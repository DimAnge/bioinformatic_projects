import argparse

file = input('Type file path: ')

#Reads Fasta File
def read_file(file_path):
    with open(file,'r') as f:
        l = [l.strip() for l in f.readlines()]
    if l[0].startswith('>') :
        seq = l[1]
        header = l[0][1:]
    else:
        seq = l[0]
    return seq

#Validates DNA Sequence
def validate_dna(seq):
    dna_nucl = ['A','T','G','C','a','t','g','c']
    for i in seq:
        if i not in dna_nucl:
            return 'Invalid Sequence. Please input a DNA sequence'
        else:
            return seq.upper()

#Calculates GC content
def gc_content(seq):
    return round(100*(seq.count('G') + seq.count('C')) / len(seq)) 

#Reverse Complement Strand
def reverse_compl(seq):
    di_repl = {'A':'T','T':'A','C':'G', 'G':'C' }
    compl = seq.translate(str.maketrans(di_repl))
    return compl[::-1]

#Nucleotide Frequency Count
def nucl_freq(seq):
    di_nu = {}
    for nucl in seq:
        if nucl in di_nu:
            di_nu[nucl] += 1
        else : 
            di_nu[nucl] = 1
    return ', '.join([f'{nucl} : {count / len(seq):.2f}' for nucl, count in di_nu.items()])

#Transcribe DNA to RNA
def transcribe(seq):
    return seq.replace('T','U')
          
seq = validate_dna(read_file(file))
print(nucl_freq(seq))   