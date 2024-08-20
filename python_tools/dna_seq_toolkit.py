import argparse

#Reads Fasta File
def read_file(file):
    seq =''
    with open(file,'r') as f:
        for i in f.readlines():
            if i.startswith('>'):
                header = i.strip
            else:
                seq += i.strip()
    return validate_dna(seq)

#Validates DNA Sequence
def validate_dna(seq):
    dna_nucl = ['A','T','G','C','a','t','g','c']
    for i in seq:
        if i not in dna_nucl:
            return 'Invalid Sequence. Please input a DNA sequence.\n'
        else:
            return seq.upper()

#Calculates GC content
def gc_content(seq):
    return round(100*(seq.count('G') + seq.count('C')) / len(seq),2) 

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
    l = f'Length = {len(seq)}'
    return ', '.join([f'{nucl} : {count / len(seq):.2f}' for nucl, count in di_nu.items()]), l

#Transcribe DNA to RNA
def transcribe(seq):
    return seq.replace('T','U')

#Translate DNA sequence to Protein sequence until it finds a stop codon (Add more genetic codes)
def translate(seq):
    codon = {'UUU': 'F',    'CUU': 'L',      'AUU': 'I'   ,  'GUU' :'V',
'UUC' :'F'    ,   'CUC' :'L'   ,   'AUC' :'I'    ,  'GUC' :'V',
'UUA': 'L'    ,   'CUA' :'L'   ,   'AUA' :'I'    ,  'GUA' :'V',
'UUG' :'L'    ,   'CUG' :'L'   ,   'AUG' :'M'    ,  'GUG' :'V',
'UCU' :'S'    ,   'CCU' :'P'   ,   'ACU' :'T'    ,  'GCU' :'A',
'UCC' :'S'    ,   'CCC' :'P'   ,   'ACC' :'T'    ,  'GCC' :'A',
'UCA' :'S'    ,   'CCA' :'P'   ,   'ACA' :'T'    ,  'GCA' :'A',
'UCG' :'S'    ,   'CCG' :'P'   ,   'ACG' :'T'    ,  'GCG' :'A',
'UAU' :'Y'    ,   'CAU' :'H'   ,   'AAU' :'N'    ,  'GAU' :'D',
'UAC' :'Y'    ,   'CAC' :'H'   ,   'AAC' :'N'    ,  'GAC' :'D',
'UAA' :'Stop' ,   'CAA' :'Q'   ,   'AAA' :'K'    ,  'GAA' :'E',
'UAG' :'Stop' ,   'CAG' :'Q'   ,   'AAG' :'K'    ,  'GAG' :'E',
'UGU' :'C'    ,   'CGU' :'R'   ,   'AGU' :'S'    ,  'GGU' :'G',
'UGC' :'C'    ,   'CGC' :'R'   ,   'AGC' :'S'    ,  'GGC' :'G',
'UGA' :'Stop' ,   'CGA' :'R'   ,   'AGA' :'R'    ,  'GGA' :'G',
'UGG' :'W'    ,   'CGG' :'R'   ,   'AGG' :'R'    ,  'GGG' :'G' }
    protein =''
    rna = transcribe(seq)
    for i in range(0, len(rna)-3, 3):
        if i == 'Stop':
            break
        else:
            protein += codon[rna[i:i+3]]
    return protein
    
#Checks for palindromes in the sequence (ex. restriction enzymes) and prints the palindrome and the position / window method        
def palindrome(seq):
    di_repl = {'A':'T','T':'A','C':'G', 'G':'C' }
    found = []
    for i in range(len(seq)):
        for j in range(i,len(seq)):
            array = seq[i:j+1]
        
            if  3 < len(array) < 13:
                compl_array = array.translate(str.maketrans(di_repl))
                compl_array = "".join(reversed(compl_array))
                
                if array == compl_array:
                    pos = i+1
                    
                    if (pos, len(compl_array)) in found:
                        continue
                    
                    else:
                        found.append((pos,compl_array))
    return found

#Finds all open reading frames
def orf_finder(seq):
    start_codon = 'ATG'
    end_codon = ['TAA','TAG','TGA']
    orfs = []
    dna = [(seq,'+'),(reverse_compl(seq),'-')]
    for seq,strand in dna:
        l = len(seq)
        for i in range(l - 2):
                temp = seq[i:i+3]
                
                if temp == start_codon:
                    orf = ''
                    found_stop = False
            
                    for j in range(i, l - 2, 3):
                        codon = seq[j:j+3]
                        
                        if codon in end_codon:
                            found_stop = True
                            break
                        orf += codon
                    
                    if found_stop:
                        start_index = i + 1
                        end_index = i + len(orf)
                        if strand == '-':
                            start_index = l - i - len(orf) +1
                            end_index = l - i 
                    
                        orfs.append((end_index, start_index, strand, orf))
        
    return orfs


#Main function/arg parser
def main():
    parser = argparse.ArgumentParser(description='DNA Sequence Toolkit. This tool takes as input a DNA sequence in FASTA or txt format and performs various functions.')
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('--output', help='Output file to save results')
    parser.add_argument('--gc', action='store_true', help='Calculate GC content')
    parser.add_argument('--rev_compl', action='store_true', help='Get reverse complement')
    parser.add_argument('--transcribe', action='store_true', help='Transcribe DNA to mRNA')
    parser.add_argument('--translate', action='store_true', help='Translate DNA to amino acid sequence')
    parser.add_argument('--count', action='store_true', help='Counts the frequency of each nucleotide')
    parser.add_argument('--palindrome', action ='store_true', help='Checks for palindromes in the sequence (ex. restriction enzymes)')
    parser.add_argument('--orf_finder', action= 'store_true', help='Finds all possible Open Reading Frames (nested orfs) for both strands.')
    
    args = parser.parse_args()
    seq = read_file(args.input)
    res = []
    
    if args.count:
        res.append(f'Nucleotide Frequency Count: {nucl_freq(seq)}')
    
    if args.gc:
        res.append(f'GC Content : {gc_content(seq)}')
        
    if args.rev_compl:
        res.append(f'Reverse Complement Strand: {reverse_compl(seq)}')
    
    if args.transcribe:
        res.append(f'mRNA: {transcribe(seq)}')
        
    if args.translate:
        res.append(f'Amino Acid Sequence: {translate(seq)}')
        
    if args.palindrome:
        palindromes = "\n".join([f"{pal[0]}  {pal[1]}" for pal in palindrome(seq)])
        res.append(f'Palindromes Found:\n {palindromes}')
        
    if args.orf_finder:
        for item in orf_finder(seq):
            head = "".join(f'>Found ORF of Length: {len(item[3])} |Strand: {item[2]} |Start: {item[0]} End: {item[1]} in original strand')
            dna = item[3]
            res.append(f'{head}\n{dna}')        
        
    if args.output:
        with open(args.output, 'w') as out:
            out.write('\n'.join(res))
    
    else:
        for res in res:
            print(res)    


if __name__ == '__main__':
    main()


