import argparse
import requests

def genbank(ids):
    url = "https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi"
    parameters = {'db': 'nuccore', 'report' : 'fasta', 'id' : '', 'withparts': 'on', 'retmode' : 'text'}
    
    for id in ids:
        parameters['id'] = id
        response = requests.get(url, params = parameters)
        if response.status_code == 200:
            filename = f"{id}.fasta"
            with open(filename, 'w') as file:
                file.write(response.text)
            print(f"Sequence {id} downloaded successfully.")
        else:
            print(f"Failed to download sequence {id}. HTTP Status Code: {response.status_code}")
            
            
def uniprot(ids):
    
    for prot_id in ids:
        url = f'https://rest.uniprot.org/uniprotkb/{prot_id}.fasta'
        response = requests.get (url)
        if response.status_code == 200:
            filename = f"{prot_id}.fasta"
            with open(filename, 'w') as file:
                file.write(response.text)
            print(f"Sequence {prot_id} downloaded successfully.")
        else:
            print(f"Failed to download sequence {prot_id}. HTTP Status Code: {response.status_code}")
            
            
def ensembl(ids):
    url = "https://rest.ensembl.org/sequence/id/"
    headers = {'Content-Type': 'text/x-fasta'}

    for seq_id in ids:
        response = requests.get(f"{url}{seq_id}", headers=headers, params={'type': 'genomic'})

        if response.status_code == 200:
            filename = f"{seq_id}.fasta"
            with open(filename, 'w') as file:
                file.write(response.text)
            print(f"Sequence {seq_id} downloaded successfully.")
        else:
            print(f"Failed to download sequence {seq_id}. HTTP Status Code: {response.status_code}")
            
def pdb(ids):
    url = "https://files.rcsb.org/download/"
    
    for id in ids:
        response = requests.get(f"{url}{id}.pdb")

        if response.status_code == 200:
            filename = f"{id}.pdb"
            with open(filename, 'w') as file:
                file.write(response.text)
            print(f"Pdb file of: {id} downloaded successfully.")
        else:
            print(f"Failed to Pdb file of {id}. HTTP Status Code: {response.status_code}")
    
    
            
def read_file(file):
    ids = []
    with open(file, 'r') as f:
        for line in f:
            ids.append(line.strip())
    return ids
            
            
def main():
    parser = argparse.ArgumentParser(description='This tool reads a file of IDs (separated with \\n) and downloads the corresponding sequence from GenBank(--genbank), Ensembl(--ensembl), Uniprot(--uniprot) in fasta format or the corresponding pdb file(--pdb).')
    parser.add_argument('input', help='Input Genbank ID (--genbank)/Uniprot ID (--uniprot) /Ensemble ID/ PDB ID.[ex] seq_fetcher <file_path> --genbank --uniprot --ensembl --pdb')
    parser.add_argument('--genbank', action='store_true', help='Downloads the sequences in fast format of the GenBank id')
    parser.add_argument('--uniprot', action='store_true', help='Downloads the sequence in fasta format of the Uniprot Accession Number')
    parser.add_argument('--ensembl', action='store_true', help='Downloads the sequences in fast format of the Esnembl id' )
    parser.add_argument('--pdb', action='store_true', help='Downloads PDB files.' )
    
    
    args = parser.parse_args()
    ids = read_file(args.input)
    
    
    if args.genbank:
        genbank(ids)
    
    if args.uniprot:
        uniprot(ids)
        
    if args.ensembl:
        ensembl(ids)
        
    if args.pdb:
        pdb(ids)
        
if __name__ == '__main__':
    main()

        