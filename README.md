# Bioinformatics Projects - Script Repository
## Overview

Welcome to my Bioinformatics and Computational Biology repository. 

Repository for writing/practising code and developing tools that will help me at my pursuit in a career in bioinformatics. The main goal is to improve my bioinformatic skills by performing analyses, such as differential expression analysis, pathway enrichment, and developing various tools in the field of computational biology.

### Mini Projects

#### Differentional Expression Analysis
* [edgeR](https://github.com/DimAnge/bioinformatic_projects/tree/main/differential_expression_analysis/edgeR)
* [DESeq2](https://github.com/DimAnge/bioinformatic_projects/tree/main/differential_expression_analysis/DESeq2)

#### Gene Set Enrichment Analysis
* [clusterProfiler](https://github.com/DimAnge/bioinformatic_projects/tree/main/gene_enrichment_analysis/cluster_profiler)

#### Python Tools
* [dna_seq_toolkit](https://github.com/DimAnge/bioinformatic_projects/blob/main/python_tools/dna_seq_toolkit.py)

 This tool takes as input a DNA sequence in FASTA or txt format and performs various functions.  
 <details> <summary>Example</summary> dna_seq_toolkit.py file_path -h --orf_finder --transcribe --translate --reverse_comp --nucl_freq --palindrome --gc </details>  
 
 * [sequence_fetcher](https://github.com/DimAnge/bioinformatic_projects/blob/main/python_tools/sequence_fetcher.py)

This tool reads a file of IDs (separated with \\n) and downloads the corresponding sequence from GenBank, Ensembl, Uniprot in fasta format or the corresponding pdb file.  
<details> <summary>Example</summary> sequence_fetcher.py  file_path  --genbank  --uniprot  --ensembl  --pdb </details>

