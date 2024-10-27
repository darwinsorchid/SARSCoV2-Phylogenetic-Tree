# Phylogenetic Tree Construction of SARS-CoV-2 Genomes

## Description
Construction of phylogenetic tree from NCBI retrieved SARS-CoV-2 genomes in GenBank format.
The tree was constructed from five viral sequences per continent or thirty sequences in total. 
The selected sequences accession numbers (file: `accessions.txt`) were chosen from:
[NCBI Virus Sequences Dashboard for SARS-CoV-2 ](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Nucleotide&VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202,%20taxid:2697049).


## MSA & Tree Contruction Algorithms
For the construction of the phylogenetic tree, multiple sequence alignment of the SARS-CoV-2 genomic sequences was performed with __Clustal Omega__ and the tree was created with the
__UPGMA__ (Unweighted Pair Group Method with Arithmetic Mean) algorithm, using a distance matrix calculated using the MSA.


