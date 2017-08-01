# guideRNA\_finder

### Description
The 2017 iGEM team has developed a modification of Howard Salis and Iman Farasats' code for classifying the binding affinity of guide RNAs against a target sequence (e.g. genome or extrachromosomal element).
This is described thoroughly in their [2016 PLoS Computational Biology paper](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004724).
 This code is used for finding the best guide RNA sequences in *Agrobacterium rhizogenes*.

### Usage
python Cas9\_Calculator/Cas9\_Calculator.py -g genome.gbk \[genome\_n.gbk] -m MODEL.mat -t target\_seq.fasta

- Genbank files are a standardized format containing annotations for genomic sequence information.
- The model file (formatted for matlab)
- The target sequence (FASTA-formatted) contains the sequence where binding sites are sought after

To test this script using data in the repository, try:
python Cas9_Calculator/Cas9_Calculator.py -g Tiplasmidsequence.gb GCA_000219665.2_ASM21966v2_genomic.gbff -m Cas9_Calculator/All_dataModel.mat -t Virb6.fa

This will print all candidate gRNAs that could be used to cut the Virb6 gene,
along with their respective deltaG,
and partition functions against the *Agrobacterium rhizogenes* genome (the organism which we are trying to make less pathogenic!).

### Future work

The output is being printed to standard output (stdout).
We are currently implementing a ranking scheme to sort the candidate guide RNAs by partition function.
These will be written to a file with the fields:
 guideRNA sequence, target sequence, position, dG_Target, % partition function

We should also have a better description for the `.mat` files.

### Contact


email: ubcigem@gmail.com

For more information on out team, please visit our [website](http://www.ubcigem.com/)
