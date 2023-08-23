# Deriving BLOSUM-like scoring matrices for nucleotides

To get this started, create the environment

```conda env create -f environment.yml```

and install [kalign](https://github.com/timolassmann/kalign) and [edirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/).

To get a BLOSUM matrix for a protein family's genes, call 
`bash src/main.bash <INTERPRO_CODE> <X>`, where `<INTERPRO_CODE>`
is the code for the family (in the case of [this family](https://www.ebi.ac.uk/interpro/entry/pfam/PF13603/), it is `PF13603`) and `<X>` is the similarity threshold for clustering the sequences (for the BLOSUM62 matrix, that is 0.62)).