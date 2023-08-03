# Make sure to have EDirect (https://www.ncbi.nlm.nih.gov/books/NBK179288/) installed.

# Most of this code is essentially taken from https://astrobiomike.github.io/unix/ncbi_eutils.
JoinIntoGroupsOf() {
    xargs -n "$@" echo |
        sed 's/ /,/g'
}

protein_family_interpro_code="IPR002547"

esearch -db protein -query IPR002547 | esummary | xtract -pattern DocumentSummary -element AccessionVersion >'IPR002547.accessions.txt'

cat 'IPR002547.accessions.txt' | JoinIntoGroupsOf 5000 | xargs -n 1 sh -c \
    'efetch -email "alexandertemper27@gmail.com" -db protein -id $0 -format fasta_cds_na' | tqdm \
    >IPR002547.faa
