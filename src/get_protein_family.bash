# Make sure to have EDirect (https://www.ncbi.nlm.nih.gov/books/NBK179288/) installed.

# Most of this code is essentially taken from https://astrobiomike.github.io/unix/ncbi_eutils.
JoinIntoGroupsOf() {
    xargs -n "$@" echo |
        sed 's/ /,/g'
}

protein_family_interpro_code="IPR045863"

esearch -db protein -query IPR045863 | esummary | xtract -pattern DocumentSummary -element AccessionVersion >'data/accession/IPR045863'

cat 'data/accession/IPR045863' | JoinIntoGroupsOf 5000 | xargs -n 1 sh -c \
    'efetch -email "alexandertemper27@gmail.com" -db protein -id $0 -format fasta_cds_na' | tqdm \
    >data/raw/IPR045863
