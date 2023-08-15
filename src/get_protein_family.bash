# Make sure to have EDirect (https://www.ncbi.nlm.nih.gov/books/NBK179288/) installed.

# Most of this code is essentially taken from https://astrobiomike.github.io/unix/ncbi_eutils.
JoinIntoGroupsOf() {
    xargs -n "$@" echo |
        sed 's/ /,/g'
}

protein_family_interpro_code=$1
accession_file_name="data/accession/${protein_family_interpro_code}"
raw_file_name="data/raw/${protein_family_interpro_code}.fasta"

esearch -db protein -query $protein_family_interpro_code | esummary | xtract -pattern DocumentSummary -element AccessionVersion \
    >$accession_file_name

cat $accession_file_name | JoinIntoGroupsOf 5000 | xargs -n 1 sh -c \
    'efetch -email "alexandertemper27@gmail.com" -db protein -id $0 -format fasta_cds_na' | tqdm \
    >$raw_file_name
