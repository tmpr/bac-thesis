JoinIntoGroupsOf() {
    xargs -n "$@" echo |
        sed 's/ /,/g'
}

python -c "from random import choices; print('\n'.join(str(n) for n in set(choices(range(100_000, 1_000_000), k=10000))));" |
    JoinIntoGroupsOf 5000 |
    xargs -n 1 sh -c \
        'efetch -email "alexandertemper27@gmail.com" -db protein -id $0 -format fasta_cds_na' | tqdm \
    >random.fasta
