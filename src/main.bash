protein_family_interpro_code=$1
raw_location="data/raw/${protein_family_interpro_code}.fasta"
preprocessed_location="data/preprocessed/${protein_family_interpro_code}.fasta"
aligned_location="data/aligned/${protein_family_interpro_code}.fasta"

bash src/get_protein_family.bash $protein_family_interpro_code

rm data/preprocessed/*

python src/preprocessing.py \
    --in_ $raw_location --out data/preprocessed/ \
    --bin-width=50 --min-bin-size 1

ls -1 data/preprocessed/ |
    grep $protein_family_interpro_code |
    xargs -I {} kalign -f fasta -i data/preprocessed/{} -o data/aligned/{} >/dev/null

python src/blosum.py $protein_family_interpro_code $2
