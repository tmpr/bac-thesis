
for protein_family_interpro_code in "IPR011014" "IPR023584" "IPR002547"; do
	for organism in "plants" "mammals" "insects" "reptiles" "animals"; do
		raw_location="data/raw/${protein_family_interpro_code}_${organism}.fasta"
		preprocessed_location="data/preprocessed/${protein_family_interpro_code}_${organism}.fasta"
		aligned_location="data/aligned/${protein_family_interpro_code}_${organism}.fasta"

		bash src/get_protein_family.bash $protein_family_interpro_code $organism

		rm data/preprocessed/*

		python src/preprocessing.py \
		    --in_ $raw_location --out data/preprocessed/ \
		    --bin-width=50 --min-bin-size 1

		rm data/aligned/*

		ls -1 data/preprocessed/ |
		    grep $protein_family_interpro_code |
		    xargs -I {} kalign -f fasta -i data/preprocessed/{} -o data/aligned/{} >/dev/null

		python src/blosum.py $protein_family_interpro_code 0.9 $organism 
	done
done

