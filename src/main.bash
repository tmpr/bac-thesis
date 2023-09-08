mkdir data blocks matrices
mkdir data/accession data/raw data/preprocessed data/aligned document/plots



for protein_family_interpro_code in "IPR011014" "IPR023584" "IPR002547"; do
	for taxon in "plants" "mammals" "insects" "reptiles" "animals" "bacteria"; do
		raw_location="data/raw/${protein_family_interpro_code}_${taxon}.fasta"
		preprocessed_location="data/preprocessed/${protein_family_interpro_code}_${taxon}.fasta"
		aligned_location="data/aligned/${protein_family_interpro_code}_${taxon}.fasta"

		bash src/get_protein_family.bash $protein_family_interpro_code $taxon

		rm data/preprocessed/*

		python src/preprocessing.py \
		    --in_ $raw_location --out data/preprocessed/ \
		    --bin-width=50 --min-bin-size 1

		rm data/aligned/*

		ls -1 data/preprocessed/ |
		    grep $protein_family_interpro_code |
		    xargs -I {} kalign -f fasta -i data/preprocessed/{} -o data/aligned/{} >/dev/null

		python src/blosum.py $protein_family_interpro_code 0.9 $taxon 
	done
done

python src/create_visualizations.py

cd document
latexmk -pdf -f -interaction=nonstopmode thesis.tex
