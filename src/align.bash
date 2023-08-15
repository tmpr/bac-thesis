ls -1 data/preprocessed | xargs -I {} kalign -f fasta -i data/preprocessed/{} -o data/aligned/{}
