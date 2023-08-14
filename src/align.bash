ls -1 preprocessed | xargs -I {} kalign -f fasta -i preprocessed/{} -o {}
