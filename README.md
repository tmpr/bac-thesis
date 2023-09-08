# Deriving BLOSUM-like scoring matrices for nucleotides

Install the following:

- [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html),
- [kalign](https://github.com/timolassmann/kalign),
- [eDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/),
- [latex](https://www.latex-project.org/get/),
- [latexmk](https://mg.readthedocs.io/latexmk.html).

Create the environment with

```
conda env create -f environment.yml
conda activate nblosum
```

Afterwards, to reproduce the paper and its results run

```
bash src/main.bash
```

Note that the NCBI database is subject to changes, thus full reproducibility like this is not guaranteed.
The pdf in the root level has been created with the data stored in `data.zip`.
