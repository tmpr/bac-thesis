from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.io as pio
from numpy.linalg import norm
from numpy.typing import NDArray
from sklearn.decomposition import PCA

# This is needed as otherwise, the plot ends up having a bug,
# see https://github.com/plotly/plotly.py/issues/3469
pio.kaleido.scope.mathjax = None


def parse_sw_matrix(txt: str) -> NDArray[np.int16]:
    rows = txt.splitlines()
    matrix = np.array(
        [[float(n) for n in row.split()[1:]] for row in rows[1:]], dtype=np.int16
    )
    return matrix


matrices = []
paths = list(Path("matrices").glob("*.matrix"))
paths = [path for path in paths if len(path.stem.split("_")) == 3]
blosumxs, protein_families, taxa, n_sequencess = [], [], [], []
for path in paths:
    if "nan" not in path.read_text():
        txt = path.read_text()
        matrix = parse_sw_matrix(txt)
        matrices.append(matrix)
        blosumx, protein_family, taxon = path.stem.split("_")
        blosumxs.append(blosumx)
        protein_families.append(protein_family)
        taxa.append(taxon)
        n_sequences = (
            Path(f"data/raw/{protein_family}_{taxon}.fasta").read_text().count(">")
        )
        n_sequencess.append(n_sequences)

vectors = np.array([matrix.reshape(-1) for matrix in matrices])
downprojected = PCA(n_components=2).fit_transform(vectors)

sizes = [15 for _ in vectors]
fig = px.scatter(
    x=downprojected[:, 0],
    y=downprojected[:, 1],
    color=taxa,
    symbol=protein_families,
    size=sizes,
    size_max=15,
    template="none",
)
fig.update_layout(
    {"legend_title": "Taxon, Interpro code", "xaxis_title": "", "yaxis_title": ""}
)
fig.write_image("document/plots/downprojection.pdf")

df = pd.DataFrame.from_records(
    list(
        zip(
            taxa,
            protein_families,
            n_sequencess,
            norm(vectors, axis=1),
            vectors.std(axis=1),
        )
    ),
    columns=[
        "Taxon",
        "Code",
        r"\# Seq.",
        r"$\sigma^2(\mathbf S)$",
        r"$\lVert \mathbf S \rVert_2$",
    ],
)
df = df.sort_values(by=r"$\sigma^2(\mathbf S)$")
df = df.round(decimals=2)
df.to_latex("document/plots/table.tex", index=False, escape=False)

nucs = list("ACGT")
df = pd.concat(
    [pd.DataFrame(m, columns=nucs, index=nucs) for m in matrices],
    keys=zip(taxa, protein_families),
)
df = df.sort_index()
df.to_latex("document/plots/matrices.tex", multirow=True)


def protein_family_list(protein_family_codes: List[str]) -> str:
    lines = [r"\begin{itemize}"]
    for protein_family_code in set(protein_family_codes):
        db = "pfam" if protein_family_code.startswith("PF") else "interpro"
        lines.append(
            f"\t \\item \\href{{https://www.ebi.ac.uk/interpro/entry/{db}/{protein_family_code}}}{{{protein_family_code}}}"
        )
    lines.append(r"\end{itemize}")
    return "\n".join(lines)


Path("document/plots/protein_families.tex").write_text(
    protein_family_list(protein_families)
)
