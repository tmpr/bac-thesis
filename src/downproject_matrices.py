from pathlib import Path

import numpy as np
from numpy.typing import NDArray
import plotly.express as px
from sklearn.decomposition import PCA
import pandas as pd
import plotly.io as pio
from numpy.linalg import norm

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
blosumxs, protein_families, organisms = [], [], []
for path in paths:
    if "nan" not in path.read_text():
        txt = path.read_text()
        matrix = parse_sw_matrix(txt)
        matrices.append(matrix)
        blosumx, protein_family, organism = path.stem.split("_")
        blosumxs.append(blosumx)
        protein_families.append(protein_family)
        organisms.append(organism)

vectors = np.array([matrix.reshape(-1) for matrix in matrices])
downprojected = PCA(n_components=2).fit_transform(vectors)

sizes = [15 for _ in vectors]
fig = px.scatter(
    x=downprojected[:, 0],
    y=downprojected[:, 1],
    color=organisms,
    symbol=protein_families,
    size=sizes,
    size_max=15,
)
fig.update_layout(
    {"legend_title": "Organisms, Interpro code", "xaxis_title": "", "yaxis_title": ""}
)
fig.write_image("document/plots/downprojection.pdf")


def bmatrix(a):
    """Returns a LaTeX bmatrix

    :a: numpy array
    :returns: LaTeX bmatrix as a string
    """
    if len(a.shape) > 2:
        raise ValueError("bmatrix can at most display two dimensions")
    lines = str(a).replace("[", "").replace("]", "").splitlines()
    rv = [r"\begin{mbatrix}"]
    rv += r"\\".join("  " + " & ".join(line.split()) for line in lines)
    rv += [r"\end{matrix}"]
    return f"${''.join(rv)}$"


df = pd.DataFrame.from_records(
    list(
        zip(
            organisms,
            protein_families,
            norm(vectors, axis=1),
            vectors.std(axis=1),
        )
    ),
    columns=[
        "Organisms",
        "Interpro code",
        r"\sigma^2(\mathbf S)",
        r"\lVert \mathbf S \rVert_2",
    ],
)
df = df.sort_values(by=r"\sigma^2(\mathbf S)")
df = df.round(decimals=2)
df.to_latex("document/plots/table.tex", bold_rows=True, escape=False)
