from pathlib import Path

import numpy as np
from numpy.linalg import norm
from numpy.typing import NDArray
import plotly.express as px
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def parse_sw_matrix(txt: str) -> NDArray[np.int16]:
    rows = txt.splitlines()
    matrix = np.array([[float(n) for n in row.split()[1:]] for row in rows[1:]])
    return matrix


matrices = []
paths = list(Path("matrices").glob("*.matrix"))
paths = [path for path in paths if len(path.stem.split("_")) == 3]
blosumxs, protein_families, organisms = [], [], []
for path in paths:
    if "nan" not in path.read_text():
        matrix = parse_sw_matrix(path.read_text())
        matrices.append(matrix)
        blosumx, protein_family, organism = path.stem.split("_")
        blosumxs.append(blosumx)
        protein_families.append(protein_family)
        organisms.append(organism)

vectors = np.array([matrix.reshape(-1) for matrix in matrices])
downprojected = PCA(n_components=2).fit_transform(vectors)

sizes = [norm(vec) for vec in vectors]
fig = px.scatter(x=downprojected[:, 0], y=downprojected[:, 1], color=organisms, symbol=protein_families, size=sizes)
fig.write_image("document/plots/downprojection.pdf")
