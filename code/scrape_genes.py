import json
import multiprocessing
import re
from multiprocessing import Pool
from time import sleep
from typing import Sequence, Tuple

import requests
from tqdm import trange


def scrape_fastas(
    position_and_refseq_nucleotide_ids: Tuple[int, Sequence[str]],
    batch_size: int = 500,
) -> list[str]:
    position, refseq_nucleotide_ids = position_and_refseq_nucleotide_ids
    sleep(position)
    base = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    ids = [("id", id_) for id_ in refseq_nucleotide_ids]
    payload = [
        ("db", "nuccore"),
        ("retmax", batch_size),
        ("usehistory", "y"),
        ("rettype", "fasta"),
    ] + ids

    fastas = []
    for retstart in trange(
        0, len(refseq_nucleotide_ids), batch_size, position=position
    ):
        fastas.append(
            requests.post(base, data=payload + [("retstart", retstart)]).text
        )
    return fastas


def scrape_genbank_ids(interpro_id: str) -> set[str]:
    response = requests.get(
        f"https://rest.uniprot.org/uniprotkb/stream?&format=json&query=%28{interpro_id}%29"
    )
    data = json.loads(response.text)
    refseq_nucleotide_ids = set()
    protein_ids = set()

    interpro_id = transmembrane_code
    for result in data["results"]:
        is_correct_interpro_id = False
        gene_id = None

        for reference in result["uniProtKBCrossReferences"]:
            if reference["database"] == "GeneID":
                protein_ids.add(reference["id"])
                for property in reference["properties"]:
                    if property["key"] == "NucleotideSequenceId":
                        gene_id = property["value"]
            if reference["database"] == "InterPro":
                if reference["id"] == interpro_id:
                    is_correct_interpro_id = True

        if is_correct_interpro_id and gene_id is not None:
            refseq_nucleotide_ids.add(gene_id)

    return refseq_nucleotide_ids


if __name__ == "__main__":
    accessions = []
    transmembrane_code = "IPR002547"

    genbank_ids = list(scrape_genbank_ids(transmembrane_code))

    refseq_nucleotide_id = list(genbank_ids)[0]
    saved_fastas = []

    # NCBI etools allow at max 3 requests per second,
    # so we have each worker start in a different second.
    with Pool() as p:
        dataset_part_size = len(genbank_ids) // multiprocessing.cpu_count()
        genbank_batches = [
            genbank_ids[i : i + dataset_part_size]
            for i in range(0, len(genbank_ids), dataset_part_size)
        ]
        with open(f"{transmembrane_code}.fasta", "w") as f:
            for fastas in p.imap_unordered(
                scrape_fastas, enumerate(genbank_batches, 1)
            ):
                for fasta in fastas:
                    saved_fastas.append(fasta)
                    f.write(fasta + "\n")

    print(len(saved_fastas))

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?db=nuccore&dbfrom=protein&linkname=protein_nuccore&id=NP_418260.1&cmd=neighbor

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=WP_000387388.1&rettype=gp&retmode=xml

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=75204808&rettype=xml

# https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=21614549&strand=1&seq_start=1,2&seq_stop=100,101&id=1866753945&rettype=fasta&retmode=text
