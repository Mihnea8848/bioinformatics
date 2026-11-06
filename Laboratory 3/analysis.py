import gzip
import io
import time
import zipfile
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path

import pandas as pd
import requests

CACHE = Path(__file__).parent / "cache"
SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/viral/assembly_summary.txt"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USDA_URL = "https://fdc.nal.usda.gov/fdc-datasets/FoodData_Central_sr_legacy_food_csv_2018-04.zip"

BASES = "UCAG"
CODONS = ["".join(c) for c in product(BASES, repeat=3)]
AMINO_ACIDS = dict(zip(CODONS, "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"))
AMINO_ACID_NAMES = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic acid", "C": "Cysteine",
    "Q": "Glutamine", "E": "Glutamic acid", "G": "Glycine", "H": "Histidine", "I": "Isoleucine",
    "L": "Leucine", "K": "Lysine", "M": "Methionine", "F": "Phenylalanine", "P": "Proline",
    "S": "Serine", "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine", "*": "Stop",
}

USDA_AMINO_ACIDS = {
    1210: "Tryptophan", 1211: "Threonine", 1212: "Isoleucine", 1213: "Leucine", 1214: "Lysine",
    1215: "Methionine", 1216: "Cystine", 1217: "Phenylalanine", 1218: "Tyrosine", 1219: "Valine",
    1220: "Arginine", 1221: "Histidine", 1222: "Alanine", 1223: "Aspartic acid", 1224: "Glutamic acid",
    1225: "Glycine", 1226: "Proline", 1227: "Serine",
}
USDA_EQUIVALENT = {"Asparagine": "Aspartic acid", "Glutamine": "Glutamic acid", "Cysteine": "Cystine"}
MIN_PROTEIN = 5


def download(url, path, max_age=None):
    path = CACHE / path
    if not path.exists() or (max_age and time.time() - path.stat().st_mtime > max_age):
        CACHE.mkdir(exist_ok=True)
        response = requests.get(url, timeout=300)
        response.raise_for_status()
        path.write_bytes(response.content)
    return path


@lru_cache
def load_assemblies():
    rows = []
    for line in download(SUMMARY_URL, "assembly_summary.txt", max_age=24 * 3600).read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[19] != "na":
            rows.append((fields[0], fields[7], fields[8] if fields[8] != "na" else "", fields[19]))
    return pd.DataFrame(rows, columns=["accession", "organism", "strain", "ftp"])


def search_entrez(query):
    params = {"db": "assembly", "term": f'"{query}"[Organism] AND latest[filter]', "retmax": 50, "retmode": "json"}
    ids = requests.get(f"{EUTILS}/esearch.fcgi", params=params, timeout=60).json()["esearchresult"]["idlist"]
    if not ids:
        return []
    params = {"db": "assembly", "id": ",".join(ids), "retmode": "json"}
    result = requests.get(f"{EUTILS}/esummary.fcgi", params=params, timeout=60).json()["result"]
    rows = []
    for uid in result["uids"]:
        item = result[uid]
        ftp = item["ftppath_refseq"] or item["ftppath_genbank"]
        if ftp:
            organism = item["organism"].removesuffix(" (viruses)")
            rows.append((item["assemblyaccession"], organism, item["speciesname"], ftp.replace("ftp://", "https://")))
    return rows


def search(query):
    query = query.replace("_", " ").strip()
    assemblies = load_assemblies()
    local = assemblies[assemblies["organism"].str.contains(query, case=False, regex=False)]
    try:
        remote = search_entrez(query)
    except requests.RequestException:
        remote = []
    found = pd.concat([local, pd.DataFrame(remote, columns=local.columns)])
    return found.drop_duplicates("accession").reset_index(drop=True)


def download_genome(accession, ftp):
    name = ftp.rstrip("/").split("/")[-1]
    path = download(f"{ftp}/{name}_genomic.fna.gz", f"{accession}.fna.gz")
    records, current = [], []
    for line in gzip.open(path, "rt"):
        if line.startswith(">"):
            if current:
                records.append("".join(current))
            current = []
        else:
            current.append(line.strip().upper().replace("T", "U"))
    records.append("".join(current))
    return records


def count_codons(records):
    counts = Counter({codon: 0 for codon in CODONS})
    for sequence in records:
        for i in range(len(sequence) - 2):
            codon = sequence[i:i + 3]
            if codon in counts:
                counts[codon] += 1
    return counts


def amino_acid_profile(counts):
    totals = Counter()
    for codon, count in counts.items():
        totals[AMINO_ACID_NAMES[AMINO_ACIDS[codon]]] += count
    return totals


@lru_cache
def load_foods():
    folder = "FoodData_Central_sr_legacy_food_csv_2018-04/"
    with zipfile.ZipFile(download(USDA_URL, "usda_sr_legacy.zip")) as archive:
        foods = pd.read_csv(archive.open(folder + "food.csv"), usecols=["fdc_id", "description", "food_category_id"])
        categories = pd.read_csv(archive.open(folder + "food_category.csv"), usecols=["id", "description"])
        nutrients = pd.read_csv(archive.open(folder + "food_nutrient.csv"), usecols=["fdc_id", "nutrient_id", "amount"])
    categories.columns = ["food_category_id", "category"]
    nutrients = nutrients[nutrients["nutrient_id"].isin(USDA_AMINO_ACIDS)]
    table = nutrients.pivot_table(index="fdc_id", columns="nutrient_id", values="amount").rename(columns=USDA_AMINO_ACIDS)
    return foods.merge(categories, on="food_category_id").set_index("fdc_id").join(table.dropna(), how="inner")


def virus_vector(amino_acids):
    vector = Counter()
    for name, count in amino_acids.items():
        if name != "Stop":
            vector[USDA_EQUIVALENT.get(name, name)] += count
    return pd.Series(vector)[list(USDA_AMINO_ACIDS.values())]


def rank_foods(vector):
    foods = load_foods()
    amounts = foods[list(USDA_AMINO_ACIDS.values())]
    protein = amounts.sum(axis=1)
    amounts = amounts[protein >= MIN_PROTEIN]
    similarity = amounts.dot(vector) / ((amounts**2).sum(axis=1) ** 0.5 * (vector**2).sum() ** 0.5)
    ranked = pd.DataFrame({
        "Food": foods.loc[amounts.index, "description"],
        "Category": foods.loc[amounts.index, "category"],
        "Match %": (similarity * 100).round(2),
        "Protein g/100 g": protein[amounts.index].round(1),
    })
    ranked["group"] = ranked["Food"].str.split(",").str[0]
    ranked = ranked.sort_values("Match %")
    eat = ranked.drop_duplicates("group").head(10).drop(columns="group")
    avoid = ranked.iloc[::-1].drop_duplicates("group").head(10).drop(columns="group")
    return eat, avoid


def analyze(accession, ftp):
    records = download_genome(accession, ftp)
    codons = count_codons(records)
    amino_acids = amino_acid_profile(codons)
    vector = virus_vector(amino_acids)
    eat, avoid = rank_foods(vector)
    return {
        "length": sum(map(len, records)),
        "segments": len(records),
        "codons": codons,
        "amino_acids": amino_acids,
        "favourite": list(vector.sort_values(ascending=False).index[:3]),
        "eat": eat,
        "avoid": avoid,
    }
