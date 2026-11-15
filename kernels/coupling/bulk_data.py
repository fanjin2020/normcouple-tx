"""GSE125583: use deposited count, never VALUE (normalized nRPKM)."""
import gzip
import io
from pathlib import Path
import re
import tarfile
import numpy as np
import pandas as pd
from common import HERE, acquire, require, save_json, read_json, sha

URLS = {
    "GSE125583_RAW.tar": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125583/suppl/GSE125583_RAW.tar",
    "GSE125583_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125583/soft/GSE125583_family.soft.gz",
    "Homo_sapiens.gene_info.gz": "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
}


def parse_soft(path):
    records, row = [], None
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        for line in stream:
            line = line.strip()
            if line.startswith("^SAMPLE = "):
                if row is not None:
                    records.append(row)
                row = {"sample_id": line.split(" = ", 1)[1], "processing": ""}
            elif row is not None and line.startswith("!Sample_characteristics_ch1 = "):
                key, value = line.split(" = ", 1)[1].split(":", 1)
                key = key.strip().lower().replace(" ", "_").replace(".", "_")
                require(key not in row, "Duplicate sample characteristic: " + key)
                row[key] = value.strip()
            elif row is not None and line.startswith("!Sample_data_processing = "):
                row["processing"] += line.split(" = ", 1)[1] + "\n"
            elif row is not None and line.startswith("!Sample_title = "):
                row["title"] = line.split(" = ", 1)[1]
    if row is not None:
        records.append(row)
    meta = pd.DataFrame(records)
    require(len(meta) == 289 and meta.sample_id.is_unique, "GSE125583 must contain 289 uniquely identified samples")
    require(meta.processing.str.contains("counts and nRPKM", case=False, regex=False).all(), "Source no longer documents the count column; inspect manually")
    required = ["case_id", "diagnosis", "age", "sex", "tissue"]
    require(set(required) <= set(meta), "Required clinical fields missing")
    require(meta.case_id.is_unique, "Repeated donors need a different model; refusing independent-sample analysis")
    return meta


def read_counts(path, sample_ids):
    columns, members = {}, []
    with tarfile.open(path, "r") as archive:
        for member in archive:
            require(member.isfile() and not member.issym() and member.size < 5_000_000, "Unexpected archive member")
            require(Path(member.name).name == member.name, "Archive paths must be plain filenames")
            match = re.fullmatch(r"(GSM\d+)_sample\d+\.tsv\.gz", member.name)
            require(match is not None, "Unexpected count filename: " + member.name)
            sid = match.group(1)
            require(sid not in columns, "Duplicate count sample")
            data = archive.extractfile(member).read()
            text = gzip.decompress(data).decode("utf-8")
            require("#count=Number of reads uniquely aligning to gene model" in text, "Count provenance header missing")
            frame = pd.read_csv(io.StringIO(text), sep="\t", comment="#", dtype={"ID_REF": str})
            require(set(frame.columns) == {"ID_REF", "VALUE", "count"}, "Unexpected count table schema")
            require(frame.ID_REF.is_unique and frame.ID_REF.str.fullmatch(r"\d+").all(), "Invalid Entrez gene identities")
            count = pd.to_numeric(frame["count"], errors="raise").to_numpy(float)
            require(np.isfinite(count).all() and (count >= 0).all() and np.array_equal(count, np.rint(count)), "count column is not nonnegative integer read counts")
            columns[sid] = pd.Series(count.astype(np.int64), index=frame.ID_REF)
            members.append({"sample_id": sid, "archive_member": member.name, "genes": len(frame), "library_size": int(count.sum())})
    require(set(columns) == set(sample_ids), "Count sample identities do not match clinical metadata")
    first = columns[sample_ids[0]].index
    require(all(s.index.equals(first) for s in columns.values()), "Gene axes differ across deposited count files")
    counts = pd.DataFrame({sid: columns[sid] for sid in sample_ids})
    require((counts.sum() > 0).all(), "An empty bulk library is present")
    return counts, pd.DataFrame(members)


def gene_mapping(path):
    info = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    info = info[info["#tax_id"].eq("9606")].copy()
    require(len(info)>20000 and info.GeneID.is_unique, "Human gene information is incomplete or duplicated")
    rows = []
    for record in info.itertuples(index=False):
        d = dict(zip(info.columns, record))
        ens = sorted(set(re.findall(r"(?:^|\|)Ensembl:(ENSG\d+)(?=\.|\||$)", d["dbXrefs"])))
        for gene in ens:
            rows.append(dict(entrez=d["GeneID"], ensembl=gene, symbol=d["Symbol"]))
    mapping = pd.DataFrame(rows).drop_duplicates()
    unique_entrez = mapping.groupby("entrez").ensembl.nunique().eq(1)
    unique_ens = mapping.groupby("ensembl").entrez.nunique().eq(1)
    mapping["one_to_one"] = mapping.entrez.map(unique_entrez) & mapping.ensembl.map(unique_ens)
    return mapping


def prepare(args, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    pins = read_json(HERE / "public_sources.json") if (HERE / "public_sources.json").exists() else {}
    snapshot = HERE / "gene_id_mapping_snapshot.csv.gz"
    requested = {k:v for k,v in URLS.items() if k != "Homo_sapiens.gene_info.gz" or not snapshot.exists()}
    paths = {name: acquire(url, args.cache, name, args.offline, pins.get(name)) for name, url in requested.items()}
    meta = parse_soft(paths["GSE125583_family.soft.gz"])
    counts, members = read_counts(paths["GSE125583_RAW.tar"], meta.sample_id.tolist())
    if snapshot.exists():
        annotation = read_json(HERE / "annotation_provenance.json")
        require(sha(snapshot)==annotation["snapshot_sha256"], "Bundled annotation snapshot changed")
        mapping = pd.read_csv(snapshot, dtype={"entrez":str,"ensembl":str,"symbol":str},keep_default_na=False)
        require(mapping.one_to_one.dtype == bool, "Invalid annotation snapshot boolean column")
    else:
        mapping = gene_mapping(paths["Homo_sapiens.gene_info.gz"])
        annotation = dict(source_sha256=sha(paths["Homo_sapiens.gene_info.gz"]),source_url=URLS["Homo_sapiens.gene_info.gz"])
    meta.to_csv(out / "all_sample_metadata.csv", index=False)
    members.to_csv(out / "count_member_inventory.csv", index=False)
    mapping.to_csv(out / "gene_id_mapping.csv", index=False)
    provenance = dict(accession="GSE125583", samples=len(meta), genes=len(counts), measurement="deposited unnormalized gene-level read counts",
        selected_column="count", rejected_normalized_column="VALUE", genome="GRCh38 as documented in GEO", annotation="Deposited Entrez IDs; pinned NCBI one-to-one Ensembl mapping",
        annotation_provenance=annotation,
        sources=[dict(name=n, url=URLS[n], sha256=sha(p), bytes=p.stat().st_size) for n,p in paths.items()])
    save_json(out / "public_input_provenance.json", provenance)
    return counts, meta, mapping, provenance


def cohort(meta):
    meta = meta.copy()
    allowed = {"Alzheimer's disease": "AD", "Alzheimer disease": "AD", "Control": "Control", "control": "Control", "Neurologically normal": "Control", "neurologically normal": "Control"}
    meta["group"] = meta.diagnosis.map(allowed)
    require(meta.group.notna().all(), "Unknown diagnosis labels; do not infer labels from sample ordering: " + str(meta.loc[meta.group.isna(), "diagnosis"].unique()))
    meta["age_numeric"] = pd.to_numeric(meta.age, errors="coerce")
    meta["included"] = meta.age_numeric.between(18, 120) & meta.sex.isin(["F", "M"]) & meta.tissue.str.lower().eq("fusiform gyrus")
    meta["exclusion"] = np.where(meta.included, "", "Missing/unsupported age, sex or tissue")
    m = meta[meta.included].copy().reset_index(drop=True)
    require(m.group.value_counts().min() >= 20 and set(m.group) == {"AD", "Control"}, "Insufficient complete-case disease groups")
    age = (m.age_numeric.to_numpy(float)-m.age_numeric.mean())/m.age_numeric.std(ddof=1)
    cov = np.column_stack([age, m.sex.eq("M").to_numpy(float)])
    x = m.group.eq("AD").to_numpy(float)
    require(np.linalg.matrix_rank(np.column_stack([np.ones(len(m)), x, cov])) == 4, "Disease/age/sex design is rank deficient")
    return m, x, cov, meta
