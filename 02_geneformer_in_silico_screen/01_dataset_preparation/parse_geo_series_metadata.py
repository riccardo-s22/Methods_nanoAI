#!/usr/bin/env python
"""Parse GSE209552_family.soft into a per-sample TSV with curated columns."""
import re
from pathlib import Path
import pandas as pd

BASE = Path("/mnt/d/TBI_SC/Geneformer_TBI")
SOFT = Path("/mnt/d/TBI_SC/GSE209552_family.soft/GSE209552_family.soft")
SRA = Path("/mnt/d/TBI_SC/SraRunTable (1).csv")
OUT = BASE / "metadata/gse209552_sample_metadata_curated.tsv"

SUFFIX_TO_REGION = {"F": "frontal", "T": "temporal", "LT": "temporal_L", "RT": "temporal_R",
                    "LF": "frontal_L", "RF": "frontal_R", "LFP": "frontoparietal_L",
                    "IFP": "inferior_frontoparietal", "LIFP": "inferior_frontoparietal_L"}

rows = []
cur = None
with open(SOFT) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE = "):
            if cur:
                rows.append(cur)
            cur = {"sample": line.split("= ", 1)[1].strip(), "characteristics": []}
        elif cur is None:
            continue
        elif line.startswith("!Sample_title = "):
            cur["title"] = line.split("= ", 1)[1].strip()
        elif line.startswith("!Sample_source_name_ch1 = "):
            cur["source_name"] = line.split("= ", 1)[1].strip()
        elif line.startswith("!Sample_library_strategy = "):
            cur["library_strategy"] = line.split("= ", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1 = "):
            cur["characteristics"].append(line.split("= ", 1)[1].strip())
        elif line.startswith("!Sample_description = "):
            cur.setdefault("descriptions", []).append(line.split("= ", 1)[1].strip())
if cur:
    rows.append(cur)

for r in rows:
    for kv in r.pop("characteristics", []):
        if ":" in kv:
            k, v = kv.split(":", 1)
            r[k.strip().lower().replace(" ", "_")] = v.strip()
    r["descriptions"] = "; ".join(r.get("descriptions", []))

md = pd.DataFrame(rows).fillna("")
md["condition"] = md["condition"].astype(str).replace({"tbi": "TBI", "control": "control"})

# Infer assay
def infer_assay(t):
    t = str(t).lower()
    if "scrnaseq" in t or "snrna" in t:
        return "snRNAseq_10x"
    if "bulk rnaseq" in t:
        return "bulk_RNAseq"
    if "hgpc" in t:
        return "hGPC_bulk_RNAseq"
    return "unknown"
md["assay"] = md["title"].map(infer_assay)

# donor + region from title and source-file naming.
def parse_donor_region(row):
    title = str(row.get("title", ""))
    m_ctrl = re.match(r"Control\s+(\d+),\s+(\w+)\s+lobe", title)
    if m_ctrl:
        return m_ctrl.group(1), m_ctrl.group(2).lower()
    m_tbi = re.match(r"TBI\s+(\d+)", title)
    if m_tbi:
        donor = m_tbi.group(1)
        for suf in sorted(SUFFIX_TO_REGION, key=len, reverse=True):
            if re.search(rf"_{suf}_", " ".join(row.get("descriptions", "").split()) + " "):
                return donor, SUFFIX_TO_REGION[suf]
        return donor, ""
    return "", ""

md[["donor", "region"]] = md.apply(lambda r: pd.Series(parse_donor_region(r)), axis=1)

# Add region from descriptions for TBI samples by file suffix _nrN_<SUF>_
def region_from_descriptions(row):
    desc = " ".join(row.get("descriptions", "").split("; "))
    m = re.search(r"_nr\d+_([A-Z]+)_", desc)
    if m:
        suf = m.group(1)
        return SUFFIX_TO_REGION.get(suf, suf)
    return row.get("region", "")
md["region"] = md.apply(region_from_descriptions, axis=1).replace("", pd.NA).fillna(md["region"])

# Pull time post injury if available
md["time_post_injury"] = md.get("time_post_injury", "")

# Pull from SRA Run Table to enrich
if SRA.exists():
    sra = pd.read_csv(SRA)
    sra_md = sra.groupby("Sample Name").first().reset_index()[["Sample Name", "AGE", "sex", "time_post_injury", "cell_type", "treatment"]]
    sra_md = sra_md.rename(columns={"Sample Name": "sample"})
    sra_md = sra_md.add_suffix("_sra")
    sra_md = sra_md.rename(columns={"sample_sra": "sample"})
    md = md.merge(sra_md, on="sample", how="left")
    for col, sra_col in [("age", "AGE_sra"), ("sex", "sex_sra"), ("time_post_injury", "time_post_injury_sra")]:
        if col not in md.columns:
            md[col] = ""
        md[col] = md[col].astype(str).where(md[col].astype(str) != "", md[sra_col].fillna("").astype(str))

# Final cleanup: keep main columns up front
front = ["sample", "title", "condition", "assay", "donor", "region", "age", "sex", "time_post_injury", "source_name", "library_strategy", "descriptions"]
others = [c for c in md.columns if c not in front]
md = md[front + others]

OUT.parent.mkdir(parents=True, exist_ok=True)
md.to_csv(OUT, sep="\t", index=False)
print(f"Wrote {OUT}")
print(md[["sample", "title", "condition", "assay", "donor", "region", "age", "sex"]].to_string(index=False))
