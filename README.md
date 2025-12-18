# In_House_E.Coli_analysis
We have **1,353** assemblies in our in-house database for E.coli. 
## Search target genes in the in-house assemblies
## Create marker gene file
nano all.gene.fasta
## make db
```
makeblastdb -in all_markers.fasta -dbtype nucl -out all_markers
```
## Run BLAST against all the assemblies
```
mkdir -p blast_out
mkdir -p blast_hits_only

for f in $(find /home/jing/E.coli_test/ecoli_all -name "*.fasta"); do
  folder=$(basename "$(dirname "$f")")
  base=$(basename "$f" .fna)
  out_file="blast_out/${folder}_${base}_hits.tsv"

  echo "Running BLAST on $folder ..."

  blastn -query all_markers.fasta \
         -subject "$f" \
         -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
         | awk '{cov=($4/$5)*100; if($3>=90 && cov>=90) print $0}' > "$out_file"

  # If file is non-empty (has hits), copy to hits-only folder
  if [ -s "$out_file" ]; then
      cp "$out_file" blast_hits_only/
      echo "✅ Hits found in $folder — copied to blast_hits_only/"
  else
      echo "❌ No hits for $folder"
  fi
done
```
## Number of isolates that have at least one of the marker genes
```
cd blast_hits_only
ls -l . | wc -l
# 1204
```
<img width="672" height="35" alt="image" src="https://github.com/user-attachments/assets/9efbc6a1-522c-4a6a-b18c-8033f5c8ef80" />

## Make presence and absence martrix
```
#!/usr/bin/env python3
import argparse
import glob
import os
import re
import pandas as pd

# ------------------------------------------------------------
# Sample ID extraction
# ------------------------------------------------------------
def extract_sample_id(path: str) -> str:
    """
    Extract sample ID like W515056 from BLAST TSV filename.

    Example:
      ecoli_all_W515056.contigs.fasta_hits.tsv  -> W515056
      ecoli_all_W515056.contigs.fasta.tsv       -> W515056
      W515056_hits.tsv                          -> W515056

    Fallback:
      strips common suffixes/extensions and returns remaining basename.
    """
    base = os.path.basename(path)

    m = re.search(r"(W\d+)", base)
    if m:
        return m.group(1)

    # Fallback stripping
    name = base
    for suf in ["_hits.tsv", ".hits.tsv", ".tsv"]:
        if name.endswith(suf):
            name = name[:-len(suf)]
            break
    # remove remaining extensions (e.g. .fasta, .contigs, etc.)
    name = re.sub(r"\.(fasta|fa|fna|contigs)(\..*)?$", "", name)
    return name


# ------------------------------------------------------------
# BLAST TSV reader (supports headerless outfmt 6)
# ------------------------------------------------------------
def read_blast_tsv(path: str) -> pd.DataFrame:
    cols11 = ["qseqid","sseqid","pident","length","qlen",
              "qstart","qend","sstart","send","evalue","bitscore"]
    cols10 = ["qseqid","sseqid","pident","length",
              "qstart","qend","sstart","send","evalue","bitscore"]

    # Try reading as headered
    try:
        df = pd.read_csv(path, sep="\t")
        if "qseqid" not in df.columns:
            raise ValueError("No header detected")
    except Exception:
        first = open(path).readline().rstrip("\n")
        if not first:
            return pd.DataFrame(columns=cols11)
        nf = len(first.split("\t"))
        if nf == 11:
            df = pd.read_csv(path, sep="\t", header=None, names=cols11)
        elif nf == 10:
            df = pd.read_csv(path, sep="\t", header=None, names=cols10)
        else:
            return pd.DataFrame(columns=cols11)

    # Ensure 11-column schema exists
    for c in cols11:
        if c not in df.columns:
            df[c] = pd.NA
    return df


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description="Build presence/absence matrix + best-hit table from BLAST outfmt 6 TSVs."
    )
    ap.add_argument("--glob", default="blast_out/*_hits.tsv",
                    help='Glob for BLAST TSVs (default: blast_out/*_hits.tsv)')
    ap.add_argument("--id", type=float, default=80.0,
                    help="Min % identity to keep (default: 80)")
    ap.add_argument("--cov", type=float, default=0.80,
                    help="Min query coverage (length/qlen) to keep when qlen exists (default: 0.80)")
    ap.add_argument("--outprefix", default="blast",
                    help="Output prefix (default: blast)")
    args = ap.parse_args()

    files = sorted(glob.glob(args.glob, recursive=True))
    if not files:
        raise SystemExit(f"❌ No files matched: {args.glob}")

    presence_rows = []
    best_rows = []

    for f in files:
        sample = extract_sample_id(f)
        df = read_blast_tsv(f)

        # empty file -> all zeros later
        if df.empty:
            presence_rows.append({"Sample": sample})
            continue

        # Filter by identity and (if available) coverage
        df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
        df["length"] = pd.to_numeric(df["length"], errors="coerce")
        df["qlen"]   = pd.to_numeric(df["qlen"], errors="coerce")

        if df["qlen"].notna().any():
            df["coverage"] = df["length"] / df["qlen"]
            df = df[(df["pident"] >= args.id) & (df["coverage"] >= args.cov)]
        else:
            df["coverage"] = 1.0
            df = df[df["pident"] >= args.id]

        if df.empty:
            presence_rows.append({"Sample": sample})
            continue

        # Use qseqid as-is (no dropping unknown)
        df["Gene"] = df["qseqid"].astype(str).str.strip().str.split().str[0]

        # Best hit per gene: max identity, then coverage, then length
        best = (
            df.sort_values(["Gene", "pident", "coverage", "length"],
                           ascending=[True, False, False, False])
              .groupby("Gene", as_index=False)
              .first()
        )

        # presence row
        row = {"Sample": sample}
        for g in best["Gene"].unique():
            row[g] = 1
        presence_rows.append(row)

        # best-hit rows
        for _, r in best.iterrows():
            best_rows.append({
                "Sample": sample,
                "Gene": r["Gene"],
                "Identity(%)": round(float(r["pident"]), 2) if pd.notna(r["pident"]) else None,
                "Coverage": round(float(r["coverage"]), 3) if pd.notna(r["coverage"]) else None,
                "Subject": r.get("sseqid", ""),
                "AlignLen": int(r["length"]) if pd.notna(r["length"]) else None,
                "Evalue": r.get("evalue", pd.NA),
                "Bitscore": r.get("bitscore", pd.NA),
            })

    # Build matrix and fill missing with 0
    mat = pd.DataFrame(presence_rows).fillna(0)

    # Ensure Sample first, keep all other columns as they appear
    cols = list(mat.columns)
    if "Sample" in cols:
        cols.remove("Sample")
        mat = mat[["Sample"] + cols]

    # Convert presence columns to int
    for c in mat.columns:
        if c != "Sample":
            mat[c] = mat[c].astype(int)

    pa_path = f"{args.outprefix}_presence_absence.tsv"
    bh_path = f"{args.outprefix}_best_hits.tsv"

    mat.to_csv(pa_path, sep="\t", index=False)
    pd.DataFrame(best_rows).to_csv(bh_path, sep="\t", index=False)

    print(f"✅ Wrote: {pa_path}")
    print(f"✅ Wrote: {bh_path}")
    print("Tip: view matrix: column -t <file> | less -S")


if __name__ == "__main__":
    main()

```
## Run the script
```
python3 1.build_presence_absence.py \
  --glob "blast_hits_only/*_hits.tsv" \
  --outprefix blast

```
# Classification EHEC or EPEC or STEC
```
#!/usr/bin/env python3
import re
import pandas as pd

IN = "blast_presence_absence.tsv"

# ---------------------------
# Helpers
# ---------------------------
def extract_id(s):
    """
    Your Sample values look like:
      ecoli_all_W515056
      ecoli_all_IEH-NGS-ECO-00478
      W507963
    We want to keep the meaningful isolate ID:
      W515056 / IEH-NGS-ECO-00478 / W507963
    If a GCA/GCF is present, keep that instead.
    """
    if pd.isna(s):
        return None
    s = str(s).strip()

    # If GCA/GCF present, prioritize it
    m = re.search(r"(GC[AF]_\d+\.\d+)", s)
    if m:
        return m.group(1)

    # Remove any path
    s = s.split("/")[-1]

    # Drop common prefixes like "ecoli_all_"
    s = re.sub(r"^ecoli_all[_\-]+", "", s)

    # Drop common suffixes
    s = re.sub(r"(\.contigs\.fasta.*|\.fasta.*|\.fa.*|\.fna.*|_hits.*|_genomic.*)$", "", s)

    return s.strip()

def ensure_binary(df, col):
    if col not in df.columns:
        raise ValueError(f"Missing column: {col}")
    df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    # anything >0 becomes 1
    df[col] = (df[col] > 0).astype(int)

# ---------------------------
# 1) Load
# ---------------------------
df = pd.read_csv(IN, sep="\t")

# normalize Sample column name
if "Sample" not in df.columns:
    cand = [c for c in df.columns if c.lower() == "sample"]
    if not cand:
        raise ValueError("No 'Sample' column found (case-insensitive).")
    df = df.rename(columns={cand[0]: "Sample"})

# create cleaned isolate id
df["sample"] = df["Sample"].apply(extract_id)

# ---------------------------
# 2) Ensure markers exist & are binary
# ---------------------------
for col in ["stx1", "stx2", "eae"]:
    ensure_binary(df, col)

# ---------------------------
# 3) Classify pathovar
# ---------------------------
df["pathovar"] = "Other"

# STEC: eae- AND (stx1+ OR stx2+)
df.loc[(df["eae"] == 0) & ((df["stx1"] == 1) | (df["stx2"] == 1)), "pathovar"] = "STEC"

# EHEC: eae+ AND (stx1+ OR stx2+)
df.loc[(df["eae"] == 1) & ((df["stx1"] == 1) | (df["stx2"] == 1)), "pathovar"] = "EHEC"

# EPEC: eae+ AND stx1- AND stx2-
df.loc[(df["eae"] == 1) & (df["stx1"] == 0) & (df["stx2"] == 0), "pathovar"] = "EPEC"

# ---------------------------
# 4) Export ID lists
# ---------------------------
ehec_list = df.loc[df["pathovar"] == "EHEC", "sample"].dropna().drop_duplicates()
epec_list = df.loc[df["pathovar"] == "EPEC", "sample"].dropna().drop_duplicates()
stec_list = df.loc[df["pathovar"] == "STEC", "sample"].dropna().drop_duplicates()

ehec_list.to_csv("EHEC_list.txt", index=False, header=False)
epec_list.to_csv("EPEC_list.txt", index=False, header=False)
stec_list.to_csv("STEC_list.txt", index=False, header=False)

# ---------------------------
# 5) Save classification table (keep original Sample too)
# ---------------------------
out_class = "pathovar_classification.tsv"
df[["Sample", "sample", "pathovar", "stx1", "stx2", "eae"]].to_csv(out_class, sep="\t", index=False)

# ---------------------------
# 6) Print summary
# ---------------------------
print(f"✅ Classification complete from {IN}")
print(f"  EHEC: {len(ehec_list)} isolates → EHEC_list.txt")
print(f"  EPEC: {len(epec_list)} isolates → EPEC_list.txt")
print(f"  STEC: {len(stec_list)} isolates → STEC_list.txt")
print(f"  Saved summary: {out_class}")

```
<img width="549" height="66" alt="image" src="https://github.com/user-attachments/assets/08c59201-6e1c-4e4a-8bbf-5eb017602cef" />

## Separate presence matrix by pathovar
```
#!/usr/bin/env python3
import pandas as pd

PA = "blast_presence_absence.tsv"
CLASS = "pathovar_classification.tsv"

# load
pa = pd.read_csv(PA, sep="\t")
cls = pd.read_csv(CLASS, sep="\t")

# ensure columns exist
if "Sample" not in pa.columns:
    cand = [c for c in pa.columns if c.lower() == "sample"]
    if not cand:
        raise ValueError("No Sample column in blast_presence_absence.tsv")
    pa = pa.rename(columns={cand[0]: "Sample"})

# merge classification onto presence/absence table
m = pa.merge(cls[["Sample", "pathovar"]], on="Sample", how="left")

# split + save (keep rows even if pathovar missing -> goes to Other)
m[m["pathovar"] == "EHEC"].drop(columns=["pathovar"]).to_csv("blast_presence_absence_EHEC.tsv", sep="\t", index=False)
m[m["pathovar"] == "EPEC"].drop(columns=["pathovar"]).to_csv("blast_presence_absence_EPEC.tsv", sep="\t", index=False)
m[m["pathovar"] == "STEC"].drop(columns=["pathovar"]).to_csv("blast_presence_absence_STEC.tsv", sep="\t", index=False)

print("✅ Wrote:")
print("  blast_presence_absence_EHEC.tsv")
print("  blast_presence_absence_EPEC.tsv")
print("  blast_presence_absence_STEC.tsv")
```

## Count all espK/espV/espN combinations with in EHEC
count_EHEC_esp.py
```
#!/usr/bin/env python3
import pandas as pd

# 1. Load your EPEC-only table
df = pd.read_csv("blast_presence_absence_EHEC.tsv", sep="\t")

# 2. Boolean presence for each gene
K = df["espK"] == 1
V = df["espV"] == 1
N = df["espN"] == 1

total = len(df)

rows = []

def add(label, mask):
    count = int(mask.sum())
    pct = round(count / total * 100, 2) if total > 0 else 0.0
    rows.append({"Pattern": label, "Count": count, "Percent(%)": pct})

# 3. Individual genes
add("espK(+)", K)
add("espV(+)", V)
add("espN(+)", N)

# 4. Simple combinations (AND / OR)
add("espK(+) AND espV(+)", K & V)
add("espK(+) AND espN(+)", K & N)
add("espV(+) AND espN(+)", V & N)
add("espK(+) AND espV(+) AND espN(+)", K & V & N)
add("espK(+) OR espV(+) OR espN(+)", K | V | N)

# 5. 8 mutually exclusive patterns (K/V/N all combos)
onlyK =  K & ~V & ~N
onlyV = ~K &  V & ~N
onlyN = ~K & ~V &  N
KV    =  K &  V & ~N
KN    =  K & ~V &  N
VN    = ~K &  V &  N
KVN   =  K &  V &  N
none  = ~K & ~V & ~N

add("espK only (K+ V- N-)", onlyK)
add("espV only (K- V+ N-)", onlyV)
add("espN only (K- V- N+)", onlyN)
add("espK & espV (K+ V+ N-)", KV)
add("espK & espN (K+ V- N+)", KN)
add("espV & espN (K- V+ N+)", VN)
add("espK & espV & espN (K+ V+ N+)", KVN)
add("none esp (K- V- N-)", none)

# 6. Make DataFrame and save / print
summary = pd.DataFrame(rows)

print(f"Total EPEC isolates: {total}\n")
print(summary)

summary.to_csv("EPEC_esp_counts_percentages.csv", index=False)
print("\nSaved to: EPEC_esp_counts_percentages.csv")
```
