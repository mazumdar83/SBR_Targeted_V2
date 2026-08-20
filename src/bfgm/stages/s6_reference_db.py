"""Stage 6: compile one searchable reference from every sequence the run produced.

The run generates sequences from two directions with incompatible headers:

* ``proteins.fasta``      UniProt, gene-first  ``>sp|P06971|FHUA_ECOLI Ferrichrome...``
* ``ko_sequences.fasta``  KEGG, KO-first       ``>eco:b0150 K02014 iron complex...``

Neither is usable as a BLAST subject on its own: the identifiers do not share a
namespace, the same protein appears in both, and neither header carries enough to trace
a hit back to a KO and a function. This stage merges them into a single deduplicated
FASTA with parseable headers, plus a mapping table so any hit resolves to KO, gene,
organism, and definition.

Header format (BLAST-safe: no whitespace in the ID field)
--------------------------------------------------------
    >bfgm|<n>|<KO>|<gene>|<source_acc> <description> [<organism>]

The ID splits on ``|`` into five fields, so ``blastp -outfmt 6`` gives a ``sseqid`` that
already contains the KO and the gene symbol without a join. The mapping TSV is there for
everything else.
"""
from __future__ import annotations

import hashlib
import re
import shutil
import subprocess
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

import pandas as pd

SAFE = re.compile(r"[^A-Za-z0-9._-]")


def _clean(s: str, fallback: str = "NA") -> str:
    s = SAFE.sub("_", str(s or "").strip())
    return s or fallback


def parse_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    if not path.exists():
        return
    header, chunks = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            header, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def _uniprot_records(run: Path) -> List[Dict]:
    """UniProt side. Prefers the anchored TSV so KO assignments come along."""
    tsv = run / "proteins_anchored.tsv"
    if not tsv.exists():
        tsv = run / "proteins.tsv"
    if not tsv.exists():
        return []
    df = pd.read_csv(tsv, sep="\t", low_memory=False)
    out = []
    for _, r in df.iterrows():
        seq = str(r.get("Sequence", "") or "")
        if not seq:
            continue
        kos = [k for k in re.findall(r"K\d{5}", str(r.get("kegg_KOs", "")))]
        gene = str(r.get("Gene Names", "") or "").split()
        out.append({
            "seq": seq,
            "ko": kos[0] if kos else "noKO",
            "all_kos": ";".join(kos),
            "gene": gene[0] if gene else str(r.get("matched_query_terms", "NA")).split(";")[0],
            "source": "uniprot",
            "source_acc": str(r.get("Entry", "NA")),
            "organism": str(r.get("Organism", "NA")),
            "taxid": str(r.get("Organism (ID)", "")),
            "description": str(r.get("Protein names", "") or "")[:140],
        })
    return out


def _kegg_records(run: Path) -> List[Dict]:
    """KEGG side. The KO and description sit in the aaseq header."""
    out = []
    for header, seq in parse_fasta(run / "ko_sequences.fasta"):
        if not seq:
            continue
        gene_id = header.split()[0]
        m = re.search(r"\b(K\d{5})\b", header)
        desc = header.split("|")[-1].strip() if "|" in header else header
        gene = "NA"
        gm = re.search(r"\(RefSeq\)\s*([A-Za-z0-9_]+)", header)
        if gm:
            gene = gm.group(1)
        out.append({
            "seq": seq,
            "ko": m.group(1) if m else "noKO",
            "all_kos": m.group(1) if m else "",
            "gene": gene,
            "source": "kegg",
            "source_acc": gene_id,
            "organism": gene_id.split(":")[0],
            "taxid": "",
            "description": desc[:140],
        })
    return out


def run(run_dir: str | Path, make_blast_db: bool = True, make_diamond_db: bool = True,
        min_length: int = 30) -> Dict:
    run_dir = Path(run_dir)
    name = run_dir.name

    records = _uniprot_records(run_dir) + _kegg_records(run_dir)
    if not records:
        raise RuntimeError("no sequences found; run `sequences` and/or `ko-sequences` first")

    # Deduplicate on the sequence itself. The same protein reached here from both
    # directions; keep one copy and retain both provenances.
    by_hash: Dict[str, Dict] = {}
    for r in records:
        if len(r["seq"]) < min_length:
            continue
        h = hashlib.sha256(r["seq"].encode()).hexdigest()
        if h in by_hash:
            prev = by_hash[h]
            if r["source"] not in prev["source"]:
                prev["source"] = prev["source"] + "+" + r["source"]
                prev["also_acc"] = r["source_acc"]
            if prev["ko"] == "noKO" and r["ko"] != "noKO":
                prev["ko"], prev["all_kos"] = r["ko"], r["all_kos"]
        else:
            r["also_acc"] = ""
            by_hash[h] = r

    rows, fasta = [], []
    for i, (h, r) in enumerate(sorted(by_hash.items(), key=lambda x: (x[1]["ko"], x[1]["gene"])), 1):
        sid = f"bfgm|{i}|{_clean(r['ko'])}|{_clean(r['gene'])}|{_clean(r['source_acc'])}"
        fasta.append(f">{sid} {r['description']} [{r['organism']}]")
        fasta.extend(r["seq"][j:j + 60] for j in range(0, len(r["seq"]), 60))
        rows.append({"seq_id": sid, "n": i, "ko": r["ko"], "all_kos": r["all_kos"],
                     "gene": r["gene"], "source": r["source"],
                     "source_acc": r["source_acc"], "also_acc": r["also_acc"],
                     "organism": r["organism"], "taxid": r["taxid"],
                     "length": len(r["seq"]), "sha256": h[:16],
                     "description": r["description"]})

    fa = run_dir / f"{name}_reference.fasta"
    fa.write_text("\n".join(fasta) + "\n")
    mp = run_dir / f"{name}_reference_map.tsv"
    pd.DataFrame(rows).to_csv(mp, sep="\t", index=False)

    stats = {"sequences": len(rows), "distinct_kos": len({r["ko"] for r in rows if r["ko"] != "noKO"}),
             "from_uniprot_only": sum(1 for r in rows if r["source"] == "uniprot"),
             "from_kegg_only": sum(1 for r in rows if r["source"] == "kegg"),
             "in_both": sum(1 for r in rows if "+" in r["source"]),
             "duplicates_collapsed": len(records) - len(rows),
             "fasta": str(fa), "map": str(mp), "blast_db": None, "diamond_db": None}

    dbdir = run_dir / "blastdb"
    if make_blast_db and shutil.which("makeblastdb"):
        dbdir.mkdir(exist_ok=True)
        subprocess.run(["makeblastdb", "-in", str(fa), "-dbtype", "prot",
                        "-out", str(dbdir / name), "-title", f"bfgm {name}",
                        "-parse_seqids"], check=True, capture_output=True)
        stats["blast_db"] = str(dbdir / name)
    if make_diamond_db and shutil.which("diamond"):
        dbdir.mkdir(exist_ok=True)
        subprocess.run(["diamond", "makedb", "--in", str(fa),
                        "-d", str(dbdir / name)], check=True, capture_output=True)
        stats["diamond_db"] = str(dbdir / name) + ".dmnd"

    (run_dir / "HOW_TO_SEARCH.md").write_text(_usage(name, stats))
    return stats


def _usage(name: str, s: Dict) -> str:
    return f"""# Searching the {name} reference

`{name}_reference.fasta` holds {s['sequences']} deduplicated protein sequences covering
{s['distinct_kos']} KOs. `{name}_reference_map.tsv` maps every sequence ID to its KO,
gene, organism and description.

## Header format

    >bfgm|<n>|<KO>|<gene>|<source_acc> <description> [<organism>]

The ID has no whitespace, so `-outfmt 6` returns a `sseqid` that already contains the KO
and gene symbol. Split on `|` to get them without a join.

## Build the database

    makeblastdb -in {name}_reference.fasta -dbtype prot -out blastdb/{name} -parse_seqids
    diamond makedb --in {name}_reference.fasta -d blastdb/{name}

## Search your genomes against it

BLAST:

    blastp -query your_proteins.faa -db blastdb/{name} \\
           -outfmt '6 qseqid sseqid pident length evalue bitscore' \\
           -evalue 1e-10 -max_target_seqs 5 -num_threads 8 -out hits.tsv

DIAMOND, which is much faster for large query sets:

    diamond blastp -q your_proteins.faa -d blastdb/{name} \\
            --outfmt 6 qseqid sseqid pident length evalue bitscore \\
            --evalue 1e-10 --max-target-seqs 5 --threads 8 -o hits.tsv

## Turn hits into a KO profile

    python - <<'PY'
    import pandas as pd
    h = pd.read_csv("hits.tsv", sep="\\t", header=None,
                    names="qseqid sseqid pident length evalue bitscore".split())
    h["ko"] = h.sseqid.str.split("|").str[2]
    best = h.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
    best = best[(best.pident >= 40) & (best.length >= 80)]
    print(best.groupby("ko").size().sort_values(ascending=False))
    PY

## Thresholds

Defaults above are deliberately permissive. For orthology calls rather than
"something in this family is present", tighten to roughly 50 percent identity over
70 percent of the subject length, and confirm with reciprocal best hits. A single
permissive hit against a large transporter family means the fold is present, not that
the specific function is.

## What is NOT in here

Systems with no KEGG orthology entry and no reviewed UniProt record produce no
sequences and therefore cannot be found by searching this reference. Check
`ko_still_uncovered.csv` and the NOT_IN_KO rows of `gene_ko_map.csv` before concluding
that a genome lacks a function. Absence from this reference is not absence from biology.
"""
