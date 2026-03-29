import re
from collections import defaultdict
from pathlib import Path

import pandas as pd


def norm(s: str) -> str:
    s = (s or "").strip()
    s = re.sub(r"\s+", " ", s)
    return s.lower()


def norm_strain(s: str) -> str:
    s = (s or "").strip()
    # Example: .../2023(H5N1)) or .../2023(H5N1) -> .../2023
    s = re.sub(r"\([HhNn0-9]+\)?$", "", s)
    s = re.sub(r"\s+", " ", s)
    return s.lower()


def main():
    repo = Path(__file__).resolve().parents[1]
    meta_path = repo / "config" / "final_metadata_50_per_country_isolates.tsv"
    fasta_path = repo / "final_sequences_from_final_isolates.fasta"

    meta = pd.read_csv(meta_path, sep="\t", dtype=str).fillna("")

    seg_by_acc = {}
    for _, r in meta.iterrows():
        acc = r.get("accession", "").strip()
        seg = r.get("segment", "").strip()
        if acc and seg.isdigit():
            seg_by_acc[acc] = int(seg)

    seg_by_fasta = {}

    acc_re = re.compile(r"\b([A-Z]{1,3}\d+\.\d+)\b")
    acc_first_re = re.compile(r"^\s*([A-Z]{1,3}\d+\.\d+)\b")
    strain_re = re.compile(r"\((A/.+)\)\s+segment\s+\d", re.IGNORECASE)

    by_strain = defaultdict(set)
    by_isolate = defaultdict(set)

    with fasta_path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            h = line[1:].strip()

            m_acc_first = acc_first_re.search(h)
            m_strain = strain_re.search(h)
            if m_acc_first and m_strain:
                acc = m_acc_first.group(1)
                strain = m_strain.group(1)
                m_seg = re.search(r"segment\s+(\d)", h, re.IGNORECASE)
                if m_seg:
                    seg_by_fasta[acc] = int(m_seg.group(1))
                by_strain[norm_strain(strain)].add(acc)
                parts = [p for p in strain.split("/") if p]
                if len(parts) >= 2:
                    by_isolate[norm(parts[-2])].add(acc)
                continue

            # Headers like: PA_CH-PD032_rectal_PQ002153.1/PA/Argentina/2023
            first = h.split("/")[0]
            m_acc = acc_re.search(first)
            if not m_acc:
                continue

            acc = m_acc.group(1)
            m_seg = re.search(r"segment\s+(\d)", h, re.IGNORECASE)
            if m_seg:
                seg_by_fasta[acc] = int(m_seg.group(1))
            prefix = first[: m_acc.start()].rstrip("_")
            if not prefix:
                continue

            by_isolate[norm(prefix)].add(acc)

            # Base key to collapse tissue/suffix variants.
            p2 = re.sub(r"^(pb2|pb1|pa|ha|np|na|ns|mp|m|nep)_", "", prefix, flags=re.IGNORECASE)
            p2 = re.sub(r"_(oral|lung|brain|tracheal|rectal|swab)$", "", p2, flags=re.IGNORECASE)
            if p2 and p2.lower() != prefix.lower():
                by_isolate[norm(p2)].add(acc)

    all_accessions = []
    segment_counts = []
    for _, r in meta.iterrows():
        acc = r.get("accession", "").strip()
        strain = r.get("strain", "").strip()
        isolate = r.get("isolate", "").strip()

        cands = set()
        if strain:
            cands |= by_strain.get(norm_strain(strain), set())
        if isolate:
            cands |= by_isolate.get(norm(isolate), set())
            iso_base = re.sub(r"^(pb2|pb1|pa|ha|np|na|ns|mp|m|nep)_", "", isolate, flags=re.IGNORECASE)
            iso_base = re.sub(r"_(oral|lung|brain|tracheal|rectal|swab)$", "", iso_base, flags=re.IGNORECASE)
            cands |= by_isolate.get(norm(iso_base), set())

        cands = {c for c in cands if acc_re.fullmatch(c)}
        if acc:
            cands.add(acc)

        # Keep only one accession per segment, preferring current-row accession.
        chosen_by_segment = {}
        for cand in sorted(cands):
            seg = seg_by_acc.get(cand, seg_by_fasta.get(cand))
            if seg is None:
                continue
            prev = chosen_by_segment.get(seg)
            if prev is None:
                chosen_by_segment[seg] = cand
                continue
            if cand == acc and prev != acc:
                chosen_by_segment[seg] = cand

        ordered = [chosen_by_segment[s] for s in sorted(chosen_by_segment)]
        all_accessions.append(";".join(ordered))
        segment_counts.append(str(len(chosen_by_segment)) if chosen_by_segment else "")

    meta["all_accessions"] = all_accessions
    meta["segment_count"] = segment_counts
    cols = [c for c in meta.columns if c != "all_accessions"] + ["all_accessions"]
    meta = meta[cols]
    meta.to_csv(meta_path, sep="\t", index=False)

    print("Rows:", len(meta))
    print("Rows with >1 accession:", sum(";" in x for x in all_accessions))

    for target in ["PV290336.1", "PV290335.1", "PQ002143.1", "PQ002152.1"]:
        row = meta.loc[meta["accession"].eq(target), ["accession", "all_accessions"]]
        if not row.empty:
            print(row.to_string(index=False))


if __name__ == "__main__":
    main()
