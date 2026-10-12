#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""repeat_join: repeat-aware sequence joining with junction evidence summaries.

    detect     identify and trim terminal repeat copies
    junctions  count supporting alignment-segment pairs and alignment-start bins
    join       concatenate a user-specified path and rotate to a marker

The selected path and chromosome topology require separate biological assessment.
Start-position bins do not establish independent biological molecules.
Run each subcommand with --help. See the repository README for interpretation.
Requires minimap2, samtools, Biopython and pysam.
"""
import argparse, collections, pathlib, subprocess, shutil, sys

__version__ = "1.0"

# ── defaults ────────────────────────────────────────────────────────────────
# These are the values used in the study this tool comes from.  They are not magic;
# see README.md ("Choosing the criteria") before changing them.
END_WINDOW = 5000   # how close to a contig end an alignment must reach to count
MIN_ALN    = 1000   # minimum aligned length on each side of a junction
MIN_MAPQ   = 20     # minimum mapping quality on each side
ADJ_GAP    = 2000   # maximum unaligned gap between the two sides within one read
START_BIN  = 200    # bin width for alignment-start diversity, not independent molecules


def need(prog):
    if shutil.which(prog) is None:
        sys.exit(f"error: '{prog}' not found on PATH")


def run(cmd):
    print(f"  $ {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True)


# ── detect ──────────────────────────────────────────────────────────────────
def cmd_detect(a):
    """Align the repeat to the other contigs and trim the copies at their ends."""
    from Bio import SeqIO
    need("minimap2")
    out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    recs = {r.id: r for r in SeqIO.parse(a.assembly, "fasta")}

    repeat = a.repeat
    if repeat is None:
        if not a.info:
            sys.exit("error: give --repeat NAME or --info assembly_info.txt")
        flagged = []
        for line in open(a.info):
            f = line.rstrip("\n").split("\t")
            if len(f) > 4 and f[4] == "Y":
                flagged.append((int(f[1]), f[0]))
        if not flagged:
            sys.exit("error: no contig flagged as a repeat in --info; use --repeat")
        if len(flagged) > 1:
            print(f"note: {len(flagged)} contigs are flagged as repeats "
                  f"({', '.join(n for _, n in sorted(flagged, reverse=True))}); "
                  f"taking the longest. Use --repeat to choose a different one.")
        repeat = max(flagged)[1]
        print(f"repeat contig taken from {a.info}: {repeat}")
    if repeat not in recs:
        sys.exit(f"error: repeat contig '{repeat}' not in {a.assembly}")

    others = [k for k in recs if k != repeat]
    rfa = out / "_repeat.fasta"; SeqIO.write(recs[repeat], rfa, "fasta")
    ofa = out / "_others.fasta"; SeqIO.write([recs[k] for k in others], ofa, "fasta")
    paf = out / "repeat_vs_contigs.paf"
    # -N/-p keep every copy: the repeat aligns to several places by definition, and
    # --secondary=no would report only the best one and miss the other contig ends.
    run(f"minimap2 -x asm5 -N 50 -p 0.1 -t {a.threads} {ofa} {rfa} > {paf} 2>/dev/null")

    hits = collections.defaultdict(list)
    for line in open(paf):
        f = line.split("\t")
        if int(f[10]) < a.min_hit:
            continue
        hits[f[5]].append((int(f[7]), int(f[8]), int(f[6])))   # tstart, tend, tlen

    manual = {}
    for spec in (a.trim or []):
        try:
            n, s, e = spec.split(":")
            manual[n] = (int(s), int(e))
        except ValueError:
            sys.exit(f"error: --trim expects NAME:START:END, got '{spec}'")

    trim, rows = {}, []
    for name in others:
        L = len(recs[name].seq)
        if name in manual:
            s, e = manual[name]
            if not 0 <= s < e <= L:
                sys.exit(f"error: --trim {name}:{s}:{e} is outside 0..{L}")
        else:
            s, e = 0, L
            for ts, te, _ in sorted(hits.get(name, [])):
                if ts <= a.end_margin and te > s:
                    s = te
                if te >= L - a.end_margin and ts < e:
                    e = ts
        trim[name] = (s, e)
        rows.append((name, L, s, L - e, e - s))

    print(f"\n{'contig':<14}{'length':>12}{'trim 5′':>10}{'trim 3′':>10}{'kept':>12}")
    for n, L, s, t3, kept in rows:
        print(f"{n:<14}{L:>12,}{s:>10,}{t3:>10,}{kept:>12,}")
    print(f"{repeat + ' (repeat)':<14}{len(recs[repeat].seq):>12,}{'-':>10}{'-':>10}"
          f"{len(recs[repeat].seq):>12,}")

    if not any(s or L - e for n, L, s, e in
               [(n, L, t[0], t[1]) for n, L, t in
                [(n, len(recs[n].seq), trim[n]) for n in others]]):
        print("\nnote: no repeat copies were found at contig ends. Either the assembler "
              "already collapsed the repeat to one copy (in which case the assembly is "
              "already a usable target — pass it straight to `junctions`), or the repeat "
              "is internal and this tool does not apply.")

    segs = []
    for name in others:
        s, e = trim[name]
        r = recs[name][s:e]; r.id = r.name = name; r.description = ""
        segs.append(r)
    r = recs[repeat]; r.id = r.name = repeat; r.description = ""
    segs.append(r)
    tfa = out / "target_dedup.fasta"
    SeqIO.write(segs, tfa, "fasta")
    total = sum(len(x.seq) for x in segs)
    print(f"\nwrote {tfa}  ({len(segs)} sequences, {total:,} bp)")
    print(f"the repeat is now present once, as '{repeat}'")
    for f in (rfa, ofa):
        f.unlink()


# ── junctions ───────────────────────────────────────────────────────────────
def count_junctions(bam_path, lengths, p):
    import pysam
    bam = pysam.AlignmentFile(bam_path, "rb")
    reads = collections.defaultdict(list)
    for al in bam.fetch(until_eof=True):
        if al.is_unmapped or al.is_secondary or al.reference_length is None \
           or al.reference_length < p.min_aln:
            continue
        qs, qe = al.query_alignment_start, al.query_alignment_end
        if al.is_reverse:                      # put coordinates back on the read
            ql = al.infer_read_length() or (qe - qs)
            qs, qe = ql - qe, ql - qs
        reads[al.query_name].append(dict(t=al.reference_name, ts=al.reference_start,
                                         te=al.reference_end, mq=al.mapping_quality,
                                         qs=qs, qe=qe))

    def which_ends(seg):
        out = []
        if seg["ts"] <= p.end_window:
            out.append("5")
        if seg["te"] >= lengths[seg["t"]] - p.end_window:
            out.append("3")
        return out

    J = collections.defaultdict(list)
    for segs in reads.values():
        segs.sort(key=lambda x: x["qs"])
        for a, b in zip(segs, segs[1:]):
            if b["qs"] - a["qe"] > p.adj_gap or min(a["mq"], b["mq"]) < p.min_mapq:
                continue
            ea, eb = which_ends(a), which_ends(b)
            if not ea or not eb:
                continue
            key = tuple(sorted([f"{a['t']}:{ea[-1]}", f"{b['t']}:{eb[-1]}"]))
            J[key].append(a["ts"])
    return J


def cmd_junctions(a):
    from Bio import SeqIO
    need("minimap2"); need("samtools")
    out = pathlib.Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    lengths = {r.id: len(r.seq) for r in SeqIO.parse(a.target, "fasta")}

    bam = out / "reads_vs_target.bam"
    if bam.exists() and not a.force:
        print(f"using existing {bam} (--force to redo)")
    else:
        run(f"minimap2 -ax map-ont -t {a.threads} --secondary=no {a.target} {a.reads} "
            f"| samtools sort -@{max(1, a.threads // 2)} -o {bam} - "
            f"&& samtools index {bam}")

    J = count_junctions(bam, lengths, a)
    print(f"\njunction evidence  (end window {a.end_window}, min aligned {a.min_aln}, "
          f"MAPQ >= {a.min_mapq}, gap <= {a.adj_gap}, start bin {a.start_bin})\n")
    print(f"{'junction':<34}{'reads':>8}{'molecules':>12}   note")
    rows = []
    for k, starts in sorted(J.items(), key=lambda x: -len(x[1])):
        indep = len({s // a.start_bin for s in starts})
        # A "junction" whose two sides are the same segment end is not a join between
        # segments; it comes from chimeric reads and from reads folding back on the
        # repeat. It is noise, and is reported only so that it is not mistaken for
        # missing output.
        self_j = k[0].split(":")[0] == k[1].split(":")[0]
        rows.append((k, len(starts), indep, self_j))
        note = "self — ignore" if self_j else ""
        print(f"{k[0] + '  <->  ' + k[1]:<34}{len(starts):>8}{indep:>12}   {note}")
    with (out / "junctions.tsv").open("w") as fh:
        fh.write("side_a\tside_b\treads\tindependent_molecules\tself_junction\n")
        for k, n, i, s in rows:
            fh.write(f"{k[0]}\t{k[1]}\t{n}\t{i}\t{'yes' if s else 'no'}\n")
    print(f"\nwrote {out / 'junctions.tsv'}")

    real = [r for r in rows if not r[3]]
    weak = [r for r in real if r[2] < a.warn_below]
    if weak:
        print(f"\nwarning: {len(weak)} of {len(real)} segment-to-segment junction(s) "
              f"have fewer than {a.warn_below} independent molecules:")
        for k, n, i, _ in weak:
            print(f"    {k[0]}  <->  {k[1]}    {n} reads but only {i} molecule(s)")
        print("  Report the independent-molecule count, not the read count: reads "
              "sharing a start\n  position are one molecule seen repeatedly, not "
              "independent evidence. A junction\n  resting on one molecule is the "
              "weakest point of the assembly and should be stated\n  as such.")


# ── join ────────────────────────────────────────────────────────────────────
def cmd_join(a):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    recs = {r.id: r for r in SeqIO.parse(a.target, "fasta")}

    seq, table = None, []
    for token in a.path.split(","):
        token = token.strip()
        rc = token.startswith("-")
        name = token[1:] if rc else token
        if name not in recs:
            sys.exit(f"error: '{name}' is not in {a.target}. "
                     f"Available: {', '.join(recs)}")
        piece = recs[name].seq.reverse_complement() if rc else recs[name].seq
        seq = piece if seq is None else seq + piece
        table.append((token, len(piece), len(seq)))

    print(f"{'segment':<16}{'length':>14}{'cumulative':>16}")
    for t, l, c in table:
        print(f"{t:<16}{l:>14,}{c:>16,}")
    print(f"\njoined length: {len(seq):,} bp")

    if a.rotate_to:
        marker = str(SeqIO.read(a.rotate_to, "fasta").seq).upper()[:a.marker_bp]
        s = str(seq).upper()
        i = (s + s).find(marker)
        if i < 0:
            i = (s + s).find(str(Seq(marker).reverse_complement()))
            if i >= 0:
                seq = seq.reverse_complement()
                s = str(seq).upper()
                i = (s + s).find(marker)
        if 0 <= i < len(seq):
            seq = seq[i:] + seq[:i]
            print(f"rotated: marker now at position 1 (was at {i + 1:,})")
        elif i >= len(seq):
            sys.exit("error: marker found only in the doubled sequence, which should "
                     "not happen; check the marker file")
        else:
            print("warning: marker not found, sequence NOT rotated")

    rec = SeqRecord(seq, id=a.name, description="circular=true")
    SeqIO.write(rec, a.out, "fasta")
    print(f"wrote {a.out}  ({len(seq):,} bp)")


# ── cli ─────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"repeat_join {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    d = sub.add_parser("detect", help="trim repeat copies left at contig ends")
    d.add_argument("--assembly", required=True, help="assembly FASTA")
    d.add_argument("--info", help="Flye assembly_info.txt, to identify the repeat")
    d.add_argument("--repeat", help="name of the repeat contig (overrides --info)")
    d.add_argument("--outdir", required=True)
    d.add_argument("--end-margin", type=int, default=END_WINDOW,
                   help=f"how close to a contig end a copy must lie (default {END_WINDOW})")
    d.add_argument("--min-hit", type=int, default=MIN_ALN,
                   help=f"minimum alignment length to call a copy (default {MIN_ALN})")
    d.add_argument("--trim", action="append", metavar="NAME:START:END",
                   help="override the detected coordinates for one contig; repeatable. "
                        "Use this to reproduce a published assembly exactly, since "
                        "different aligners place the ends of a repeat copy a few tens "
                        "of bases apart")
    d.add_argument("--threads", type=int, default=8)
    d.set_defaults(func=cmd_detect)

    j = sub.add_parser("junctions", help="count reads spanning repeat boundaries")
    j.add_argument("--target", required=True, help="FASTA with the repeat present once")
    j.add_argument("--reads", required=True, help="reads FASTQ (the assembly input)")
    j.add_argument("--outdir", required=True)
    j.add_argument("--end-window", type=int, default=END_WINDOW)
    j.add_argument("--min-aln", type=int, default=MIN_ALN)
    j.add_argument("--min-mapq", type=int, default=MIN_MAPQ)
    j.add_argument("--adj-gap", type=int, default=ADJ_GAP)
    j.add_argument("--start-bin", type=int, default=START_BIN)
    j.add_argument("--warn-below", type=int, default=3,
                   help="warn about junctions with fewer independent molecules")
    j.add_argument("--threads", type=int, default=8)
    j.add_argument("--force", action="store_true", help="redo the mapping")
    j.set_defaults(func=cmd_junctions)

    o = sub.add_parser("join", help="concatenate along a path and rotate")
    o.add_argument("--target", required=True)
    o.add_argument("--path", required=True,
                   help="comma-separated segment names, '-' to reverse complement, "
                        "e.g. 'contig_1,contig_3,-contig_2,contig_3'")
    o.add_argument("--out", required=True)
    o.add_argument("--name", default="chromosome", help="FASTA header id")
    o.add_argument("--rotate-to", help="FASTA of a marker gene, e.g. the dnaA CDS")
    o.add_argument("--marker-bp", type=int, default=600,
                   help="how much of the marker to match (default 600)")
    o.set_defaults(func=cmd_join)

    a = p.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
