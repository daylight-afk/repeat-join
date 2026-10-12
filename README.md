# repeat_join

A custom command-line tool for repeat-aware sequence joining, junction alignment evidence summaries, and marker-based rotation.

This repository contains only the tool, installation instructions and synthetic tests. It contains no experimental reads, study result tables, manuscript files or project working history.

## Install

Use Python 3 in an isolated environment:

```sh
python3 -m venv .venv
. .venv/bin/activate
python3 -m pip install -r repeat_join/requirements.txt
python3 repeat_join/repeat_join.py --help
```

The `detect` command additionally requires minimap2 on PATH. `junctions` requires minimap2 and samtools. `join` uses Biopython. The requirements file gives minimum Python package versions, not a full environment lockfile.

## Commands

- `detect`: align the repeat against other contigs and trim terminal repeat copies; manual trim coordinates can override the inferred trim.
- `junctions`: map reads to the resulting target and summarize adjacent alignment-segment evidence at contig ends.
- `join`: concatenate the user-specified oriented segment path and optionally rotate to a marker sequence.

```sh
python3 repeat_join/repeat_join.py detect --help
python3 repeat_join/repeat_join.py junctions --help
python3 repeat_join/repeat_join.py join --help
```

Synthetic example (no experimental data or external aligners required):

```sh
python3 repeat_join/repeat_join.py join \
  --target tests/fixtures/segments.fasta \
  --path 'U1,R,-U2,R' \
  --name synthetic \
  --out synthetic_joined.fasta
```

A minus sign reverses and complements a segment. Repeated names include the segment multiple times. `--trim NAME:START:END` coordinates use zero-based, end-exclusive Python slicing. Real data require a justified path and trim coordinates; the example is synthetic.

## Interpretation and limitations

Defaults: 5,000 bp contig-end window; minimum alignment reference span 1,000 bp; MAPQ >=20 on both segments; query gap <=2,000 bp; 200 bp start-position bins. Secondary and unmapped alignments are excluded; supplementary alignments are eligible. The span threshold uses `reference_length`, not just M/=/X bases.

The legacy output column `independent_molecules` counts alignment-start bins; it does not establish biological independence. The `reads` column counts adjacent alignment-segment pairs, not necessarily unique reads. Column names remain unchanged for compatibility.

The legacy `self_junction` flag marks pairs on the same contig, even when ends differ. It does not validate same-contig circular closure. `join` adds a `circular=true` label but does not prove circularity or chromosome completeness. A missing rotation marker produces a warning; inspect the log and output. Existing mapping files may be reused unless `--force` is supplied.

The original command wrapper constructs shell commands, so use trusted paths without whitespace or shell metacharacters. Run in a dedicated working directory. This tool supports a user-selected reconstruction; it does not automatically establish a correct biological assembly.

## Tests

```sh
python3 -m unittest discover -s tests -v
```

Six synthetic tests cover joining/reverse complementation, marker rotation, unknown-segment rejection, supplementary alignment pairs/start bins, and alignment filtering. They do not validate complete experimental genome assemblies. Packaging checks passed with Biopython 1.84 and pysam 0.22.1.

## License

No explicit software reuse license has been assigned. Public source availability alone does not grant an additional reuse license.
