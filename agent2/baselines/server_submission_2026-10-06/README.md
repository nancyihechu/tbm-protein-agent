# AlphaFold Server jobs completed October 6, 2026

All three web-server jobs completed and were downloaded through the signed-in
Chrome AlphaFold Server session. Agent 2 validated the downloads, fixed primary
model selection, then evaluated those models against experiment. See the
[complete nine-comparison report](../../results/with_server_2026-10-06/REPORT.md)
and [portable evaluation configuration](../../config_with_server.json).

| Target | Actual server job | Sequence length | Seed | Primary file index | Exported ranking score |
|---|---|---:|---:|---:|---:|
| T1147 | [41cb32c284a7a5d7](https://alphafoldserver.com/fold/41cb32c284a7a5d7) | 103 | 1 | 0 | 0.83 |
| T1133 | [6d32ec7dd204db3e](https://alphafoldserver.com/fold/6d32ec7dd204db3e) | 585 | 1 | 0 | 0.87 |
| T1183 | [156d096f94d69fb](https://alphafoldserver.com/fold/156d096f94d69fb) | 200 | 1 | 0 | 0.92 |

Ranking scores are server confidence/ranking values, not experimental accuracy.
Each archive contains five predictions. Selection used the highest exported
`ranking_score` before reference scoring, with the lowest file index breaking
rounded ties. T1147 files 0 and 1 both export 0.83; file 0 was selected. The
importer does not read experimental coordinates.

## Inputs and provenance

The jobs were submitted by pasting each exact sequence into the browser form.
Each exported request confirms the expected job name, one unmodified protein
copy, seed 1, enabled PDB templates and no user-supplied MSA or custom template.
The request sequence and selected model's full sequence match Nancy's FASTA
exactly. Actual exported requests use schema version 3.

`jobs.json`, individual job JSONs and `manifest.json` preserve preparation before
submission; their `PREPARED_NOT_SUBMITTED` status is historical. They were not
uploaded successfully, and their version-1 settings are not proof of actual
server inputs. The authoritative records are the downloaded `*_job_request.json`
files and [import provenance](imported/import_provenance.json). The separate
[completion record](completion.json) records the later evaluated state.

Template-date evidence is deliberately distinguished:

- T1147 explicitly exports `maxTemplateDate: "2021-09-30"`.
- T1133 and T1183 used the server's default template setting and omit
  `maxTemplateDate` from their exported requests. Their resolved cutoff is not
  independently verified by those files.
- The submission UI displayed the default cutoff as `29/09/2021`. For the first
  job, entering ISO date `2021-09-30` also displayed `29/09/2021`, while the
  download retained `2021-09-30`. This observed display discrepancy is consistent
  with date/timezone formatting, but its cause was not established. It does not
  justify inserting an inferred date into the other exported requests.
- These October 2026 runs are retrospective. Template restrictions do not prove
  absence of target-related information in model training; the local and server
  runs are not a controlled comparison of deployment environments.

`imported/<target>/` preserves the original ZIP, submitted FASTA, every extracted
file (models, confidence/PAE data, MSAs, templates and terms), and per-target
provenance with SHA-256 checksums. [Completion screenshot](server_completed_all.jpg)
shows all three completed jobs. The original six-comparison evaluation and
Nancy's files remain unchanged.

## Reproduce the import and evaluation

From the repository root, use `import_server_outputs.py --help`. Provide exactly
three `--archive TARGET=ZIP` and `--job-url TARGET=URL` arguments, or an
`--inputs-json` file. The importer validates all archives before creating a fresh
output folder and refuses an existing destination. It writes a copied evaluation
configuration in that new folder; it does not modify the base config.

For the preserved import, use:

```bash
python -m agent2 evaluate --config agent2/config_with_server.json --tmscore /path/to/TMscore --require-complete
```

The published evaluation has status `COMPLETE`: nine scored rows, no missing
metrics or pending server inputs, and common masks of 92, 399 and 195 C-alpha
positions for T1147, T1133 and T1183, respectively. See the report for domain
normalization and the T1133 partial-coverage limitation.

## Output notice

AlphaFold Server outputs and derivatives are provided under and subject to the
[AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms).
Each extracted output folder and original archive contains `terms_of_use.md`.
Retain it with the output. Original predictions are unmodified. Agent 2's derived
files extract/renumber paired C-alpha atoms, apply rigid superpositions, and
calculate metrics and figures.

Citation: Abramson et al., [Accurate structure prediction of biomolecular
interactions with AlphaFold 3](https://doi.org/10.1038/s41586-024-07487-w),
Nature (2024). The [official server JSON documentation](https://github.com/google-deepmind/alphafold/blob/main/server/README.md)
describes the request fields; actual downloaded files remain authoritative for
these jobs.
