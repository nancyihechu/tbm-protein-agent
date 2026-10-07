# Agent 2 experimental-reference evaluation

Status: **PARTIAL_ASSIGNMENT**

Owner: Matthew. Orchestration: ChatGPT (Astra), with deterministic scoring tools.

## Input and assignment status

| Target | Reference | CASP domain | Classification | AF3 server |
|---|---|---|---|---|
| T1147 | 8EM5:A | 12-103 | TBM-easy | NOT_PROVIDED |
| T1133 | 8DYS:A | 4-427 | TBM-easy | NOT_PROVIDED |
| T1183 | 8IFX:B | 1-195 | TBM-easy | NOT_PROVIDED |

## Accuracy against experiment

| Target | Method | Status | Paired CA/domain length | TM-score | GDT-TS (0-1) | lDDT-CA | RMSD (A) |
|---|---|---|---:|---:|---:|---:|---:|
| T1147 | tbm | SCORED | 92/92 | 0.6338 | 0.6005 | 0.5290 | 7.8043 |
| T1147 | alphafold3_local | SCORED | 92/92 | 0.9827 | 0.9946 | 0.9764 | 0.4726 |
| T1147 | alphafold3_server | NOT_PROVIDED | pending | pending | pending | pending | pending |
| T1133 | tbm | SCORED | 399/424 | 0.8535 | 0.6946 | 0.7295 | 3.8127 |
| T1133 | alphafold3_local | SCORED | 399/424 | 0.9266 | 0.9021 | 0.9433 | 1.0119 |
| T1133 | alphafold3_server | NOT_PROVIDED | pending | pending | pending | pending | pending |
| T1183 | tbm | SCORED | 195/195 | 0.7495 | 0.6090 | 0.6445 | 5.2934 |
| T1183 | alphafold3_local | SCORED | 195/195 | 0.9866 | 0.9808 | 0.9626 | 0.6179 |
| T1183 | alphafold3_server | NOT_PROVIDED | pending | pending | pending | pending | pending |

## Interpretation and limits

- TBM is Agent 1's DOPE-selected model; the local AF3 model is Agent 1's top-ranked output. Selection precedes experimental scoring.
- Each target uses one fixed intersection of observed reference and supplied model C-alpha residues within its CASP domain. No score-based residue pruning.
- TM-score and GDT-TS use the declared full CASP domain length. Missing positions remain in that denominator. These are our reproducible calculations, not official CASP assessor scores.
- RMSD uses a Kabsch fit over the fixed paired C-alpha mask. lDDT-CA uses that mask and is not all-atom lDDT or pLDDT. See METRICS.md for definitions.
- T1133's TBM covers 414/585 target residues; its experimental reference and CASP domain also cover only part of the full sequence. High regional agreement does not establish full-length accuracy.
- Local AF3 models were generated on Hellbender. Only rows labeled alphafold3_server satisfy the web-server comparison requirement; inspect their input status above.
- All methods use the same mask per run. Adding a new method can change the intersection; compare methods within one run.
- Existing Agent 1 scores compare TBM with AF3; those remain separate from this experimental-reference evaluation.

## Remaining deliverables

- T1147: obtain AlphaFold web-server output using the identical FASTA; record job provenance and evaluate it.
- T1133: obtain AlphaFold web-server output using the identical FASTA; record job provenance and evaluate it.
- T1183: obtain AlphaFold web-server output using the identical FASTA; record job provenance and evaluate it.

## Evidence

- [CASP15 target assignments](https://predictioncenter.org/casp15/targetlist.cgi)
- [CASP15 domain definitions and classification](https://predictioncenter.org/casp15/domains_summary.cgi)
- Source checksums, software versions and configuration are saved with this report.
