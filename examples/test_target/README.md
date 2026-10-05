# Demo run: test_target (CASP15 T1112)

The first end-to-end run of the pipeline. `test_target.fasta` has the same
sequence as T1112 in `targets/candidates_backup.fasta` (460 residues,
O93732.1, *Methanothermus fervidus*).

| Step | Output |
|---|---|
| Template search (`run_template_search.slurm`) | `hhsearch_output/test_target.a3m`, `.hhr` |
| Template chains | `templates/4ehx_A.pdb`, `templates/4itl_A.pdb` |
| Alignment | `alignments/test_target_4ehx_A.pir` |
| MODELLER, residues 103-460 | `models/first_pass/` (5 models, `run_modeller.log`) |
| DOPE profile figures | `models/first_pass/figures/` |
| AlphaFold3 (`run_af3.slurm`, module `alphafold3/alphafold3_deepmind`, package 3.0.0) | `af3/output/test_target/` |
| TBM vs AF3 comparison and presentation figures | `comparison_results/` |

Key numbers (from `comparison_results/presentation/summary.json`):
template 4EHX_A (14% identity, query 103-448), best model
`test_target.B99990003.pdb` (DOPE -29215.4, z-DOPE 1.16), AF3 mean pLDDT
95.62 and pTM 0.93, TM-score vs AF3 0.554 residue-matched and 0.643 with
TM-align (RMSD 4.92 A over 307 aligned residues).

This demo used an earlier version of `build_pir.py` that kept the template
overhangs in the alignment; current runs keep only the aligned range. The
SLURM scripts here are the demo versions; `SCRIPTS/` has the per-target ones.
