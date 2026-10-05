# TBM vs AlphaFold3: model quality comparison

- TBM model: `../examples/test_target/models/first_pass/test_target.B99990003.pdb` (MODELLER, template 4EHX_A, residues 103-460)
- AF3 model: `../examples/test_target/af3/output/test_target/test_target_model.cif` (local AlphaFold3 3.0.0, full length 1-460, top-ranked of 5 samples)

## Structural agreement, residue-matched (TMscore, residues 103-460)

Same sequence, so each TBM residue is compared only with the same AF3 residue.

| Metric | Value |
|---|---|
| TM-score | 0.554 |
| RMSD (optimal superposition, 358 residues) | 8.85 Å |
| GDT-TS | 0.325 |

## Structural agreement, sequence-independent (TM-align)

| Metric | Value |
|---|---|
| TM-score (normalized by TBM, 358 res) | 0.643 |
| TM-score (normalized by AF3, 460 res) | 0.520 |
| RMSD over aligned residues | 4.92 Å |
| Aligned length | 307 |
| Aligned pairs within 5 Å | 221 |
| Aligned pairs at the same residue number | 137 of 307 |
| Residue-matched CA RMSD, all 358 residues (same superposition) | 9.20 Å |

## Model quality, residues 103-460

| Metric | TBM | AlphaFold3 |
|---|---|---|
| DOPE (lower = better) | -29215.4 | -44235.6 |
| z-DOPE (lower = better; < -1 native-like) | 1.160 | -1.811 |
| TM-score vs template 4EHX_A | 0.742 | 0.629 |
| Mean pLDDT (all / 103-460) | n/a | 95.6 / 95.9 |
| pTM / ranking score | n/a | 0.93 / 0.93 |

## Regions

| Region | AF3 pLDDT | TBM DOPE | AF3 DOPE | CA deviation (Å) |
|---|---|---|---|---|
| 1-102 | 95.1 | not modeled | not scored | n/a |
| 182-215 | 90.1 | -0.023 | -0.0396 | 9.18 |
| 224-230 | 96.2 | -0.0117 | -0.0239 | 7.98 |
| 331-345 | 96.8 | -0.03 | -0.0384 | 10.61 |
