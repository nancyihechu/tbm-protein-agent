# Agent 2 scoring conventions

The study uses the official CASP15 domain intervals but scores the downloaded RCSB
coordinates ourselves. It is a retrospective course comparison, not a recreation
of the blind CASP15 assessment or a claim about training-set independence.

## Correspondence and coverage

Biopython reads one explicit author chain and coordinate model (zero-based model
index). Target sequences are checked against every mapped amino acid. The supplied
manifest's author-number mapping was verified for both predictions and references.
Sequence alignment or explicit identity-checked maps are available for new inputs.
Insertion codes, missing C-alpha atoms, and ambiguity must be handled explicitly.

Before any fit, the workflow freezes the intersection of observed reference and
all supplied model C-alpha positions within the specified CASP domain. All methods
for a target share that mask, without score-based pruning. Every target position
and exclusion reason appears in `residue_mapping.csv`. Adding a server baseline
can change the intersection; interpret paired scores from the same run.

| Target | CASP domain | TM/GDT normalization length | Primary paired CA count |
|---|---|---:|---:|
| T1147 | 12-103 | 92 | 92 |
| T1133 | 4-427 | 424 | 399 |
| T1183 | 1-195 | 195 | 195 |

T1133 has 15 unresolved positions inside the reference domain and 10 additional
observed positions outside the TBM's 5-418 modeled interval. Scores normalized by
424 therefore retain those missing positions in their denominator. This is an
explicit conservative normalization choice; it is not the same as normalizing by
399 matched residues, by 409 experimentally observed residues, or by 585 full
target residues. Full-target model coverage is reported separately.

## Metric definitions

| Field | Calculation | Units / direction |
|---|---|---|
| `rmsd_ca_angstrom` | Least-squares Kabsch rigid fit over all fixed paired C-alpha positions; proper rotation only, no reflection or outlier rejection | Angstrom; lower is better |
| `tm_score` | Official TMscore optimization using sequence-corresponding C-alpha PDBs and `-l` equal to the declared CASP domain length; parse the user-normalized output | 0-1; higher is better |
| `gdt_ts` | Official TMscore optimized GDT-TS, using the 1, 2, 4, 8 angstrom thresholds; convert its observed-input normalization to the same declared domain length | 0-1, not percent; higher is better |
| `lddt_ca_common_mask` | Fraction of reference C-alpha distances preserved, averaged over strict absolute-error thresholds of 0.5, 1, 2, 4 angstrom; distinct-residue reference pairs closer than 15 angstrom, within the frozen mask; no fit | 0-1; higher is better |

The lDDT value is **C-alpha only on the shared observed mask**. It is not all-atom
lDDT, pLDDT, a stereochemical check, or a whole-sequence score. Unobserved positions
outside the mask do not contribute contacts; evaluate it together with coverage.
Zero eligible contacts produces null, never a fabricated perfect score.

TMscore receives derived C-alpha PDBs renumbered to target positions, preserving
the shared correspondence. PDB coordinates are rounded to 0.001 angstrom. The
wrapper validates both input lengths and the matched residue count. Its GDT-TS
output uses native input length even with `-l`, so the wrapper rescales by
`matched_count / domain_length`. The source, command, executable hash, stdout,
stderr and any errors are retained. A missing tool or unrecognized output yields
null TM/GDT, not a Kabsch-derived approximation labeled as optimized scores.

The overlay is the **Kabsch** fit used for RMSD. It does not claim to depict the
separate TM-score or GDT optimizing fit. It displays only evaluated C-alpha traces;
curves break across missing positions. Original all-atom models remain available
unchanged in the Agent 1 folders.

## Interpretation

DOPE and GA341 are model-quality estimates; AF3 pLDDT and pTM are confidence
estimates. They are not experimental accuracy and do not replace these comparisons.
Nancy's existing TBM-versus-AF3 scores remain a separate agreement analysis.

Primary predictions are selected before reference scoring. The three proteins
give a descriptive comparison. Partial models, unresolved reference regions,
different database availability and retrospective target exposure limit broader
claims. Template accessions are checked against the native accession; this check
does not prove absence of homologous or learned information.

The October 6 web-server runs are retrospective predictions. An explicit template
cutoff constrains template search; it does not establish training independence.
Server and local runs can also differ in model release, databases, sampling and
other defaults. Their score differences cannot isolate an effect of deployment.
When a downloaded request omits a template date, record the server default and
the observed UI evidence separately instead of asserting an exported date.

## Primary sources

- [CASP15 domain definitions](https://predictioncenter.org/casp15/domains_summary.cgi)
- [Official TMscore source and output](https://github.com/pylelab/USalign/blob/master/TMscore.cpp)
- [Official score normalization implementation](https://github.com/pylelab/USalign/blob/master/TMscore.h)
- [OpenStructure lDDT definition](https://openstructure.org/docs/2.12/mol/alg/lddt/)
- [Mariani et al., lDDT paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3799472/)

Dependency versions, code hashes, reference hashes, and the configuration snapshot
are recorded per run. Scientific computations and missing/invalid-input behavior
are tested under `agent2/tests/`.
