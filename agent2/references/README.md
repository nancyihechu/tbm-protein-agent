# Experimental reference provenance

These three unmodified mmCIF files were downloaded from RCSB PDB on **2026-10-06** for Agent 2, **ChatGPT (Astra)**. The associated model source commit is `d65f2b5ae2acb150704bf20347604c3de3b896c3`. Paths in `agent2/config.json` are relative to the repository root.

Reference identities and domain intervals were selected from the [official CASP15 target list](https://predictioncenter.org/casp15/targetlist.cgi) and [official CASP15 domain definition summary](https://predictioncenter.org/casp15/domains_summary.cgi). Chain selection used target sequence identity and a deterministic rule, independently of any prediction-versus-experiment accuracy score: the first author-chain copy for 8EM5 (A), the sole protein chain for 8DYS (A), and the TsaB chain for 8IFX (B). Model index `0` means the first coordinate model, PDB/mmCIF model number `1`.

| Target | Reference and selected author chain | X-ray resolution | Official evaluation unit | Interval and normalization length | CASP class |
|---|---|---:|---|---|---|
| T1147 | [8EM5](https://www.rcsb.org/structure/8EM5), A | 1.95 A | T1147-D1 | 12-103; 92 residues | TBM-easy |
| T1133 | [8DYS](https://www.rcsb.org/structure/8DYS), A | 1.80 A | T1133-D1 | 4-427; 424 residues | TBM-easy |
| T1183 | [8IFX](https://www.rcsb.org/structure/8IFX), B | 2.00 A | T1183-D1 | 1-195; 195 residues | TBM-easy |

The CASP classification applies to the official evaluation units. It is separate from template identity or screening labels assigned by this repository. The normalization lengths above are the inclusive official domain lengths, not the number of observed coordinates or the size of a prediction/reference intersection.

## Sequence and residue-number verification

Each deposited polymer sequence was compared with the repository target FASTA through the RCSB polymer-entity API. The downloaded files were then checked directly: every selected-chain `ATOM` C-alpha residue in model 1 matches the FASTA residue indexed by its `auth_seq_id`. No sequence mismatch was found. The configured mapping is therefore `author_number_equals_target`; evaluation must use author residue numbers rather than assuming mmCIF label numbering.

- **T1147:** [CASP target page](https://predictioncenter.org/casp15/target.cgi?id=102&view=all), [RCSB entity 1](https://data.rcsb.org/rest/v1/core/polymer_entity/8EM5/1), [chain A mapping](https://data.rcsb.org/rest/v1/core/polymer_entity_instance/8EM5/A). The deposited 103-residue sequence exactly equals `targets/T1147/T1147.fasta`; it includes the initial `GPLGS`. Identical protein copies are author chains A-F. Chain A contains 92 observed C-alpha residues, positions **12-103**; positions **1-11** have no coordinates. Its author and label sequence numbering both equal target numbering.
- **T1133:** [CASP target page](https://predictioncenter.org/casp15/target.cgi?id=74&view=all), [RCSB entity 1](https://data.rcsb.org/rest/v1/core/polymer_entity/8DYS/1), [chain A mapping](https://data.rcsb.org/rest/v1/core/polymer_entity_instance/8DYS/A). The 603-residue deposited sequence contains the 18-residue tag `MHHHHHHSSGRENLYFQG` followed by the exact 585-residue `targets/T1133/T1133.fasta`. Tag author numbers are **-17 through 0**. For target residues, **label_seq_id = target position + 18**, while **auth_seq_id = target position**. There are 409 observed C-alpha residues at target positions **4-98, 105-241, 251-427**. Target positions **1-3, 99-104, 242-250, 428-585** have no coordinates. The official domain spans 424 positions and includes 15 internal positions without coordinates. The TBM model covers only **5-418**, giving **399** observed reference residues in the common range; it must not be represented as a full-length 585-residue prediction.
- **T1183:** [CASP target page](https://predictioncenter.org/casp15/target.cgi?id=153&view=all), [RCSB entity 2](https://data.rcsb.org/rest/v1/core/polymer_entity/8IFX/2), [chain B mapping](https://data.rcsb.org/rest/v1/core/polymer_entity_instance/8IFX/B). TsaB author chain **B** exactly matches the 200-residue `targets/T1183/T1183.fasta`; author chain A is the different protein TsaD. There are 195 observed C-alpha residues at positions **1-195**; **196-200** have no coordinates. Author and label numbering both equal target numbering. The current CASP domain-summary row has an inconsistent `Res` column of 190; the explicit official evaluation interval is **1-195**, and the target page/FASTA length is 200. This configuration uses that explicit interval.

The best TBM PDB files named in each target's `metadata.json` were also inspected directly. They all use author chain **A** and preserve target residue numbering: T1147 has 103 C-alpha residues numbered 1-103, T1133 has 414 numbered 5-418, and T1183 has 200 numbered 1-200. Configured AF3 chain A refers to the local predictions; `af3_server` is unset.

## Alternate-coordinate policy

The T1133 reference selects alternate location **A** explicitly, before scoring.
Its C-alpha atoms at residues 5, 174, and 255 have tied A/B occupancies of 0.5.
The other references use the unique highest-occupancy C-alpha alternate at each
residue. In T1147 this selects A at 67 (0.8 versus 0.2) and B at 103 (0.51 versus
0.49). The selected alternate and occupancy are retained in the residue-mapping
table. Neither policy is selected by prediction accuracy.

## T1183 reference discrepancy

The current official CASP15 target and domain tables link **8IFX**, so that entry is fixed as the primary reference here. The structure-provider paper [Protein target highlights in CASP15: Analysis of models by structure providers](https://pmc.ncbi.nlm.nih.gov/articles/PMC10792529/) instead labels T1183 as **8IEY**. [8IEY](https://www.rcsb.org/structure/8IEY) and [8IFX](https://www.rcsb.org/structure/8IFX) are related Aquifex aeolicus TsaD-TsaB crystal structures; 8IFX is in complex with ADP. In both entries, entity 2/author chain B exactly matches the 200-residue target sequence, the resolution is 2.00 A, and positions 196-200 are unobserved. The choice of 8IFX follows the current official CASP mapping and was not chosen by comparing prediction accuracy. 8IEY was not downloaded into this reference set.

The files here are the identified current RCSB experimental coordinates, with the hashes below. Identity to the original CASP evaluation coordinate files has not been established. The [official public CASP target archive](https://predictioncenter.org/download_area/CASP15/targets/casp15.targets.TS-domains.public_12.20.2022.tar.gz) has not been inspected, so these references must not be described as byte-identical CASP-native files.

## Downloads and SHA-256

Download date for every file: **2026-10-06**. Reference files were retained without modification.

| Repository-relative file | Download URL | SHA-256 |
|---|---|---|
| `agent2/references/8EM5.cif` | [RCSB 8EM5 mmCIF](https://files.rcsb.org/download/8EM5.cif) | `41bced682b3a5b88754a7771eda126766bcca695036d6d2be744035f7dccc89b` |
| `agent2/references/8DYS.cif` | [RCSB 8DYS mmCIF](https://files.rcsb.org/download/8DYS.cif) | `4213eab7013a47d1fdd2c1392dbedf4013ccb3191ead916188e5f7d203b290e3` |
| `agent2/references/8IFX.cif` | [RCSB 8IFX mmCIF](https://files.rcsb.org/download/8IFX.cif) | `8a4a3eb36c50debb971e28365c07247271a862261c442107b8b2cad58f5157ec` |

Structure data citations: [8EM5 DOI](https://doi.org/10.2210/pdb8EM5/pdb), [8DYS DOI](https://doi.org/10.2210/pdb8DYS/pdb), [8IFX DOI](https://doi.org/10.2210/pdb8IFX/pdb). Relevant publications are [Cuthbert et al., 2024, MmpS5](https://doi.org/10.1093/mtomcs/mfae011) and [Lu et al., 2024, TsaD-TsaB](https://doi.org/10.1016/j.jbc.2024.107962). The RCSB 8DYS entry lists its structure citation as "To be published."
