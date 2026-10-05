"""Copy template mmCIFs from the local PDB library and write chain PDB files.

For each PDB ID, the mmCIF is copied into the output directory and the
requested protein chain is written as <id>_<chain>.pdb (protein residues
only, first altloc kept, no waters or ligands) for use with MODELLER.
"""

import argparse
import os
import shutil

from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

MMCIF_DIR = "/cluster/VAST/cietest/alphafold_database/pdb_mmcif/mmcif_files"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M",
}


def cif_path(pdb_id):
    return os.path.join(MMCIF_DIR, f"{pdb_id.lower()}.cif")


def _as_list(value):
    return value if isinstance(value, list) else [value]


def read_chain_residues(cif_file, chain):
    """Return the SEQRES sequence and observed residues of an author chain.

    Returns (seqres, observed) where seqres is the one-letter entity
    sequence and observed maps label_seq_id (1-based SEQRES position) to
    (auth_seq_id, ins_code, one_letter) for residues with coordinates in
    the first model.
    """
    d = MMCIF2Dict(cif_file)

    group = _as_list(d["_atom_site.group_PDB"])
    auth_asym = _as_list(d["_atom_site.auth_asym_id"])
    label_seq = _as_list(d["_atom_site.label_seq_id"])
    auth_seq = _as_list(d["_atom_site.auth_seq_id"])
    icode = _as_list(d["_atom_site.pdbx_PDB_ins_code"])
    comp = _as_list(d["_atom_site.label_comp_id"])
    entity = _as_list(d["_atom_site.label_entity_id"])
    model = _as_list(d["_atom_site.pdbx_PDB_model_num"])

    observed = {}
    entity_id = None
    first_model = model[0]
    for i in range(len(group)):
        if auth_asym[i] != chain or model[i] != first_model or label_seq[i] in (".", "?"):
            continue
        entity_id = entity[i]
        pos = int(label_seq[i])
        if pos not in observed:
            ins = "" if icode[i] in (".", "?") else icode[i]
            observed[pos] = (int(auth_seq[i]), ins, THREE_TO_ONE.get(comp[i], "X"))

    if entity_id is None:
        raise ValueError(f"No polymer residues for chain {chain} in {cif_file}")

    ent = _as_list(d["_entity_poly_seq.entity_id"])
    num = _as_list(d["_entity_poly_seq.num"])
    mon = _as_list(d["_entity_poly_seq.mon_id"])
    residues = {}
    for e, n, m in zip(ent, num, mon):
        if e == entity_id:
            residues.setdefault(int(n), THREE_TO_ONE.get(m, "X"))
    seqres = "".join(residues[n] for n in sorted(residues))
    return seqres, observed


class ChainProteinSelect(Select):
    """Keep standard protein residues of one chain and the first altloc."""

    def __init__(self, chain):
        self.chain = chain

    def accept_model(self, model):
        return model.id == 0

    def accept_chain(self, chain):
        return chain.id == self.chain

    def accept_residue(self, residue):
        hetflag = residue.id[0]
        return hetflag == " " or residue.get_resname() == "MSE"

    def accept_atom(self, atom):
        if atom.element in ("H", "D"):
            return False
        if atom.is_disordered() or atom.get_altloc() != " ":
            keep = atom.get_altloc() in ("A", "1")
            if keep:
                atom.set_altloc(" ")
            return keep
        return True


def write_chain_pdb(pdb_id, chain, out_dir):
    """Copy the mmCIF and write <id>_<chain>.pdb; return the PDB path."""
    src = cif_path(pdb_id)
    shutil.copy2(src, os.path.join(out_dir, os.path.basename(src)))

    structure = MMCIFParser(QUIET=True).get_structure(pdb_id.lower(), src)
    for residue in structure[0][chain]:
        if residue.get_resname() == "MSE":
            residue.resname = "MET"
            residue.id = (" ", residue.id[1], residue.id[2])
            for atom in residue:
                if atom.element == "SE":
                    atom.fullname, atom.element = " SD ", "S"

    out = os.path.join(out_dir, f"{pdb_id.lower()}_{chain}.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save(out, ChainProteinSelect(chain))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("templates", nargs="+",
                        help="Templates as PDBID_CHAIN, e.g. 4EHX_A")
    parser.add_argument("-o", "--out-dir", default="templates")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    for tmpl in args.templates:
        pdb_id, chain = tmpl.split("_")
        out = write_chain_pdb(pdb_id, chain, args.out_dir)
        seqres, observed = read_chain_residues(cif_path(pdb_id), chain)
        missing = len(seqres) - len(observed)
        print(f"{tmpl}: {out}  SEQRES={len(seqres)}  observed={len(observed)}  "
              f"missing={missing}")


if __name__ == "__main__":
    main()
