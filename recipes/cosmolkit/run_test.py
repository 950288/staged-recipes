"""Exercise the installed native extension without the source checkout or RDKit."""

import cosmolkit as ck
import numpy as np


mol = ck.Molecule.from_smiles("CCO")
assert mol.num_atoms() == 3
assert mol.to_smiles() == "CCO"

hydrogenated = mol.with_hydrogens()
assert hydrogenated.num_atoms() == 9
assert mol.num_atoms() == 3

restored = ck.Molecule.from_binary(mol.to_binary())
assert restored.num_atoms() == 3
assert restored.to_smiles() == "CCO"

fingerprint = mol.fingerprint_morgan()
bits = fingerprint.to_numpy()
assert bits.dtype == np.uint8
assert bits.shape == (2048,)
assert np.isin(bits, [0, 1]).all()
assert bits.any()

query = ck.QueryGraph.from_smarts("CO")
assert mol.substruct_matches(query)

print("COSMolKit: parsing, COW hydrogens, native archives, NumPy, and search passed")
