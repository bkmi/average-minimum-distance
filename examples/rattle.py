import amd
import numpy as np

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Structure

from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
from pymatgen.io.cif import CifWriter

structure = Structure.from_file("CaCO3.cif")
adaptor = AseAtomsAdaptor()
atoms = adaptor.get_atoms(structure)
rattle_levels = [0.0, 0.01, 0.05, 0.1, 0.5, 1.0, 1.5, 2.0]  # Angstroms, for example
rattled_structures = []
cifs = []
crystals = []
for level in rattle_levels:
    atoms_copy = atoms.copy()
    atoms_copy.rattle(stdev=level, seed=42)  # seed for reproducibility
    rattled_structure = adaptor.get_structure(atoms_copy)
    rattled_structures.append(rattled_structure)
    cif_str = str(CifWriter(rattled_structure))
    cifs.append(cif_str)
    crystals.extend(list((amd.StrCifReader(cif_str)))) # read crystals

pdds = [amd.PDD(crystal, 100) for crystal in crystals] # calculate PDDs (k=100)
d = [amd.EMD(pdds[0], pdd) for pdd in pdds[1:]] # (Earth mover's) distance between 1st and 2nd PDDs

matcher = StructureMatcher()
rmsds = []
print(f"{'Rattle level':<20} {'RMSD with original':<25} {'EMD distance':<15}")
for i, (rattled_structure, distance) in enumerate(zip(rattled_structures[1:], d)):
    try:
        rmsd, _ = matcher.get_rms_dist(rattled_structures[0], rattled_structure)
    except:
        rmsd = np.nan
    rmsds.append(rmsd)
    print(f"{rattle_levels[i+1]:<20.2f} {rmsds[-1]:<25} {distance:<15}")
