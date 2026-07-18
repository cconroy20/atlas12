"""
matcher.py -- match Kurucz molecules.dat entries to BC16 species names.

Reconstruction of the April 2026 module deleted after commit bb386c9.
The match set is validated by reproducing bb386c9's molecules.dat: the
49 rows tagged "BC16 fit" there (41 neutral diatomics + 8 molecular
ions) define exactly which entries must match.

match_kurucz_to_bc16(entries, bc16) -> [(entry, bc16_name | None, D0)]

Matching is by unordered atom multiset (handles Kurucz "SH" vs BC16
"HS") plus charge state.  Entries with the atom-pool marker (anions
like C2-, whose components include pseudo-element 100) never match:
BC16's anion data uses a different reaction convention.
"""

import re

# Element symbols indexed by atomic number (subset is enough: BC16
# diatomics only involve Z <= 92 with standard symbols)
from kurucz_molec import ELEMENT_SYMBOLS

_SYM_TO_Z = {s: z for z, s in enumerate(ELEMENT_SYMBOLS) if s}


def _bc16_atoms(name):
    """BC16 species name -> (sorted atom Z list, charge) or None."""
    charge = name.count('+')
    core = name.rstrip('+')
    if core.endswith('-'):
        return None                       # anions: not matched
    atoms = []
    for sym, count in re.findall(r'([A-Z][a-z]?)(\d?)', core):
        if not sym:
            continue
        if sym not in _SYM_TO_Z:
            return None
        atoms += [_SYM_TO_Z[sym]] * (int(count) if count else 1)
    return sorted(atoms), charge


def match_kurucz_to_bc16(entries, bc16):
    """One (entry, name-or-None, D0-or-None) tuple per input entry."""
    # Index BC16 species by (atom multiset, charge)
    index = {}
    for name in bc16:
        key = _bc16_atoms(name)
        if key is not None:
            index[(tuple(key[0]), key[1])] = name

    out = []
    for e in entries:
        name = None
        if (e['n_real_atoms'] == 2 and not e['has_atom_pool_marker']):
            key = (tuple(sorted(e['real_atoms'])), e['ion'])
            name = index.get(key)
        d0 = bc16[name].D0_adopted if name is not None else None
        out.append((e, name, d0))
    return out
