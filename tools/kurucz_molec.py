"""
Parse Kurucz molecules.dat and evaluate ln K_eq(T).

Functional form decoded from NMOLEC in atlas12_modules.f90 (line ~6150):

  ln K_eq(T) = E1/(kT_eV) - E2 + E3*T - E4*T^2 + E5*T^3 - E6*T^4
               - 1.5 * (N_atoms - 2*ION - 1) * ln T

where:
  E1 = dissociation energy D0 [eV]
  E2..E6 = polynomial coefficients (units absorb Kurucz conventions)
  N_atoms = number of constituent atoms in the molecule
  ION = ionization state (number of electrons removed)

The last term is the Sackur-Tetrode translational contribution;
N_atoms - 2*ION - 1 counts the net reactant-minus-product translational
degrees of freedom in the dissociation reaction.

Species code encoding (from READMOL):
  Pairs of digits from left to right give atomic numbers.
  The decimal portion (.NN) gives the ionization state.

Examples:
   100.00      -> H I       (atom 1, +100 marker=neutral pseudo-atom, ION=0)
   101.00      -> H2        (atoms 1,1)
   101.01      -> H2+       (atoms 1,1; ION=1)
   606.00      -> C2
   607.00      -> CN
   608.00      -> CO
   822.00      -> TiO
   10108.00    -> H2O  (atoms 1,1,8)
   60808.00    -> CO2
   1010106.00  -> CH4  (atoms 1,1,1,6)
"""

import numpy as np

K_BOLTZMANN_EV = 8.617333262e-5  # Boltzmann constant in eV/K

# Element symbols indexed by atomic number
ELEMENT_SYMBOLS = [
    '',   'H',  'He', 'Li', 'Be', 'B',  'C',  'N',  'O',  'F',  'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P',  'S',  'Cl', 'Ar', 'K',  'Ca',
    'Sc', 'Ti', 'V',  'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y',  'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',
    'Sb', 'Te', 'I',  'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
    'Lu', 'Hf', 'Ta', 'W',  'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',
    'Pa', 'U'
]


def decode_species_code(code_real):
    """
    Decode a Kurucz real-valued species code into a list of atomic numbers
    plus an ionization state.

    Returns
    -------
    atoms : list of int
        Atomic numbers of constituent atoms (marker 100 = neutral pseudo-atom
        is preserved if present, indicating single atom rather than molecule).
    ion : int
        Number of electrons removed (0 = neutral molecule or atom).
    is_true_molecule : bool
        True if this is a polyatomic species (2+ real atoms).
    """
    # Integer part holds atom codes, fractional part holds ionization
    int_part = int(round(code_real * 100)) // 100  # integer part as int
    frac_part = int(round(code_real * 100)) - int_part * 100  # .NN as integer
    ion = frac_part

    # Peel off pairs of digits from the right
    atoms = []
    n = int_part
    while n > 0:
        pair = n % 100
        atoms.append(pair)
        n //= 100
    atoms.reverse()

    # True molecule if it has at least 2 real atoms (element 100 is marker)
    real_atoms = [a for a in atoms if 1 <= a <= 99]
    is_true_molecule = len(real_atoms) >= 2

    return atoms, ion, is_true_molecule


def species_formula(atoms, ion):
    """Build a human-readable formula from the atomic composition.

    Uses canonical astrochemistry ordering: for diatomic hydrides,
    metal/non-metal appears first (CH, OH, FeH), matching standard
    astronomical notation. For diatomics without H, the convention
    is heavier element last (CO, SiO, TiO). For triatomics and larger,
    H atoms appear first (H2O, H2S, H2CO), carbon-based molecules
    follow chemical formula convention (CO2, C2H2).
    """
    real_atoms = [a for a in atoms if 1 <= a <= 99]

    if not real_atoms:
        return '?'

    # Count each element
    counts = {}
    for a in real_atoms:
        counts[a] = counts.get(a, 0) + 1

    # Canonical ordering rules for common cases
    unique = sorted(counts.keys())

    def fmt(z, n):
        sym = ELEMENT_SYMBOLS[z] if z < len(ELEMENT_SYMBOLS) else f'Z{z}'
        return sym if n == 1 else f'{sym}{n}'

    ordered = []

    if len(unique) == 1:
        # Homonuclear: just one element
        z = unique[0]
        ordered = [fmt(z, counts[z])]
    elif len(unique) == 2:
        z1, z2 = unique  # z1 < z2
        n1, n2 = counts[z1], counts[z2]

        if z1 == 1 and n1 >= 2:
            # Hydrogen-bearing with 2+ H atoms: H first (H2O, H2S, H2CO, NH3, CH4)
            ordered = [fmt(z1, n1), fmt(z2, n2)]
        elif z1 == 1 and n1 == 1:
            # Single hydride: heavy element first (CH, OH, NH, FeH, SiH, MgH)
            ordered = [fmt(z2, n2), fmt(z1, n1)]
        elif z1 == 6 and z2 == 8:
            # Carbon-oxygen: CO, CO2
            ordered = [fmt(z1, n1), fmt(z2, n2)]
        elif z1 == 7 and z2 == 8:
            # Nitrogen-oxygen: NO, N2O, NO2
            ordered = [fmt(z1, n1), fmt(z2, n2)]
        elif z1 == 6 and z2 == 7:
            # CN, C2N, CN2
            ordered = [fmt(z1, n1), fmt(z2, n2)]
        elif z2 == 8:
            # Anything-oxide: heavy first for metals (SiO, TiO, FeO, CaO, MgO)
            # but chemistry convention is variable. Use z1 first (lighter).
            # For AlO, SiO, TiO, etc. the astronomy convention is heavy-first
            # (TiO, SiO, FeO) — put z2=8 last.
            ordered = [fmt(z1, n1), fmt(z2, n2)]
        elif z2 == 16:
            # Sulfides: MgS, AlS, SiS, FeS, CaS -- heavy first, S last
            ordered = [fmt(z1, n1), fmt(z2, n2)]
        else:
            # Default: lighter first
            ordered = [fmt(z1, n1), fmt(z2, n2)]
    else:
        # 3+ element types: list in order (H first if present, then by Z)
        if 1 in unique:
            ordered.append(fmt(1, counts[1]))
            for z in unique:
                if z != 1:
                    ordered.append(fmt(z, counts[z]))
        else:
            ordered = [fmt(z, counts[z]) for z in unique]

    formula = ''.join(ordered)

    if ion == 1:
        formula += '+'
    elif ion > 1:
        formula += f'{ion}+'

    return formula


# Canonical chemistry-convention names for diatomics, matching BC16 Table 1
# naming conventions. Astronomy convention (which matches BC16 and is the
# standard in stellar spectroscopy) places the heavier element first for
# oxides and metal hydrides/sulfides: TiO not OTi, CaH not HCa.
# Some pairings (CO, CN, NO, HF, HCl, SO) are universally lighter-first.
# This table overrides species_formula's defaults for the cases that matter.
_CANONICAL_DIATOMIC = {
    # Hydrogen halides: H first
    (1, 9):   'HF',
    (1, 17):  'HCl',
    # Oxides: heavy-Z first (matching BC16: MgO, AlO, SiO, CaO, TiO, VO, FeO)
    (8, 12):  'MgO',
    (8, 13):  'AlO',
    (8, 14):  'SiO',
    (8, 20):  'CaO',
    (8, 22):  'TiO',
    (8, 23):  'VO',
    (8, 26):  'FeO',
    # Sulfides: SO is the only oxide-like flip; metal sulfides use BC16 form
    (8, 16):  'SO',
    (16, 20): 'CaS',
    (16, 26): 'FeS',
    # C/N/Si/Al diatomics: BC16 alphabetical-ish convention (lower-Z first
    # except for CN, CO, etc. which are universal)
    (6, 13):  'AlC',
    (6, 14):  'SiC',
    (7, 13):  'AlN',
    (7, 14):  'SiN',
    (7, 12):  'MgN',     # not in BC16, but 'MgN' matches the convention
    # Note: CO, CN, NO, CS, NS keep species_formula default (lighter-first),
    # which is correct for these.
}


# Canonical chemistry-convention names for polyatomics, indexed by the
# sorted multiset of atomic numbers. These override species_formula's
# default ordering for entries where the universal convention disagrees
# with simple "H-first then by Z" sorting (e.g. NH3 vs H3N, CaOH vs HOCa).
_CANONICAL_POLYATOMIC = {
    (1, 7, 7):           'N2H',     # not in our set, placeholder
    (1, 1, 7):           'NH2',
    (1, 1, 1, 7):        'NH3',
    (1, 1, 1, 1, 6):     'CH4',
    (1, 1, 1, 6):        'CH3',
    (1, 1, 6):           'CH2',
    (1, 6, 7):           'HCN',
    (1, 6, 8):           'HCO',
    (1, 7, 8):           'HNO',
    (1, 8, 8):           'HO2',
    (1, 8, 11):          'NaOH',
    (1, 8, 12):          'MgOH',
    (1, 8, 20):          'CaOH',
    (1, 1, 1, 1, 14):    'SiH4',
    (1, 1, 1, 1, 1):     'H3+',     # ionic, special-cased below
    (6, 6, 6):           'C3',
    (6, 6, 7, 7):        'C2N2',
    (6, 6, 8):           'C2O',
    (6, 6, 14):          'C2Si',
    (6, 8, 8):           'CO2',
    (6, 8, 16):          'OCS',
    (6, 14, 14):         'CSi2',
    (6, 16, 16):         'CS2',
    (7, 7, 8):           'N2O',
    (7, 8, 8):           'NO2',
    (8, 8, 14):          'SiO2',
    (8, 8, 16):          'SO2',
    (1, 1, 6, 6):        'C2H2',
    (1, 1, 6):           'CH2',
    (1, 1, 8):           'H2O',
    (1, 1, 16):          'H2S',
    (1, 6):              'CH',
    (1, 6, 6):           'C2H',
}


# Roman numeral converter for atomic ionization labels (HI, HII, FeXXVI, ...)
def _roman(n):
    """Convert integer 1..99 to a Roman numeral string."""
    table = [(50, 'L'), (40, 'XL'), (10, 'X'), (9, 'IX'),
             (5, 'V'),  (4, 'IV'), (1, 'I')]
    if n < 1:
        return ''
    out = []
    for value, sym in table:
        while n >= value:
            out.append(sym)
            n -= value
    return ''.join(out)


def species_label(atoms, ion):
    """
    Build a human-readable label for an atom or molecule.

    Atoms (a single real-element atom): astronomy Roman-numeral notation
        H + neutral  -> 'HI'
        H + ion=1    -> 'HII'
        Fe + ion=2   -> 'FeIII'

    Anions (last component is the element-100 "extra electron" marker
    in Kurucz's convention): get a '-' suffix.
        [1, 100]    -> 'H-'
        [1, 1, 100] -> 'H2-'
        [6, 7, 100] -> 'CN-'

    Molecules (>=2 real atoms): standard chemistry-convention names.
        Looked up in _CANONICAL_POLYATOMIC for triatomics+, otherwise
        fall through to species_formula() which already handles diatomics.
        Singly-ionized molecules get a '+' suffix (OH+, H2+, CN+).
        2+ ionized molecules get '<n>+' suffix.
    """
    real_atoms = [a for a in atoms if 1 <= a <= 99]
    if not real_atoms:
        return '?'

    # Detect anion: a trailing element-100 marker means an extra electron
    # is attached. The actual chemical species is the rest.
    is_anion = (len(atoms) > 0 and atoms[-1] == 100)

    # Single real-atom stub (no anion marker): Roman-numeral atom notation.
    if len(real_atoms) == 1 and not is_anion:
        z = real_atoms[0]
        sym = ELEMENT_SYMBOLS[z] if z < len(ELEMENT_SYMBOLS) else f'Z{z}'
        return f'{sym}{_roman(ion + 1)}'

    # Single real-atom anion: e.g. [1, 100] -> H-
    if len(real_atoms) == 1 and is_anion:
        z = real_atoms[0]
        sym = ELEMENT_SYMBOLS[z] if z < len(ELEMENT_SYMBOLS) else f'Z{z}'
        return f'{sym}-'

    # Molecular species: try canonical lookup first
    key = tuple(sorted(real_atoms))
    if key in _CANONICAL_POLYATOMIC:
        base = _CANONICAL_POLYATOMIC[key]
    elif key in _CANONICAL_DIATOMIC:
        base = _CANONICAL_DIATOMIC[key]
    else:
        # Fall back to species_formula for entries not in either table,
        # passing ion=0 to get the bare formula (we add the suffix below).
        base = species_formula(real_atoms, 0)

    if is_anion:
        return base + '-'
    if ion == 1:
        return base + '+'
    elif ion > 1:
        return f'{base}{ion}+'
    return base


def parse_molecules_dat(path):
    """
    Parse Kurucz molecules.dat.

    Supports both:
      Old format (cols 1-18 = code, 19-25 = E1, 26-80 = E2..E6).
      New format (cols 1-10 = human label, 11-28 = code, 29-35 = E1,
                  36-90 = E2..E6).

    Format is auto-detected per-line by checking whether cols 1-18 of
    the trimmed line contain a parseable float.

    Returns list of dicts, one per entry.
    """
    entries = []
    with open(path, 'r') as f:
        for line_no, line in enumerate(f, 1):
            # Strip trailing newline; preserve internal whitespace
            stripped = line.rstrip('\n')
            # Skip blank lines
            if not stripped.strip():
                continue
            # Skip Fortran-style comment lines (after optional leading whitespace)
            if stripped.lstrip().startswith('!'):
                continue

            # Auto-detect format: try cols 1-18 first.
            # If they parse as a float, it's the old format (code in 1-18).
            # If they don't, it's the new format with a 10-char label
            # followed by the code in cols 11-28.
            label = ''
            old_format = True
            head18 = stripped[0:18] if len(stripped) >= 18 else stripped.ljust(18)
            try:
                float(head18.strip())
                old_format = True
            except ValueError:
                # Not a number in cols 1-18 -> must be new format with label
                old_format = False
                label = stripped[0:10].strip()

            if old_format:
                # Old layout: code 1-18, E1 19-25, E2..E6 in 26-80
                if len(stripped) < 18:
                    stripped = stripped.ljust(18)
                code_str = stripped[0:18]
                e1_str = stripped[18:25] if len(stripped) > 18 else ''
                e2_str = stripped[25:36] if len(stripped) > 25 else ''
                e3_str = stripped[36:47] if len(stripped) > 36 else ''
                e4_str = stripped[47:58] if len(stripped) > 47 else ''
                e5_str = stripped[58:69] if len(stripped) > 58 else ''
                e6_str = stripped[69:80] if len(stripped) > 69 else ''
            else:
                # New compact layout (90 cols total):
                #   label 1-10, code 11-22, gap 23, E1 24-30, E2..E6 in 31-90
                #   (5 fields of 12 chars each: E12.4)
                if len(stripped) < 22:
                    stripped = stripped.ljust(22)
                code_str = stripped[10:22]
                e1_str = stripped[23:30] if len(stripped) > 23 else ''
                e2_str = stripped[30:42] if len(stripped) > 30 else ''
                e3_str = stripped[42:54] if len(stripped) > 42 else ''
                e4_str = stripped[54:66] if len(stripped) > 54 else ''
                e5_str = stripped[66:78] if len(stripped) > 66 else ''
                e6_str = stripped[78:90] if len(stripped) > 78 else ''

            try:
                code = float(code_str)
            except ValueError:
                continue

            # Legacy terminator (still recognized)
            if code == 0.0:
                break

            def parse_float(s, default=0.0):
                s = s.strip()
                if not s:
                    return default
                try:
                    return float(s)
                except ValueError:
                    return default

            e1 = parse_float(e1_str)
            e2 = parse_float(e2_str)
            e3 = parse_float(e3_str)
            e4 = parse_float(e4_str)
            e5 = parse_float(e5_str)
            e6 = parse_float(e6_str)

            atoms, ion, is_true_molecule = decode_species_code(code)
            # Apply the same remapping READMOL does: a decoded atomic number
            # of 0 (from a "00" digit pair) represents element 100, the
            # neutral-atom pool marker.
            atoms = [100 if a == 0 else a for a in atoms]
            real_atoms = [a for a in atoms if 1 <= a <= 99]

            # Standard molecule: all components are real atoms (no element-100
            # "neutral atom pool" markers embedded). These are the entries
            # that correspond to dissociation equilibria comparable to
            # modern references (B&C 2016, ExoMol). Entries with an embedded
            # element-100 marker (e.g. code 10100.00 = H + H + atom_pool)
            # are auxiliary equations for the chemical equilibrium solver
            # and have a different physical interpretation.
            has_atom_pool_marker = any(a == 100 for a in atoms)
            # is_true_molecule needs re-evaluation now
            is_true_molecule = len(real_atoms) >= 2
            is_standard_molecule = is_true_molecule and not has_atom_pool_marker

            formula = species_formula(atoms, ion)

            # N_atoms in the dissociation reaction as used by NMOLEC:
            # NCOMP = number of entries in KCOMPS (atom IDs decoded from
            # the species code, with 0-pairs mapped to the atom-pool
            # sentinel 100, plus ION electron markers (101) for ionized
            # species).
            # In the ln T Sackur-Tetrode term: -1.5*(NCOMP - 2*ION - 1).
            # NCOMP = len(atoms) + ion (atom-pool sentinels are counted),
            # so NCOMP - 2*ION - 1 = len(atoms) - ION - 1.
            # Verify against F90 NCOMP-based runtime:
            #   H2     (atoms=[1,1],   ION=0): n_trans=1
            #   H2+    (atoms=[1,1],   ION=1): n_trans=0
            #   H2O    (atoms=[1,1,8], ION=0): n_trans=2
            #   C-     (atoms=[6,100], ION=0): n_trans=1  (atom + electron)
            #   C2-    (atoms=[6,6,100], ION=0): n_trans=2  (atom+atom+electron)
            #   atomic C (atoms=[6],   ION=0): n_trans=0  (no dissociation)
            n_real_atoms = len(real_atoms)
            n_trans = len(atoms) - ion - 1

            entries.append({
                'line': line_no,
                'code': code,
                'atoms': atoms,
                'real_atoms': real_atoms,
                'ion': ion,
                'n_real_atoms': n_real_atoms,
                'n_trans': n_trans,
                'is_true_molecule': is_true_molecule,
                'is_standard_molecule': is_standard_molecule,
                'has_atom_pool_marker': has_atom_pool_marker,
                'formula': formula,
                'E': (e1, e2, e3, e4, e5, e6),
            })

    return entries


def kurucz_lnKeq(entry, T):
    """
    Evaluate Kurucz's ln K_eq(T) for a molecule entry.

    Matches NMOLEC's expression exactly:
      ln K = E1/kT - E2 + E3 T - E4 T^2 + E5 T^3 - E6 T^4 - 1.5*n_trans*ln T

    Parameters
    ----------
    entry : dict
        One molecule entry from parse_molecules_dat.
    T : float or ndarray
        Temperature in Kelvin.
    """
    T = np.asarray(T, dtype=float)
    E1, E2, E3, E4, E5, E6 = entry['E']
    n_trans = entry['n_trans']

    if E1 == 0.0:
        # Entries with E1=0 are atoms (not real molecules) - no polynomial eval
        return np.full_like(T, np.nan, dtype=float)

    TkeV = K_BOLTZMANN_EV * T
    lnK = (E1 / TkeV
           - E2
           + E3 * T
           - E4 * T**2
           + E5 * T**3
           - E6 * T**4
           - 1.5 * n_trans * np.log(T))
    return lnK


if __name__ == '__main__':
    entries = parse_molecules_dat('molecules.dat')
    atoms_only = [e for e in entries if not e['is_true_molecule']]
    std_mols = [e for e in entries if e['is_standard_molecule']]
    aux_eqs = [e for e in entries if e['is_true_molecule'] and e['has_atom_pool_marker']]

    print(f'Total entries:                    {len(entries)}')
    print(f'Atom-only entries:                {len(atoms_only)}  (E1=0, ionization bookkeeping)')
    print(f'Standard molecules:               {len(std_mols)}  (comparable to BC16/ExoMol)')
    print(f'Auxiliary eq. entries:            {len(aux_eqs)}  (have atom-pool marker, different physics)')
    print()
    print('Standard molecules (candidates for BC16/ExoMol comparison):')
    print(f'  {"code":>12s}  {"formula":<10s}  {"N":>2s}  {"ION":>3s}  {"D0[eV]":>8s}')
    for e in std_mols:
        E1 = e['E'][0]
        print(f'  {e["code"]:12.2f}  {e["formula"]:<10s}  '
              f'{e["n_real_atoms"]:2d}  {e["ion"]:3d}  {E1:8.4f}')

    print()
    print('Auxiliary equation entries (contain element-100 marker):')
    print(f'  {"code":>12s}  {"formula":<10s}  {"D0[eV]":>8s}  {"E2":>7s}')
    for e in aux_eqs:
        E1, E2 = e['E'][0], e['E'][1]
        print(f'  {e["code"]:12.2f}  {e["formula"]:<10s}  {E1:8.4f}  {E2:7.2f}')
