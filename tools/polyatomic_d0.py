"""
polyatomic_d0.py -- verified atomization energies for ExoMol polyatomics.

Reconstruction of the April 2026 module deleted after commit bb386c9.
The D0 values are recovered from the April fits by exact-digit
inversion (a1 is degenerate with D0/k_B, so the printed a1 of each
bb386c9 row pins the assembly D0 to ~5e-6 eV); provenance strings are
16 "ExoMol fit" rows of bb386c9's data/molecules.dat (D0 column +
trailing comments); atoms lists follow the species codes.  Keyed by
ExoMol .pf basename.  H3+ has a local .pf but no entry -- the April
pipeline handled neutrals only (ion == 0).
"""

POLYATOMIC_D0 = {
    'H2O':  {'exomol_name': 'H2O',  'atoms': ['H', 'H', 'O'],
             'D0': 9.51129,  'source': 'CCCBDB (2006Rus/Pin)', 'verified': True},
    'H2S':  {'exomol_name': 'H2S',  'atoms': ['H', 'H', 'S'],
             'D0': 7.50995,  'source': 'CCCBDB (CODATA)',      'verified': True},
    'HCN':  {'exomol_name': 'HCN',  'atoms': ['H', 'C', 'N'],
             'D0': 13.117, 'source': 'CCCBDB (Gurvich)',     'verified': True},
    'NaOH': {'exomol_name': 'NaOH', 'atoms': ['H', 'O', 'Na'],
             'D0': 7.85612,  'source': 'CCCBDB (JANAF)',       'verified': True},
    'CaOH': {'exomol_name': 'CaOH', 'atoms': ['H', 'O', 'Ca'],
             'D0': 8.396095,  'source': 'CCCBDB (Gurvich)',     'verified': True},
    'C3':   {'exomol_name': 'C3',   'atoms': ['C', 'C', 'C'],
             'D0': 13.67275, 'source': 'ATcT1.220',            'verified': True},
    'CO2':  {'exomol_name': 'CO2',  'atoms': ['C', 'O', 'O'],
             'D0': 16.561065, 'source': 'CCCBDB (CODATA)',      'verified': True},
    'OCS':  {'exomol_name': 'OCS',  'atoms': ['O', 'C', 'S'],
             'D0': 14.24776, 'source': 'CCCBDB (Gurvich)',     'verified': True},
    'N2O':  {'exomol_name': 'N2O',  'atoms': ['N', 'N', 'O'],
             'D0': 11.435935, 'source': 'CCCBDB (Gurvich)',     'verified': True},
    'SiO2': {'exomol_name': 'SiO2', 'atoms': ['O', 'O', 'Si'],
             'D0': 12.78443, 'source': 'ATcT1.220',            'verified': True},
    'SO2':  {'exomol_name': 'SO2',  'atoms': ['O', 'O', 'S'],
             'D0': 11.01515, 'source': 'CCCBDB (CODATA)',      'verified': True},
    'CH3':  {'exomol_name': 'CH3',  'atoms': ['H', 'H', 'H', 'C'],
             'D0': 12.533725, 'source': 'ATcT1.220',            'verified': True},
    'NH3':  {'exomol_name': 'NH3',  'atoms': ['H', 'H', 'H', 'N'],
             'D0': 11.998715, 'source': 'CCCBDB (CODATA)',      'verified': True},
    'C2H2': {'exomol_name': 'C2H2', 'atoms': ['H', 'H', 'C', 'C'],
             'D0': 16.856445, 'source': 'CCCBDB (Gurvich)',     'verified': True},
    'CH4':  {'exomol_name': 'CH4',  'atoms': ['H', 'H', 'H', 'H', 'C'],
             'D0': 17.016, 'source': 'CCCBDB (Gurvich)',     'verified': True},
    'SiH4': {'exomol_name': 'SiH4', 'atoms': ['H', 'H', 'H', 'H', 'Si'],
             'D0': 13.132155, 'source': 'ATcT1.220',            'verified': True},
    # --- July 2026 additions (D0 ledger in the adding commit; "JANAF
    # chain" = raw JANAF dHf(0 K) of molecule and atoms, janaf_loader) ---
    'KOH':  {'exomol_name': 'KOH',  'atoms': ['H', 'O', 'K'],
             'D0': 8.096,  'source': 'JANAF chain',              'verified': True},
    'LiOH': {'exomol_name': 'LiOH', 'atoms': ['H', 'O', 'Li'],
             'D0': 8.836,  'source': 'JANAF chain',              'verified': True},
    'SiH2': {'exomol_name': 'SiH2', 'atoms': ['H', 'H', 'Si'],
             'D0': 6.270,  'source': 'Berkowitz & Ruscic 1987',  'verified': True},
    'HBO':  {'exomol_name': 'HBO',  'atoms': ['H', 'B', 'O'],
             'D0': 12.701, 'source': 'JANAF dHf298 + ATcT B(g)', 'verified': False},
}
