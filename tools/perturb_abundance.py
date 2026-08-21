#!/usr/bin/env python3
"""Raise one element's abundance in an ATLAS12 .atm by a given amount.

Used to measure a curve-of-growth slope with the code itself: synthesise
twice at FIXED structure, differing only in A(X), and dEW/dlogA follows.

THE TRAP THIS EXISTS TO AVOID: an .atm carries the abundances TWICE -- as
`ABUNDANCE CHANGE` cards and again as an `ABUNDANCE TABLE` block -- and
READIN applies them in file order, so the TABLE (which comes last) wins.
Editing only the CHANGE cards changes nothing at all, silently.  Both are
rewritten here.  The TABLE is fixed-format, 5 entries of 20 columns per line
(7 label, F7.3 value, F6.3 offset), so it is edited by column, not by regex.

  python3 tools/perturb_abundance.py in.atm out.atm --z 12 --delta 0.1
"""
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("infile"); ap.add_argument("outfile")
    ap.add_argument("--z", type=int, required=True)
    ap.add_argument("--delta", type=float, required=True)
    a = ap.parse_args()

    out, in_table, nchange, ntable = [], False, 0, 0
    for ln in open(a.infile):
        s = ln.rstrip("\n")
        if s.strip().startswith("ABUNDANCE TABLE"):
            in_table = True
            out.append(s); continue
        if in_table:
            # the table ends at the first card that is not a table row
            if s.strip().startswith("READ") or "ABUNDANCE" in s:
                in_table = False
            else:
                new = ""
                for i in range(0, len(s), 20):
                    f = s[i:i+20]
                    if len(f) >= 14:
                        lab = f[:7]
                        try: z = int("".join(c for c in lab if c.isdigit()))
                        except ValueError: z = -1
                        if z == a.z and z > 2:
                            v = float(f[7:14]) + a.delta
                            f = lab + f"{v:7.3f}" + f[14:]
                            ntable += 1
                    new += f
                out.append(new); continue
        if s.strip().startswith("ABUNDANCE CHANGE") and " CHANGE" in s:
            head, _, rest = s.partition("CHANGE")
            t = rest.split()
            for i in range(0, len(t)-1, 2):
                if int(t[i]) == a.z and a.z > 2:
                    v = float(t[i+1]) + a.delta
                    # keep the original column layout: same width field
                    old = f"{t[i]:>2s}{t[i+1]:>7s}"
                    new = f"{t[i]:>2s}{v:7.2f}"
                    if old in s:
                        s = s.replace(old, new, 1); nchange += 1
        out.append(s)

    open(a.outfile, "w").write("\n".join(out) + "\n")
    print(f"Z={a.z} {a.delta:+.2f} dex: {ntable} table entr(y/ies), "
          f"{nchange} change card entr(y/ies) rewritten -> {a.outfile}")
    if ntable == 0:
        raise SystemExit("ERROR: nothing changed in the ABUNDANCE TABLE -- "
                         "that block is the one that takes effect")


if __name__ == "__main__":
    main()
