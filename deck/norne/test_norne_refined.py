"""Check the face CoSTA API on Norne, with a vertically refined truth.

The truth is Norne with every layer split in two (refine_vertical.py),
restricted mass consistently to the original grid. The wells differ
between the grids, so face CoSTA also gets a cell correction for what
face fluxes cannot represent, which should be localized at the wells.

Usage: python test_norne_refined.py <directory with NORNE_ATW2013.DATA>
"""

import os
import sys
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from costa_check import Restriction, run  # noqa: E402
from costa_remote import Remote  # noqa: E402
import refine_vertical  # noqa: E402

DECK = 'NORNE_ATW2013.DATA'
NSTEPS = 8
R = 2

# CoSTA assumes accurate solves, far below flow's default tolerances.
# Only the states in memory are used, so no ECL output files are written.
ARGS = ['--enable-ecl-output=false', '--linear-solver-reduction=1e-12',
        '--tolerance-cnv=1e-10', '--tolerance-cnv-relaxed=1e-10',
        '--tolerance-mb=1e-14', '--tolerance-mb-relaxed=1e-14',
        '--newton-max-iterations=30']


def main(norne):
    with tempfile.TemporaryDirectory() as tmp:
        refine_vertical.main(os.path.join(norne, DECK), tmp, R)
        truth = Remote('BlackOil', os.path.join(tmp, DECK), ARGS)
        face = Remote('BlackOil', os.path.join(norne, DECK), ARGS)
        rhs = Remote('BlackOil', os.path.join(norne, DECK), ARGS)
        try:
            restriction = Restriction(truth, face, (1, 1, R))
            # Ambiguous upwinding where a phase potential vanishes gives
            # flux errors up to about 1e-3, see test_norne.py.
            ok, _ = run('Norne, refined truth', truth, face, rhs, NSTEPS, restriction,
                        hybrid=True, flux_tol=1e-3, require_switching=False)
        finally:
            for sim in (truth, face, rhs):
                sim.close()
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main(os.path.abspath(sys.argv[1])))
