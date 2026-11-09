"""Check the face CoSTA API on Norne, with a refined truth.

The truth is Norne refined, by default with every layer split in two
(refine_vertical.py), restricted mass consistently to the original grid.
The wells differ between the grids, so face CoSTA also gets a cell
correction for what face fluxes cannot represent, which should be
localized at the wells.

Usage: python test_norne_refined.py <directory with NORNE_ATW2013.DATA>
           [--steps N] [--fine DECK --ratio RX RY RZ]

--fine takes a refined deck made elsewhere, e.g. by opmt regrid, nested in
the original grid with the given refinement ratios.
"""

import argparse
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

# CoSTA assumes accurate solves, below flow's default tolerances. Tighter
# than this, the refined model cannot be run for long: near round-off the
# Newton iterations stall and the steps are chopped until they fail.
# Only the states in memory are used, so no ECL output files are written.
ARGS = ['--enable-ecl-output=false', '--linear-solver-reduction=1e-10',
        '--tolerance-cnv=1e-8', '--tolerance-cnv-relaxed=1e-6',
        '--tolerance-mb=1e-12', '--tolerance-mb-relaxed=1e-10',
        '--newton-max-iterations=30']


def main(norne, nsteps=NSTEPS, fine=None, ratio=(1, 1, R)):
    with tempfile.TemporaryDirectory() as tmp:
        if fine is None:
            refine_vertical.main(os.path.join(norne, DECK), tmp, R)
            fine = os.path.join(tmp, DECK)
        truth = Remote('BlackOil', fine, ARGS)
        face = Remote('BlackOil', os.path.join(norne, DECK), ARGS)
        rhs = Remote('BlackOil', os.path.join(norne, DECK), ARGS)
        try:
            restriction = Restriction(truth, face, ratio)
            # Ambiguous upwinding where a phase potential vanishes gives
            # flux errors up to about 1e-3, see test_norne.py.
            ok, _ = run('Norne, refined truth', truth, face, rhs, nsteps, restriction,
                        hybrid=True, flux_tol=1e-3, require_switching=False)
        finally:
            for sim in (truth, face, rhs):
                sim.close()
    return 0 if ok else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('norne', help='directory with NORNE_ATW2013.DATA')
    parser.add_argument('--steps', type=int, default=NSTEPS)
    parser.add_argument('--fine', help='refined deck nested in the original')
    parser.add_argument('--ratio', type=int, nargs=3, default=(1, 1, R),
                        metavar=('RX', 'RY', 'RZ'))
    a = parser.parse_args()
    raise SystemExit(main(os.path.abspath(a.norne), a.steps,
                          a.fine and os.path.abspath(a.fine), tuple(a.ratio)))
