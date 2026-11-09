"""Check the face CoSTA API on synthetic blackoil decks based on SPE1.

Two scenarios, with decks from make_decks.py:

 A. Operator error: truth FINE, model FINE_PERTURBED (perturbed transmissibilities).
    The exact labels satisfy delta* = -D tau*, so face CoSTA must recover the
    truth state and fluxes, while RHS CoSTA only recovers the state.
 B. Grid error: truth FINE restricted to the COARSE grid. States are restricted
    mass consistently: the coarse state holds the summed fine component
    amounts. The labels are then consistent except at the producer, whose
    well model differs between the grids. Face fluxes cannot represent that,
    so face CoSTA also gets the cell correction sigma* = delta* + D tau*,
    which must be localized, and must then recover the state and fluxes.

Gas appears and redissolves, so the corrections are solved with primary
variable switching. "switched" counts cells where gas is present in the
corrected state but not in the prediction p0 (the Newton initial guess), or
vice versa, and "over the step" the same relative to the start of the step.

Usage: python test_blackoil.py <path to SPE1CASE2.DATA>
"""

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from costa_check import Restriction, run  # noqa: E402
from costa_remote import Remote  # noqa: E402
import make_decks  # noqa: E402

NSTEPS = 6

# CoSTA assumes accurate solves, far below flow's default tolerances.
# Only the states in memory are used, so no ECL output files are written.
ARGS = ['--enable-ecl-output=false', '--linear-solver-reduction=1e-10',
        '--tolerance-cnv=1e-8', '--tolerance-cnv-relaxed=1e-8',
        '--tolerance-mb=1e-12', '--tolerance-mb-relaxed=1e-12',
        '--newton-max-iterations=30']


def robustness(truth, face, ut_prev):
    """A correction far too large to converge must not fail the step. It is
    applied partially, with the failed cells reported, and the returned
    fluxes conserve mass for the applied fraction."""
    print('\nC. robustness')
    ut = np.array(truth.predict(ut_prev))
    dt = truth.attr('dt')
    qstar = np.array(truth.fluxes({}, ut))
    p0 = np.array(face.predict(ut_prev, dt))
    tau = 1e4 * (qstar - np.array(face.fluxes({}, ut)))
    face.residual({}, ut_prev, ut)
    u, q = map(np.array, face.correct_face({}, p0, tau))
    info = face.attr('correction')
    div = np.array(face.divergence(tau))
    balance = np.array(face.residual({}, ut_prev, u)) + info['fraction'] * div
    bal = np.linalg.norm(balance) / np.linalg.norm(info['fraction'] * div) \
        if info['fraction'] > 0 else np.linalg.norm(balance)
    print(f"  fraction {info['fraction']}, {info['attempts']} attempts, "
          f"failed cells {[c for c, _, _ in info['failed_cells']]}, balance {bal:.2e}")
    truth.correct({}, ut, np.zeros(truth.ndof))
    ok = info['fraction'] < 1 and len(info['failed_cells']) > 0 and bal < 1e-6
    print(f'C. robustness: {"PASSED" if ok else "FAILED"}')
    return ok


def main(spe1):
    with tempfile.TemporaryDirectory() as tmp:
        make_decks.main(spe1, tmp)
        decks = {n: os.path.join(tmp, n) for n in ('COARSE', 'FINE', 'FINE_PERTURBED')}

        sims = [Remote('BlackOil', decks[n], ARGS)
                for n in ('FINE', 'FINE_PERTURBED', 'FINE_PERTURBED')]
        try:
            ok_a, ut = run('A. operator error', *sims, NSTEPS)
            ok_a &= robustness(sims[0], sims[1], ut)
        finally:
            for sim in sims:
                sim.close()

        sims = [Remote('BlackOil', decks[n], ARGS) for n in ('FINE', 'COARSE', 'COARSE')]
        try:
            r = make_decks.R
            restriction = Restriction(sims[0], sims[1], (r, r, r))
            ok_b, _ = run('B. grid error', *sims, NSTEPS, restriction, hybrid=True)
        finally:
            for sim in sims:
                sim.close()
    return 0 if ok_a and ok_b else 1


if __name__ == '__main__':
    raise SystemExit(main(os.path.abspath(sys.argv[1])))
