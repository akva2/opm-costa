"""Check the face CoSTA API on a nested single-phase fine/coarse pair.

Uses exact labels from the fine simulation (teacher forcing), so we can
verify the properties from the face CoSTA description:

 * The labels are consistent: delta* = -D tau* (up to accumulation mismatch).
 * correct_face() with tau* recovers both p* and Q*.
 * correct_face() conserves mass for any tau: D q = b - accumulation.
 * correct() (RHS CoSTA) recovers p*, but not Q*.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from costa_remote import Remote  # noqa: E402

R = 4                  # coarsening ratio in both directions
NXC, NYC = 8, 4        # coarse grid
NXF, NYF = NXC * R, NYC * R
NSTEPS = 3
P_INIT = 200e5          # initial PRESSURE in deck

# CoSTA assumes accurate coarse solves, far below flow's default tolerances.
# Only the states in memory are used, so no ECL output files are written.
ARGS = ['--enable-ecl-output=false', '--linear-solver-reduction=1e-12',
        '--tolerance-cnv=1e-10', '--tolerance-cnv-relaxed=1e-10',
        '--tolerance-mb=1e-14', '--tolerance-mb-relaxed=1e-14',
        '--newton-max-iterations=20']


def restrict_pressure(pf):
    """Average fine pressures over each coarse cell."""
    return pf.reshape(NYC, R, NXC, R).mean(axis=(1, 3)).ravel()


def restrict_flux(fine_faces, coarse_faces, qf):
    """Sum fine face fluxes over each coarse face."""
    def coarse_cell(c):
        i, j = c % NXF, c // NXF
        return (j // R) * NXC + i // R
    index = {tuple(c): n for n, c in enumerate(coarse_faces['cartesian'])}
    q = np.zeros(len(index))
    for (m, p), flux in zip(fine_faces['cartesian'], qf):
        key = (coarse_cell(m), coarse_cell(p))
        if key[0] != key[1]:
            q[index[key]] += flux
    return q


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def perr(p, pstar):
    """Pressure error relative to the pressure change from the initial state."""
    return np.linalg.norm(p - pstar) / np.linalg.norm(pstar - P_INIT)


def main():
    fine = Remote('OnePhase', 'FINE', ARGS)
    face = Remote('OnePhase', 'COARSE', ARGS)    # face CoSTA instance
    rhs = Remote('OnePhase', 'COARSE', ARGS)     # RHS CoSTA instance
    ffaces, cfaces = fine.faces(), face.faces()
    print(f'fine: {fine.ndof} cells, {fine.nfaces} faces; '
          f'coarse: {face.ndof} cells, {face.nfaces} faces')

    uf = np.zeros(fine.ndof)
    pstar_prev = np.array(face.initial_condition({}))
    ok = True
    for step in range(NSTEPS):
        # Steps are finished by a correction, a zero one for the fine grid
        uf = np.array(fine.correct({}, fine.predict(uf), np.zeros(fine.ndof)))
        qf = np.array(fine.fluxes({}, uf))
        pstar = restrict_pressure(uf)
        qstar = restrict_flux(ffaces, cfaces, qf)

        p0 = np.array(face.predict(pstar_prev))
        p0_rhs = np.array(rhs.predict(pstar_prev))
        q0 = np.array(face.fluxes({}, p0))

        # Exact labels
        tau = qstar - np.array(face.fluxes({}, pstar))
        delta = np.array(face.residual({}, pstar_prev, pstar))
        div_tau = np.array(face.divergence(tau))
        rhs.residual({}, pstar_prev, pstar)  # keep instances in the same state

        # Face CoSTA
        p_face, q_face = map(np.array, face.correct_face({}, p0, tau))
        # Conservation: accumulation + D q - b = R(p) + D tau
        balance = np.array(face.residual({}, pstar_prev, p_face)) + div_tau

        # RHS CoSTA
        p_rhs = np.array(rhs.correct({}, p0_rhs, delta))
        q_rhs = np.array(rhs.fluxes({}, p_rhs))

        print(f'step {step}:')
        print(f'  |delta* + D tau*| / |delta*|   = {rel(-div_tau, delta):.3e}')
        print(f'  uncorrected : p err {perr(p0, pstar):.3e}'
              f'  q err {rel(q0, qstar):.3e}')
        print(f'  RHS CoSTA   : p err {perr(p_rhs, pstar):.3e}'
              f'  q err {rel(q_rhs, qstar):.3e}')
        print(f'  face CoSTA  : p err {perr(p_face, pstar):.3e}'
              f'  q err {rel(q_face, qstar):.3e}'
              f'  balance {np.linalg.norm(balance) / np.linalg.norm(div_tau):.3e}')
        ok &= rel(q_face, qstar) < 1e-6 and rel(q_rhs, qstar) > 1e-3
        pstar_prev = pstar

    for sim in (fine, face, rhs):
        sim.close()
    print('PASSED' if ok else 'FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main())
