"""Check the face CoSTA API on synthetic blackoil decks based on SPE1.

Two scenarios, with decks from make_decks.py:

 A. Operator error: truth FINE, model FINE_PERTURBED (perturbed transmissibilities).
    The exact labels satisfy delta* = -D tau*, so face CoSTA must recover the
    truth state and fluxes, while RHS CoSTA only recovers the state.
 B. Grid error: truth FINE restricted to the COARSE grid. The labels are only
    approximately consistent (nonlinear storage, producer well), so face CoSTA
    must conserve mass and improve the fluxes over both alternatives.

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
from costa_remote import Remote  # noqa: E402
import make_decks  # noqa: E402

NSTEPS = 6
NEQ = 3  # equations per cell

# CoSTA assumes accurate solves, far below flow's default tolerances.
ARGS = ['--linear-solver-reduction=1e-10',
        '--tolerance-cnv=1e-8', '--tolerance-cnv-relaxed=1e-8',
        '--tolerance-mb=1e-12', '--tolerance-mb-relaxed=1e-12',
        '--newton-max-iterations=30']


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def flux_error(q, qstar):
    """Max relative error over components with nonzero reference flux."""
    q, qstar = q.reshape(-1, NEQ), qstar.reshape(-1, NEQ)
    return max(rel(q[:, i], qstar[:, i]) for i in range(NEQ)
               if np.linalg.norm(qstar[:, i]) > 0)


def state_errors(u, ustar, uprev, names):
    """Pressure error relative to the change over the step, saturations as
    max absolute error, Rs max relative error where oil is present."""
    n = len(names)
    u, ustar, uprev = (x.reshape(-1, n) for x in (u, ustar, uprev))
    var = dict(zip(names, range(n)))
    so = 1 - ustar[:, var['sw']] - ustar[:, var['sg']]
    errors = {}
    for name, i in var.items():
        if name == 'p':
            errors[name] = rel(u[:, i], ustar[:, i]) * np.linalg.norm(ustar[:, i]) / \
                np.linalg.norm(ustar[:, i] - uprev[:, i])
        elif name.startswith('s'):
            errors[name] = np.abs(u[:, i] - ustar[:, i]).max()
        else:
            m = so > 1e-8
            errors[name] = np.abs(u[m, i] - ustar[m, i]).max() / np.abs(ustar[m, i]).max()
    return errors


def gas_present(u, names):
    return u.reshape(-1, len(names))[:, names.index('sg')] > 0


class Restriction:
    """Transfer from the fine to the coarse grid, both Cartesian and all active."""
    def __init__(self, r, names, fine_faces, coarse_faces):
        self.r, self.names = r, names
        nx, ny, nz = make_decks.NX, make_decks.NY, make_decks.NZ
        index = {tuple(c): n for n, c in enumerate(coarse_faces['cartesian'])}

        def coarse(c):
            fx, fy = nx * r, ny * r
            i, j, k = c % fx, (c // fx) % fy, c // (fx * fy)
            return i // r + nx * (j // r + ny * (k // r))
        self.parent = np.array([coarse(c) for c in range(nx * ny * nz * r**3)])
        self.face_map = []  # (fine face, coarse face)
        for f, (m, p) in enumerate(fine_faces['cartesian']):
            cm, cp = coarse(m), coarse(p)
            if cm != cp:
                self.face_map.append((f, index[(cm, cp)]))
        self.ncoarse = nx * ny * nz
        self.ncoarse_faces = len(index)

    def state(self, u):
        n = len(self.names)
        u = u.reshape(-1, n)
        var = dict(zip(self.names, range(n)))
        result = np.zeros((self.ncoarse, n))
        count = np.bincount(self.parent, minlength=self.ncoarse)
        for name, i in var.items():
            if name == 'rs':
                continue
            result[:, i] = np.bincount(self.parent, u[:, i], self.ncoarse) / count
        # Rs weighted by oil saturation (all cells have equal pore volume)
        so = 1 - u[:, var['sw']] - u[:, var['sg']]
        wsum = np.bincount(self.parent, so, self.ncoarse)
        rs = np.bincount(self.parent, so * u[:, var['rs']], self.ncoarse)
        mean = np.bincount(self.parent, u[:, var['rs']], self.ncoarse) / count
        result[:, var['rs']] = np.where(wsum > 0, rs / np.maximum(wsum, 1e-300), mean)
        return result.ravel()

    def flux(self, q):
        q = q.reshape(-1, NEQ)
        result = np.zeros((self.ncoarse_faces, NEQ))
        for f, c in self.face_map:
            result[c] += q[f]
        return result.ravel()


def run(title, truth, face, rhs, restrict_state, restrict_flux, exact):
    names = face.state_vars
    print(f'\n{title}: state {names}')
    ut_prev = np.array(truth.initial_condition({}))
    ustar_prev = restrict_state(ut_prev)
    ok, switched_total = True, 0
    for step in range(NSTEPS):
        # Steps are finished by a correction, a zero one for the truth
        ut = np.array(truth.correct({}, truth.predict(ut_prev), np.zeros(truth.ndof)))
        ustar = restrict_state(ut)
        qstar = restrict_flux(np.array(truth.fluxes({}, ut)))

        p0 = np.array(face.predict(ustar_prev))
        p0_rhs = np.array(rhs.predict(ustar_prev))
        q0 = np.array(face.fluxes({}, p0))

        # Labels
        tau = qstar - np.array(face.fluxes({}, ustar))
        delta = np.array(face.residual({}, ustar_prev, ustar))
        div_tau = np.array(face.divergence(tau))
        rhs.residual({}, ustar_prev, ustar)  # same old state as face

        results = [('uncorrected', p0, q0)]
        try:
            u_face, q_face = map(np.array, face.correct_face({}, p0, tau))
            balance = np.array(face.residual({}, ustar_prev, u_face)) + div_tau
            results.append(('face CoSTA', u_face, q_face))
        except RuntimeError as e:
            print(f'step {step}: face CoSTA failed: {e}')
            face.correct({}, p0, np.zeros(face.ndof))  # finish the step
            u_face = None
        try:
            u_rhs = np.array(rhs.correct({}, p0_rhs, delta))
            results.append(('RHS CoSTA', u_rhs, np.array(rhs.fluxes({}, u_rhs))))
        except RuntimeError as e:
            print(f'step {step}: RHS CoSTA failed: {e}')
            rhs.correct({}, p0_rhs, np.zeros(rhs.ndof))  # finish the step
            u_rhs = None

        print(f'step {step}: |delta* + D tau*| / |delta*| = {rel(-div_tau, delta):.2e}')
        errors = {}
        for name, u, q in results:
            errs = state_errors(u, ustar, ustar_prev, names)
            errors[name] = flux_error(q, qstar)
            print(f'  {name:11s}: ' + '  '.join(f'{k} {v:.2e}' for k, v in errs.items())
                  + f'  q err {errors[name]:.2e}')
        if u_face is None or u_rhs is None:
            ok = False
            continue
        switched = int((gas_present(u_face, names) != gas_present(p0, names)).sum())
        over_step = int((gas_present(u_face, names) !=
                         gas_present(ustar_prev, names)).sum())
        switched_total += switched
        bal = np.linalg.norm(balance) / np.linalg.norm(div_tau)
        print(f'  face CoSTA balance {bal:.2e}, switched {switched} cells, '
              f'{over_step} over the step')
        if exact:
            ok &= errors['face CoSTA'] < 1e-5 and errors['RHS CoSTA'] > 1e-3
        else:
            ok &= bal < 1e-6 and errors['face CoSTA'] < min(errors['uncorrected'],
                                                             errors['RHS CoSTA'])
        ut_prev, ustar_prev = ut, ustar
    if exact:
        ok &= switched_total > 0
    print(f'{title}: {"PASSED" if ok else "FAILED"}')
    return ok


def main(spe1):
    with tempfile.TemporaryDirectory() as tmp:
        make_decks.main(spe1, tmp)
        decks = {n: os.path.join(tmp, n) for n in ('COARSE', 'FINE', 'FINE_PERTURBED')}

        sims = [Remote('BlackOil', decks[n], ARGS)
                for n in ('FINE', 'FINE_PERTURBED', 'FINE_PERTURBED')]
        try:
            ok_a = run('A. operator error', *sims, lambda u: u, lambda q: q, exact=True)
        finally:
            for sim in sims:
                sim.close()

        sims = [Remote('BlackOil', decks[n], ARGS) for n in ('FINE', 'COARSE', 'COARSE')]
        try:
            restriction = Restriction(make_decks.R, sims[1].state_vars,
                                      sims[0].faces(), sims[1].faces())
            ok_b = run('B. grid error', *sims, restriction.state, restriction.flux,
                       exact=False)
        finally:
            for sim in sims:
                sim.close()
    return 0 if ok_a and ok_b else 1


if __name__ == '__main__':
    raise SystemExit(main(os.path.abspath(sys.argv[1])))
