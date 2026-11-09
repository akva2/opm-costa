"""Check the face CoSTA API on the Norne blackoil case.

There is no nested coarse Norne grid, so the model error is introduced
in the operator instead: the "truth" is the original deck, the model has
perturbed inter-cell transmissibilities. Wells, PVT and pore volumes are
identical, so the exact labels satisfy delta* = -D tau*, and

 * correct_face() with tau* recovers both the truth state and fluxes,
 * correct_face() conserves mass: accumulation + D q - b = 0,
 * correct() (RHS CoSTA) recovers the truth state, but not its fluxes.

Usage: python test_norne.py <directory with NORNE_ATW2013.DATA>
"""

import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from costa_remote import Remote  # noqa: E402

DECK = 'NORNE_ATW2013.DATA'
NSTEPS = 5

# Transmissibility perturbation of the model, added to the EDIT section
PERTURBATION = '''MULTIPLY
 TRANX 0.5 /
 TRANY 2.0 /
/
'''

# CoSTA assumes accurate solves, far below flow's default tolerances.
ARGS = ['--linear-solver-reduction=1e-10',
        '--tolerance-cnv=1e-8', '--tolerance-cnv-relaxed=1e-8',
        '--tolerance-mb=1e-12', '--tolerance-mb-relaxed=1e-12',
        '--newton-max-iterations=30']

NEQ = 3  # equations per cell


def make_decks(norne, tmp):
    """Write the truth and model decks to tmp, sharing norne/INCLUDE."""
    os.symlink(os.path.join(norne, 'INCLUDE'), os.path.join(tmp, 'INCLUDE'))
    with open(os.path.join(norne, DECK)) as f:
        lines = f.read().split('\n')
    with open(os.path.join(tmp, 'TRUTH.DATA'), 'w') as f:
        f.write('\n'.join(lines))
    edit = next(i for i, line in enumerate(lines) if line.strip() == 'EDIT')
    lines.insert(edit + 1, PERTURBATION)
    with open(os.path.join(tmp, 'MODEL.DATA'), 'w') as f:
        f.write('\n'.join(lines))
    return os.path.join(tmp, 'TRUTH'), os.path.join(tmp, 'MODEL')


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def state_errors(u, ustar, uprev, names):
    """Errors per state variable: pressure relative to the change over the
    step, saturations as max absolute error, others max relative error
    where the corresponding phase is present."""
    n = len(names)
    u, ustar, uprev = (x.reshape(-1, n) for x in (u, ustar, uprev))
    errors = {}
    sat = dict(zip(names, ustar.T))
    so = 1 - sat['sw'] - sat['sg']
    # Rs/Rv are undefined without oil/gas, and then arbitrary
    present = {'rs': so > 1e-8, 'rv': sat['sg'] > 1e-8}
    for i, name in enumerate(names):
        if name in present:
            m = present[name]
            errors[name] = np.abs(u[m, i] - ustar[m, i]).max() / np.abs(ustar[m, i]).max()
        elif name == 'p':
            errors[name] = np.linalg.norm(u[:, i] - ustar[:, i]) / \
                np.linalg.norm(ustar[:, i] - uprev[:, i])
        elif name.startswith('s'):
            errors[name] = np.abs(u[:, i] - ustar[:, i]).max()
        else:
            errors[name] = np.abs(u[:, i] - ustar[:, i]).max() / np.abs(ustar[:, i]).max()
    return errors


def flux_error(q, qstar):
    """Max relative error over components with nonzero reference flux."""
    q, qstar = q.reshape(-1, NEQ), qstar.reshape(-1, NEQ)
    return max(rel(q[:, i], qstar[:, i]) for i in range(NEQ)
               if np.linalg.norm(qstar[:, i]) > 0)


def main(norne):
    with tempfile.TemporaryDirectory() as tmp:
        truth_deck, model_deck = make_decks(norne, tmp)
        truth = Remote('BlackOil', truth_deck, ARGS)
        face = Remote('BlackOil', model_deck, ARGS)   # face CoSTA instance
        rhs = Remote('BlackOil', model_deck, ARGS)    # RHS CoSTA instance
        try:
            return run(truth, face, rhs)
        finally:
            for sim in (truth, face, rhs):
                sim.close()


def run(truth, face, rhs):
    tf, mf = truth.faces(), face.faces()
    assert np.array_equal(tf['cells'], mf['cells'])
    trans_ratio = np.asarray(mf['trans']) / np.asarray(tf['trans'])
    print(f'{face.ndof // NEQ} cells, {face.nfaces} faces, state {face.state_vars}, model/truth '
          f'transmissibility ratios {np.unique(np.round(trans_ratio, 6))}')

    ustar_prev = np.array(truth.initial_condition({}))
    ok = True
    for step in range(NSTEPS):
        # Steps are finished by a correction, a zero one for the truth
        ustar = np.array(truth.correct({}, truth.predict(ustar_prev),
                                       np.zeros(truth.ndof)))
        qstar = np.array(truth.fluxes({}, ustar))

        p0 = np.array(face.predict(ustar_prev))
        p0_rhs = np.array(rhs.predict(ustar_prev))
        q0 = np.array(face.fluxes({}, p0))

        # Exact labels
        tau = qstar - np.array(face.fluxes({}, ustar))
        delta = np.array(face.residual({}, ustar_prev, ustar))
        div_tau = np.array(face.divergence(tau))
        rhs.residual({}, ustar_prev, ustar)  # same old state as face

        # Face CoSTA
        try:
            u_face, q_face = map(np.array, face.correct_face({}, p0, tau))
            balance = np.array(face.residual({}, ustar_prev, u_face)) + div_tau
        except RuntimeError as e:
            print(f'step {step}: face CoSTA failed: {e}')
            face.correct({}, p0, np.zeros(face.ndof))  # finish the step
            u_face = None

        # RHS CoSTA
        try:
            u_rhs = np.array(rhs.correct({}, p0_rhs, delta))
            q_rhs = np.array(rhs.fluxes({}, u_rhs))
        except RuntimeError as e:
            print(f'step {step}: RHS CoSTA failed: {e}')
            rhs.correct({}, p0_rhs, np.zeros(rhs.ndof))  # finish the step
            u_rhs = None

        print(f'step {step}:')
        print(f'  |delta* + D tau*| / |delta*| = {rel(-div_tau, delta):.3e}')
        results = [('uncorrected', p0, q0)]
        if u_rhs is not None:
            results.append(('RHS CoSTA', u_rhs, q_rhs))
        if u_face is not None:
            results.append(('face CoSTA', u_face, q_face))
        for name, u, q in results:
            errs = state_errors(u, ustar, ustar_prev, face.state_vars)
            print(f'  {name:11s}: ' + '  '.join(f'{k} {v:.2e}' for k, v in errs.items())
                  + f'  q err {flux_error(q, qstar):.3e}')
        if u_face is not None:
            print(f'  face balance |R + D tau| / |D tau| = '
                  f'{np.linalg.norm(balance) / np.linalg.norm(div_tau):.3e}')
        ok &= u_face is not None and u_rhs is not None and \
            flux_error(q_face, qstar) < 1e-5 and flux_error(q_rhs, qstar) > 1e-3
        ustar_prev = ustar

    print('PASSED' if ok else 'FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main(os.path.abspath(sys.argv[1])))
