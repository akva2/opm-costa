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
NSTEPS = 8
MAX_CHOPS = 5
DAY = 86400.0

# Transmissibility perturbation of the model, added to the EDIT section
PERTURBATION = '''MULTIPLY
 TRANX 0.5 /
 TRANY 2.0 /
/
'''

# CoSTA assumes accurate solves, far below flow's default tolerances.
# The CNV criterion scales with dt, so short substeps need it tight.
# Only the states in memory are used, so no ECL output files are written.
ARGS = ['--enable-ecl-output=false', '--linear-solver-reduction=1e-12',
        '--tolerance-cnv=1e-10', '--tolerance-cnv-relaxed=1e-10',
        '--tolerance-mb=1e-14', '--tolerance-mb-relaxed=1e-14',
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
        # The truth chooses the step length, and the models use the same.
        # If the face correction fails, the step is predicted again shorter.
        dt = 0.0
        for chop in range(MAX_CHOPS):
            try:
                ustar = np.array(truth.predict(ustar_prev, dt))
            except RuntimeError as e:
                print(f'step {step}: truth failed: {e}')
                dt /= 3
                continue
            dt, time = truth.attr('dt'), truth.attr('time')
            qstar = np.array(truth.fluxes({}, ustar))

            p0 = np.array(face.predict(ustar_prev, dt))
            q0 = np.array(face.fluxes({}, p0))

            # Exact labels
            tau = qstar - np.array(face.fluxes({}, ustar))
            delta = np.array(face.residual({}, ustar_prev, ustar))
            div_tau = np.array(face.divergence(tau))
            try:
                u_face, q_face = map(np.array, face.correct_face({}, p0, tau, strict=True))
                break
            except RuntimeError as e:
                print(f'step {step}: face CoSTA failed for dt {dt / DAY:.3g} days: {e}')
                dt /= 3
        else:
            print(f'FAILED, could not correct step {step}')
            return 1
        balance = np.array(face.residual({}, ustar_prev, u_face)) + div_tau

        # RHS CoSTA, same step
        p0_rhs = np.array(rhs.predict(ustar_prev, dt))
        rhs.residual({}, ustar_prev, ustar)  # same old state as face
        try:
            u_rhs = np.array(rhs.correct({}, p0_rhs, delta))
            if rhs.attr('correction')['fraction'] < 1:
                print(f"step {step}: RHS CoSTA applied {rhs.attr('correction')['fraction']}")
            q_rhs = np.array(rhs.fluxes({}, u_rhs))
        except RuntimeError as e:
            print(f'step {step}: RHS CoSTA failed: {e}')
            rhs.correct({}, p0_rhs, np.zeros(rhs.ndof))  # finish the step
            u_rhs = None

        # Steps are finished by a correction, a zero one for the truth
        truth.correct({}, ustar, np.zeros(truth.ndof))

        print(f'step {step} at {time / DAY:.4g} days, dt {dt / DAY:.3g} days:')
        print(f'  |delta* + D tau*| / |delta*| = {rel(-div_tau, delta):.3e}')
        results = [('uncorrected', p0, q0), ('face CoSTA', u_face, q_face)]
        if u_rhs is not None:
            results.append(('RHS CoSTA', u_rhs, q_rhs))
        for name, u, q in results:
            errs = state_errors(u, ustar, ustar_prev, face.state_vars)
            print(f'  {name:11s}: ' + '  '.join(f'{k} {v:.2e}' for k, v in errs.items())
                  + f'  q err {flux_error(q, qstar):.3e}')
        bal = np.linalg.norm(balance) / np.linalg.norm(div_tau)
        print(f'  face balance |R + D tau| / |D tau| = {bal:.3e}')
        # Where a phase potential vanishes across a face, the upwind direction
        # is ambiguous, and truth and corrected model can pick different
        # neighbouring solutions. This gives flux errors up to about 1e-3.
        ok &= u_rhs is not None and bal < 1e-6 and \
            flux_error(q_face, qstar) < 1e-3 and flux_error(q_rhs, qstar) > 1e-2
        ustar_prev = ustar

    print('PASSED' if ok else 'FAILED')
    return 0 if ok else 1


if __name__ == '__main__':
    raise SystemExit(main(os.path.abspath(sys.argv[1])))
