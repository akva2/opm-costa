"""Checks of face CoSTA corrections against a truth, shared by the tests.

A truth simulator runs alongside two models, one corrected with face CoSTA
and one with RHS CoSTA, using exact labels computed from the truth. The
truth can be the same grid as the models, or a nested finer grid, which is
then restricted to the model grid.
"""

import numpy as np

DAY = 86400.0
NEQ = 3  # equations per cell


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def flux_error(q, qstar):
    """Max relative error over components with nonzero reference flux."""
    q, qstar = q.reshape(-1, NEQ), qstar.reshape(-1, NEQ)
    return max(rel(q[:, i], qstar[:, i]) for i in range(NEQ)
               if np.linalg.norm(qstar[:, i]) > 0)


def state_errors(u, ustar, uprev, names):
    """Pressure error relative to the change over the step, saturations as
    max absolute error, others max relative error where the corresponding
    phase is present."""
    n = len(names)
    u, ustar, uprev = (x.reshape(-1, n) for x in (u, ustar, uprev))
    var = dict(zip(names, range(n)))
    sw = ustar[:, var['sw']] if 'sw' in var else 0.0
    sg = ustar[:, var['sg']] if 'sg' in var else 0.0
    # Rs/Rv are undefined without oil/gas, and then arbitrary
    present = {'rs': 1 - sw - sg > 1e-8, 'rv': sg > 1e-8}
    errors = {}
    for name, i in var.items():
        if name == 'p':
            errors[name] = np.linalg.norm(u[:, i] - ustar[:, i]) / \
                np.linalg.norm(ustar[:, i] - uprev[:, i])
        elif name.startswith('s'):
            errors[name] = np.abs(u[:, i] - ustar[:, i]).max()
        else:
            m = present.get(name, np.ones(len(u), dtype=bool))
            errors[name] = np.abs(u[m, i] - ustar[m, i]).max() / \
                max(np.abs(ustar[m, i]).max(), 1e-300) if m.any() else 0.0
    return errors


def gas_present(u, names):
    return u.reshape(-1, len(names))[:, names.index('sg')] > 0


class Restriction:
    """Transfer from a fine grid to a coarse grid it is nested in.

    Cells are matched by Cartesian indices, with the refinement ratio in
    each direction. States are restricted mass consistently: the coarse
    state holds the summed component amounts of its children. A coarse cell
    with inactive children gets its amounts scaled up, as if they had the
    state of the active ones. Fluxes are summed over the fine connections
    between the children of two coarse cells.
    """
    def __init__(self, fine, coarse, ratio):
        self.fine, self.coarse = fine, coarse
        self.names = coarse.state_vars
        cdims = np.array(coarse.attr('cartesian_dims'))
        fdims = np.array(fine.attr('cartesian_dims'))
        ratio = np.array(ratio)
        assert np.array_equal(cdims * ratio, fdims), (cdims, fdims, ratio)

        # coarse active index for each Cartesian index, -1 if inactive
        ccart = np.array(coarse.attr('cartesian_indices'))
        lookup = np.full(int(np.prod(cdims)), -1)
        lookup[ccart] = np.arange(len(ccart))

        def parent(cart):
            cart = np.asarray(cart)
            i, j, k = cart % fdims[0], (cart // fdims[0]) % fdims[1], \
                cart // (fdims[0] * fdims[1])
            return lookup[i // ratio[0] + cdims[0] * (j // ratio[1] +
                                                      cdims[1] * (k // ratio[2]))]

        self.ncoarse = len(ccart)
        self.parent = parent(fine.attr('cartesian_indices'))
        self.orphans = int((self.parent < 0).sum())
        children = np.bincount(self.parent[self.parent >= 0], minlength=self.ncoarse)
        if (children == 0).any():
            raise ValueError(f'{int((children == 0).sum())} coarse cells without children')
        self.scale = int(np.prod(ratio)) / children
        self.incomplete = int((children < np.prod(ratio)).sum())

        index = {tuple(c): n for n, c in enumerate(coarse.faces()['cartesian'])}
        self.ncoarse_faces = len(index)
        fine_faces = np.asarray(fine.faces()['cells'])
        pm, pp = self.parent[fine_faces[:, 0]], self.parent[fine_faces[:, 1]]
        self.face_map, self.unmatched = [], 0
        for f, (m, p) in enumerate(zip(pm, pp)):
            if m == p:
                continue
            if m < 0 or p < 0:
                self.unmatched += 1
                continue
            key = (int(ccart[m]), int(ccart[p]))
            if key in index:
                self.face_map.append((f, index[key], 1.0))
            elif key[::-1] in index:
                self.face_map.append((f, index[key[::-1]], -1.0))
            else:
                self.unmatched += 1
        self.face_map = np.array(self.face_map)

    def summary(self):
        return (f'{self.orphans} fine cells without active parent, '
                f'{self.incomplete} coarse cells with inactive children, '
                f'{self.unmatched} fine connections without coarse connection')

    def sum(self, x, n):
        """Sum n values per fine cell over each coarse cell."""
        x = x.reshape(-1, n)
        m = self.parent >= 0
        return np.column_stack([np.bincount(self.parent[m], x[m, i], self.ncoarse)
                                for i in range(n)])

    def state(self, u):
        """The coarse state holding the fine component amounts, starting
        from cell averages."""
        n = len(self.names)
        m = self.parent >= 0
        guess = self.sum(u, n) / np.bincount(self.parent[m], minlength=self.ncoarse)[:, None]
        amounts = self.sum(np.array(self.fine.storage({}, u)), NEQ) * self.scale[:, None]
        return np.array(self.coarse.state_from_storage({}, amounts.ravel(), guess.ravel()))

    def flux(self, q):
        q = q.reshape(-1, NEQ)
        f, c, sign = self.face_map[:, 0].astype(int), self.face_map[:, 1].astype(int), \
            self.face_map[:, 2]
        result = np.zeros((self.ncoarse_faces, NEQ))
        np.add.at(result, c, q[f] * sign[:, None])
        return result.ravel()


def run(title, truth, face, rhs, nsteps, restriction=None, hybrid=False,
        flux_tol=1e-5, require_switching=True, max_chops=5):
    """Correct nsteps steps of the models with exact labels from the truth.

    The truth chooses the step lengths, and the models use the same. If the
    face correction fails, the step is predicted again shorter. With hybrid,
    face CoSTA also gets the cell correction sigma* = delta* + D tau* for
    what face fluxes cannot represent, such as differing well models.
    Returns whether the checks passed, and the last truth state.
    """
    restrict_state = restriction.state if restriction else (lambda u: u)
    restrict_flux = restriction.flux if restriction else (lambda q: q)
    names = face.state_vars
    print(f'\n{title}: state {names}')
    if restriction:
        print(f'  restriction: {restriction.summary()}')
    ut_prev = np.array(truth.initial_condition({}))
    ustar_prev = restrict_state(ut_prev)
    ok, switched_total = True, 0
    for step in range(nsteps):
        dt = 0.0
        for chop in range(max_chops):
            try:
                ut = np.array(truth.predict(ut_prev, dt))
            except RuntimeError as e:
                print(f'step {step}: truth failed: {e}')
                dt /= 3
                continue
            dt, time = truth.attr('dt'), truth.attr('time')
            ustar = restrict_state(ut)
            qstar = restrict_flux(np.array(truth.fluxes({}, ut)))

            p0 = np.array(face.predict(ustar_prev, dt))
            q0 = np.array(face.fluxes({}, p0))

            # Labels
            tau = qstar - np.array(face.fluxes({}, ustar))
            delta = np.array(face.residual({}, ustar_prev, ustar))
            div_tau = np.array(face.divergence(tau))

            # Cell correction for what face fluxes cannot represent
            sigma = delta + div_tau if hybrid else np.zeros(0)
            try:
                u_face, q_face = map(np.array,
                                     face.correct_face({}, p0, tau, sigma, strict=True))
                break
            except RuntimeError as e:
                print(f'step {step}: face CoSTA failed for dt {dt / DAY:.3g} days: '
                      f'{str(e)[:120]}')
                dt /= 3
        else:
            print(f'{title}: FAILED, could not correct step {step}')
            return False, ut_prev
        balance = np.array(face.residual({}, ustar_prev, u_face)) + div_tau
        if hybrid:
            balance -= sigma
        results = [('uncorrected', p0, q0), ('face CoSTA', u_face, q_face)]

        p0_rhs = np.array(rhs.predict(ustar_prev, dt))
        rhs.residual({}, ustar_prev, ustar)  # same old state as face
        u_rhs = np.array(rhs.correct({}, p0_rhs, delta))
        if rhs.attr('correction')['fraction'] < 1:
            print(f"step {step}: RHS CoSTA applied {rhs.attr('correction')['fraction']}")
        results.append(('RHS CoSTA', u_rhs, np.array(rhs.fluxes({}, u_rhs))))

        # Steps are finished by a correction, a zero one for the truth
        truth.correct({}, ut, np.zeros(truth.ndof))

        print(f'step {step} at {time / DAY:.4g} days, dt {dt / DAY:.3g} days: '
              f'|delta* + D tau*| / |delta*| = {rel(-div_tau, delta):.2e}')
        if hybrid:
            mag = np.abs(sigma.reshape(-1, NEQ)).sum(axis=1)
            cells = np.argsort(mag)[::-1]
            share = np.cumsum(mag[cells]) / mag.sum()
            n99 = int(np.searchsorted(share, 0.99)) + 1
            print(f'  sigma*: 99% in {n99} of {len(mag)} cells, largest at cells '
                  f'{cells[:3].tolist()}')
        errors = {}
        for name, u, q in results:
            errs = state_errors(u, ustar, ustar_prev, names)
            errors[name] = flux_error(q, qstar)
            print(f'  {name:11s}: ' + '  '.join(f'{k} {v:.2e}' for k, v in errs.items())
                  + f'  q err {errors[name]:.2e}')
        bal = np.linalg.norm(balance) / np.linalg.norm(div_tau)
        msg = f'  face CoSTA balance {bal:.2e}'
        if 'sg' in names:
            switched = int((gas_present(u_face, names) != gas_present(p0, names)).sum())
            over_step = int((gas_present(u_face, names) !=
                             gas_present(ustar_prev, names)).sum())
            switched_total += switched
            msg += f', switched {switched} cells, {over_step} over the step'
        print(msg)
        ok &= bal < 1e-6 and errors['face CoSTA'] < flux_tol and \
            errors['RHS CoSTA'] > 10 * flux_tol
        ut_prev, ustar_prev = ut, ustar
    if require_switching:
        ok &= switched_total > 0
    print(f'{title}: {"PASSED" if ok else "FAILED"}')
    return ok, ut_prev
