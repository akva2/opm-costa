"""Check the face CoSTA API on a nested single-phase fine/coarse pair.

Uses exact labels from the fine simulation (teacher forcing), so we can
verify the properties from the face CoSTA description:

 * The labels are consistent: delta* = -D tau* (up to accumulation mismatch).
 * correct_face() with tau* recovers both p* and Q*.
 * correct_face() conserves mass for any tau: D q = b - accumulation.
 * correct() (RHS CoSTA) recovers p*, but not Q*.
"""

import multiprocessing as mp

import numpy as np

R = 4                  # coarsening ratio in both directions
NXC, NYC = 8, 4        # coarse grid
NXF, NYF = NXC * R, NYC * R
NSTEPS = 3
P_INIT = 200e5          # initial PRESSURE in deck

# CoSTA assumes accurate coarse solves, far below flow's default tolerances.
ARGS = ['--linear-solver-reduction=1e-12',
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


class Remote:
    """Run a simulator in a separate process.

    OPM keeps global state, so only one simulator can live in a process.
    Method calls are forwarded to the process.
    """
    def __init__(self, deck):
        self._conn, child = mp.Pipe()
        self._proc = mp.get_context('spawn').Process(target=Remote._serve,
                                                     args=(child, deck))
        self._proc.start()
        self.ndof, self.nfaces = self._recv()

    def _recv(self, timeout=300):
        if not self._conn.poll(timeout):
            self._proc.kill()
            raise RuntimeError('simulator process did not respond')
        return self._conn.recv()

    @staticmethod
    def _serve(conn, deck):
        import opm_CoSTA
        sim = opm_CoSTA.OnePhase(deck, False, ARGS)
        conn.send((sim.ndof, sim.nfaces))
        while (msg := conn.recv()) is not None:
            name, args = msg
            conn.send(getattr(sim, name)(*args))

    def __getattr__(self, name):
        if name.startswith('_'):
            raise AttributeError(name)
        def call(*args):
            self._conn.send((name, args))
            return self._recv()
        return call

    def close(self):
        self._conn.send(None)
        self._proc.join()


def rel(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)


def perr(p, pstar):
    """Pressure error relative to the pressure change from the initial state."""
    return np.linalg.norm(p - pstar) / np.linalg.norm(pstar - P_INIT)


def main():
    fine = Remote('FINE')
    face = Remote('COARSE')    # face CoSTA instance
    rhs = Remote('COARSE')     # RHS CoSTA instance
    ffaces, cfaces = fine.faces(), face.faces()
    print(f'fine: {fine.ndof} cells, {fine.nfaces} faces; '
          f'coarse: {face.ndof} cells, {face.nfaces} faces')

    uf = np.zeros(fine.ndof)
    pstar_prev = np.full(face.ndof, P_INIT)
    ok = True
    for step in range(NSTEPS):
        uf = np.array(fine.predict(uf))
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
