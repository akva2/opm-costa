"""Run opm_CoSTA simulators in separate processes.

OPM keeps global state, so only one simulator can live in a process.
"""

import multiprocessing as mp


class Remote:
    """Proxy for an opm_CoSTA simulator running in its own process.

    Method calls are forwarded to the process. ndof, nfaces and state_vars
    are read once.
    """
    def __init__(self, variant, deck, args=(), timeout=1800):
        self._timeout = timeout
        self._conn, child = mp.Pipe()
        self._proc = mp.get_context('spawn').Process(target=_serve,
                                                     args=(child, variant, deck,
                                                           list(args)))
        self._proc.start()
        self.ndof, self.nfaces, self.state_vars = self._recv()

    def _recv(self):
        if not self._conn.poll(self._timeout):
            self._proc.kill()
            raise RuntimeError('simulator process did not respond')
        result = self._conn.recv()
        if isinstance(result, Exception):
            raise result
        return result

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


def _serve(conn, variant, deck, args):
    try:
        import opm_CoSTA
        sim = getattr(opm_CoSTA, variant)(deck, False, args)
        conn.send((sim.ndof, sim.nfaces, sim.state_vars))
    except Exception as e:
        conn.send(e)
        return
    while (msg := conn.recv()) is not None:
        name, args = msg
        try:
            conn.send(getattr(sim, name)(*args))
        except Exception as e:
            conn.send(e)
