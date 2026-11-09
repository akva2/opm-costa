"""Generate nested blackoil test decks from SPE1 CASE2.

Writes to an output directory:
 * COARSE.DATA    - the SPE1 grid (10x10x3), gas injector replaced by a source
 * PERTURBED.DATA - COARSE with perturbed inter-cell transmissibilities
 * FINE.DATA      - COARSE refined R x R x R, with the source split over children
 * FINE_PERTURBED.DATA - FINE with perturbed inter-cell transmissibilities

Gas injection into undersaturated oil (no DRSDT) makes free gas appear and
redissolve, so cells switch primary variables between Rs and Sg.

Usage: python make_decks.py <path to SPE1CASE2.DATA> <output directory>
"""

import os
import sys

R = 2                                  # refinement factor in each direction
NX, NY, NZ = 10, 10, 3                 # coarse grid
DXY = 1000.0                           # ft
DZ = [20.0, 30.0, 50.0]                # ft
PERM = [500.0, 50.0, 200.0]            # mD
TOP = 8325.0                           # ft
INJ_CELL = (1, 1, 1)                   # coarse cell with gas source (1-based)
GAS_RATE = 100000e3 * 0.0533           # 100 MMscf/d as lb/day (DENSITY gas)
TSTEP = '30*10'                        # days

PERTURBATION = '''MULTIPLY
 TRANX 0.25 /
 TRANY 4.0 /
 TRANZ 0.1 /
/
'''


def sections(deck):
    """Split a deck into sections, keyed by section name."""
    names = ('RUNSPEC', 'GRID', 'EDIT', 'PROPS', 'REGIONS', 'SOLUTION',
             'SUMMARY', 'SCHEDULE')
    result, current = {}, None
    for line in deck.split('\n'):
        word = line.split('--')[0].strip()
        if word in names:
            current = word
            result[current] = []
        elif current:
            result[current].append(line)
    return {k: '\n'.join(v) for k, v in result.items()}


def grid(r):
    nx, ny, nz = NX * r, NY * r, NZ * r
    n = nx * ny
    dz = ' '.join(f'{n}*{DZ[k // r] / r:g}' for k in range(nz))
    perm = ' '.join(f'{n}*{PERM[k // r]:g}' for k in range(nz))
    return f'''GRID
INIT
DX
{n * nz}*{DXY / r:g} /
DY
{n * nz}*{DXY / r:g} /
DZ
{dz} /
TOPS
{n}*{TOP:g} /
PORO
{n * nz}*0.3 /
PERMX
{perm} /
PERMY
{perm} /
PERMZ
{perm} /
'''


def schedule(r):
    i, j, k = INJ_CELL
    children = [(ii, jj, kk)
                for kk in range((k - 1) * r + 1, k * r + 1)
                for jj in range((j - 1) * r + 1, j * r + 1)
                for ii in range((i - 1) * r + 1, i * r + 1)]
    source = '\n'.join(f'{ii} {jj} {kk} GAS {GAS_RATE / len(children):.10g} /'
                       for ii, jj, kk in children)
    p = NX * r, NY * r, NZ * r
    return f'''SCHEDULE
WELSPECS
 'PROD' 'G1' {p[0]} {p[1]} 8400 'OIL' /
/
COMPDAT
 'PROD' {p[0]} {p[1]} {p[2]} {p[2]} 'OPEN' 1* 1* 0.5 /
/
WCONPROD
 'PROD' 'OPEN' 'ORAT' 20000 4* 1000 /
/
SOURCE
{source}
/
TSTEP
{TSTEP} /
END
'''


def deck(base, r, perturb=False):
    s = sections(base)
    runspec = s['RUNSPEC'].replace('10 10 3', f'{NX * r} {NY * r} {NZ * r}')
    edit = 'EDIT\n' + PERTURBATION if perturb else ''
    return '\n'.join(['RUNSPEC', runspec, grid(r), edit,
                      'PROPS', s['PROPS'], 'SOLUTION', s['SOLUTION'],
                      'SUMMARY', schedule(r)])


def main(base_file, out):
    with open(base_file) as f:
        base = f.read()
    os.makedirs(out, exist_ok=True)
    for name, r, perturb in (('COARSE', 1, False), ('PERTURBED', 1, True),
                             ('FINE', R, False), ('FINE_PERTURBED', R, True)):
        with open(os.path.join(out, name + '.DATA'), 'w') as f:
            f.write(deck(base, r, perturb))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
