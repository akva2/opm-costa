"""Refine a corner-point deck vertically, splitting every layer into R layers.

The refined grid is nested in the original one: lateral geometry, faults
and pillars are unchanged, and each cell is split into R cells of equal
thickness. Written for the Norne deck, it handles the keywords used there:

 * SPECGRID, DIMENS: nz scaled by R
 * ZCORN: corners interpolated linearly within each original layer
 * per-cell arrays (ACTNUM, PORO, PERMX, ...): repeated for each sublayer
 * MULTZ arrays and boxes: the value is kept on the bottom sublayer only,
   since it multiplies the transmissibility to the layer below. Faces
   inside an original layer get 1.
 * EQUALS, MULTIPLY, ADD, COPY, BOX: K ranges remapped
 * FAULTS: K ranges remapped, Z faces only between original layers
 * COMPDAT: one connection per sublayer, connection factor and Kh split
 * PINCH, MINPV: thresholds divided by R, so the same cells are pinched
   out or removed
 * WELLDIMS: maximum number of connections per well scaled by R
 * SUMMARY: removed, as block quantities would need remapping

Anything else is copied unchanged. Keywords addressing cells which are not
handled raise an error.

Usage: python refine_vertical.py <deck> <output directory> [R]
"""

import os
import re
import sys

import numpy as np

KEYWORD = re.compile(r'^([A-Z][A-Z0-9_]{0,7})\s*(--.*)?$')
SECTIONS = ('RUNSPEC', 'GRID', 'EDIT', 'PROPS', 'REGIONS', 'SOLUTION',
            'SUMMARY', 'SCHEDULE')
ARRAYS = ('ACTNUM', 'PORO', 'NTG', 'PERMX', 'PERMY', 'PERMZ', 'MULTX', 'MULTY',
          'MULTZ', 'FLUXNUM', 'FIPNUM', 'EQLNUM', 'SATNUM', 'IMBNUM', 'PVTNUM',
          'MULTNUM', 'SWATINIT', 'ZCORN', 'COORD')
BOX_RECORDS = {'EQUALS': 2, 'MULTIPLY': 2, 'ADD': 2, 'COPY': 2}  # first box item
# Keywords for transmissibility multipliers on the face to the layer below
BOTTOM_FACE = ('MULTZ',)


def tokens(text):
    """Split deck data into tokens, without comments. Text after a / ending
    a record is a comment too."""
    result = []
    for line in text.split('\n'):
        line = line.split('--')[0]
        for tok in re.findall(r"'[^']*'|/|[^\s/]+", line):
            result.append(tok)
            if tok == '/':
                break
    return result


def expand(record):
    """Expand repeat counts n*value and defaults n* into single items."""
    result = []
    for tok in record:
        m = re.fullmatch(r'(\d+)\*(.*)', tok)
        if m and not tok.startswith("'"):
            result += [m.group(2) if m.group(2) else '1*'] * int(m.group(1))
        else:
            result.append(tok)
    return result


def records(toks):
    """Split tokens of a record keyword into records, without the empty
    record ending the keyword."""
    result, current = [], []
    for tok in toks:
        if tok == '/':
            if current:
                result.append(current)
            current = []
        else:
            current.append(tok)
    return result


def compress(values):
    """Format an array with repeat counts."""
    out, i = [], 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[j + 1] == values[i]:
            j += 1
        n = j - i + 1
        out.append(f'{n}*{values[i]}' if n > 1 else values[i])
        i = j + 1
    lines = [' '.join(out[i:i + 8]) for i in range(0, len(out), 8)]
    return '\n'.join(lines) + '\n/\n'


def number(value):
    return f'{value:.10g}'


class Refiner:
    def __init__(self, deck, out, r):
        self.root = os.path.dirname(os.path.abspath(deck))
        self.out, self.r = out, r
        self.dims = None
        self.box = None     # current BOX in original indices, 1-based inclusive
        self.section = None

    def fine_k(self, k1, k2):
        """Fine layer range for original layers k1..k2."""
        return (k1 - 1) * self.r + 1, k2 * self.r

    def run(self, deck):
        self.process(os.path.abspath(deck))

    def process(self, path):
        rel = os.path.relpath(path, self.root)
        target = os.path.join(self.out, rel)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        with open(path, errors='replace') as f:
            lines = f.read().split('\n')
        with open(target, 'w') as f:
            f.write(self.transform(lines, os.path.dirname(path)))

    def transform(self, lines, directory):
        out, i = [], 0
        while i < len(lines):
            m = KEYWORD.match(lines[i].rstrip())
            if not m:
                if self.section != 'SUMMARY':
                    out.append(lines[i])
                i += 1
                continue
            kw = m.group(1)
            # data lines up to the next keyword
            j = i + 1
            while j < len(lines) and not KEYWORD.match(lines[j].rstrip()):
                j += 1
            data = '\n'.join(lines[i + 1:j])
            if kw in SECTIONS:
                self.section = kw
                out.append(lines[i])
                if kw != 'SUMMARY':
                    out.append(data)
            elif self.section == 'SUMMARY':
                pass
            else:
                out.append(self.keyword(kw, data, directory, lines[i]))
            i = j
        return '\n'.join(out)

    def keyword(self, kw, data, directory, line):
        toks = tokens(data)
        if kw == 'INCLUDE':
            name = toks[0].strip("'")
            self.process(os.path.normpath(os.path.join(self.root, name)))
            return f"INCLUDE\n '{name}' /\n"
        if kw in ('SPECGRID', 'DIMENS'):
            rec = expand(records(toks)[0])
            self.dims = [int(x) for x in rec[:3]]
            rec[2] = str(self.dims[2] * self.r)
            return f"{kw}\n {' '.join(rec)} /\n"
        if kw == 'BOX':
            rec = [int(x) for x in expand(records(toks)[0])]
            self.box = rec
            k1, k2 = self.fine_k(rec[4], rec[5])
            return f'BOX\n {rec[0]} {rec[1]} {rec[2]} {rec[3]} {k1} {k2} /\n'
        if kw == 'ENDBOX':
            self.box = None
            return line + '\n' + data
        if kw == 'WELLDIMS':
            # The maximum number of connections per well grows with R
            rec = expand(records(toks)[0])
            rec[1] = str(int(rec[1]) * self.r)
            return f"WELLDIMS\n {' '.join(rec)} /\n"
        if kw in ('PINCH', 'MINPV'):
            rec = expand(records(toks)[0])
            rec[0] = number(float(rec[0]) / self.r)
            return f"{kw}\n {' '.join(rec)} /\n"
        if kw in ARRAYS:
            return f'{kw}\n' + self.array(kw, expand(records(toks)[0]))
        if kw in BOX_RECORDS:
            return f'{kw}\n' + self.box_records(kw, records(toks)) + '/\n'
        if kw == 'FAULTS':
            return 'FAULTS\n' + self.faults(records(toks)) + '/\n'
        if kw == 'COMPDAT':
            return 'COMPDAT\n' + self.compdat(records(toks)) + '/\n'
        if kw == 'WELOPEN':
            for rec in records(toks):
                rec = expand(rec)
                if any(x not in ('0', '1*') for x in rec[2:5]):
                    raise ValueError(f'WELOPEN for connections not handled: {rec}')
        if kw in ('COMPSEGS', 'WELSEGS', 'COMPLUMP', 'WPIMULT', 'CARFIN', 'LGR'):
            raise ValueError(f'{kw} not handled')
        return line + '\n' + data

    def array(self, kw, values):
        nx, ny, nz = self.dims
        r = self.r
        if kw == 'COORD':
            return compress(values)
        if kw == 'ZCORN':
            z = np.array(values, dtype=float).reshape(nz, 2, 2 * ny, 2 * nx)
            top, bottom = z[:, 0], z[:, 1]
            fine = np.empty((nz * r, 2, 2 * ny, 2 * nx))
            for s in range(r):
                fine[s::r, 0] = top + (bottom - top) * s / r
                fine[s::r, 1] = top + (bottom - top) * (s + 1) / r
            return compress([number(v) for v in fine.ravel()])
        if self.box:
            i1, i2, j1, j2, k1, k2 = self.box
            shape = (k2 - k1 + 1, j2 - j1 + 1, i2 - i1 + 1)
        else:
            shape = (nz, ny, nx)
        a = np.array(values, dtype=object).reshape(shape)
        if kw in BOTTOM_FACE:
            fine = np.full((shape[0] * r,) + shape[1:], '1', dtype=object)
            fine[r - 1::r] = a
        else:
            fine = np.repeat(a, r, axis=0)
        return compress(list(fine.ravel()))

    def box_records(self, kw, recs):
        first = BOX_RECORDS[kw]
        out = []
        for rec in recs:
            rec = expand(rec)
            target = rec[1] if kw == 'COPY' else rec[0]
            if len(rec) < first + 6 or '1*' in rec[first + 4:first + 6]:
                # whole grid or current box, which is remapped already
                out.append(' '.join(rec) + ' /')
                continue
            k1, k2 = int(rec[first + 4]), int(rec[first + 5])
            if target.strip("'") in BOTTOM_FACE:
                for k in range(k1, k2 + 1):
                    rec[first + 4] = rec[first + 5] = str(k * self.r)
                    out.append(' '.join(rec) + ' /')
            else:
                rec[first + 4], rec[first + 5] = map(str, self.fine_k(k1, k2))
                out.append(' '.join(rec) + ' /')
        return '\n'.join(out) + '\n'

    def faults(self, recs):
        out = []
        for rec in recs:
            rec = expand(rec)
            k1, k2 = int(rec[5]), int(rec[6])
            face = rec[7].strip("'").upper()
            if face in ('Z', 'Z+'):
                ranges = [(k * self.r, k * self.r) for k in range(k1, k2 + 1)]
            elif face == 'Z-':
                ranges = [((k - 1) * self.r + 1,) * 2 for k in range(k1, k2 + 1)]
            else:
                ranges = [self.fine_k(k1, k2)]
            for a, b in ranges:
                rec[5], rec[6] = str(a), str(b)
                out.append(' '.join(rec) + ' /')
        return '\n'.join(out) + '\n'

    def compdat(self, recs):
        out = []
        for rec in recs:
            rec = expand(rec)
            rec += ['1*'] * (14 - len(rec))
            if '1*' in rec[3:5] or rec[3] == '0':
                raise ValueError(f'COMPDAT without layers not handled: {rec}')
            k1, k2 = int(rec[3]), int(rec[4])
            for idx in (7, 9):   # connection factor and Kh
                if rec[idx] != '1*':
                    rec[idx] = number(float(rec[idx]) / self.r)
            for k in range(k1, k2 + 1):
                a, b = self.fine_k(k, k)
                for kf in range(a, b + 1):
                    rec[3] = rec[4] = str(kf)
                    out.append(' '.join(rec) + ' /')
        return '\n'.join(out) + '\n'


def main(deck, out, r=2):
    Refiner(deck, out, r).run(deck)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 2)
