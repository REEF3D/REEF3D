#!/usr/bin/env python3
# Layering of REEF3D::MAST (src/mast/): checks the include rules of src/mast/README.md, no build needed.
#   * a folder includes only from the folders it may use (table ALLOWED below = README "Layers")
#   * quoted includes outside src/mast/reef3d/ start with "mast/" (no REEF3D headers there)
#   * <mpi.h> only in src/mast/reef3d/
#   * file names are unique across the folders
#   python3 mast_layering_check.py     exit 1 with the violations
# Architect: Hans Bihs
import os, re, sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.normpath(os.path.join(here, '..', '..', 'src', 'mast'))

PHYSICS = ['core', 'structure', 'linear', 'wind', 'aero', 'control', 'hydro', 'panel',
           'mooring', 'pto', 'soil', 'ice']
ALLOWED = {
    'core':      ['core'],
    'structure': ['core'],
    'wind':      ['core'],
    'control':   ['core'],
    'hydro':     ['core'],
    'linear':    ['core', 'structure'],
    'mooring':   ['core', 'structure'],
    'pto':       ['core', 'structure'],
    'soil':      ['core', 'structure'],
    'ice':       ['core', 'structure'],
    'aero':      ['core', 'structure', 'wind'],
    'panel':     ['core', 'hydro'],
    'io':        PHYSICS + ['io'],
    'loadcase':  PHYSICS + ['io', 'loadcase'],
    'post':      PHYSICS + ['io', 'post'],
    'app':       PHYSICS + ['io', 'loadcase', 'post', 'app'],
    'reef3d':    PHYSICS + ['io', 'loadcase', 'post', 'reef3d'],
}

inc = re.compile(r'^\s*#\s*include\s*([<"])([^>"]+)[>"]')

def strip_comments(text):
    text = re.sub(r'/\*.*?\*/', lambda m: '\n'*m.group(0).count('\n'), text, flags=re.S)
    return re.sub(r'//[^\n]*', '', text)

errors = []
names = {}
nfiles = 0
if not os.path.isdir(root):
    print('no src/mast')
    sys.exit(1)
for d in sorted(os.listdir(root)):
    path = os.path.join(root, d)
    if not os.path.isdir(path):
        continue
    if d not in ALLOWED:
        errors.append('folder %s/ is not in the layer table (add it to ALLOWED and to src/mast/README.md)' % d)
        continue
    for dirpath, _, files in os.walk(path):
        for fn in sorted(files):
            if not fn.endswith(('.h', '.hpp', '.cpp')):
                continue
            nfiles += 1
            rel = os.path.relpath(os.path.join(dirpath, fn), root)
            if fn in names:
                errors.append('%s: file name also used in %s' % (rel, names[fn]))
            names[fn] = rel
            for line in strip_comments(open(os.path.join(dirpath, fn), errors='replace').read()).split('\n'):
                m = inc.match(line)
                if not m:
                    continue
                kind, target = m.group(1), m.group(2).strip()
                if kind == '<':
                    if target.startswith('mpi') and d != 'reef3d':
                        errors.append('%s: <%s> outside reef3d/' % (rel, target))
                    continue
                if target.startswith('mast/'):
                    parts = target.split('/')
                    if len(parts) < 3 or parts[1] not in ALLOWED:
                        errors.append('%s: "%s" is not a MAST folder' % (rel, target))
                    elif parts[1] not in ALLOWED[d]:
                        errors.append('%s: %s/ may not include from %s/ ("%s")' % (rel, d, parts[1], target))
                elif d != 'reef3d':
                    errors.append('%s: "%s" (REEF3D headers only in reef3d/; MAST includes start with "mast/")' % (rel, target))

for e in errors:
    print('FAIL  ' + e)
print('%d files in %d folders, %d violations' % (nfiles, len(ALLOWED), len(errors)))
sys.exit(1 if errors else 0)
