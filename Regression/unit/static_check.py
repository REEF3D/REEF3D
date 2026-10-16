#!/usr/bin/env python3
# Static data in the sources (the contract of the mesh refinement: every grid has its own lexer,
# fdm and solver objects, so data that belongs to one grid must not be static).  Lists every
# non-constant static or thread_local variable (class members, function-local and file-scope
# statics) and every non-constant variable at namespace scope in src/*.h, src/*.cpp and the
# REEF3D::CORAL folders below src/coral/, and compares
# them with static_allowlist.txt (one per line: file | declaration | why it is safe).
#   python3 static_check.py            check: exit 1 with the new ones (and the stale entries)
#   python3 static_check.py --list     print the current list in the allowlist format
# No build needed.  A new static: make it a member of the object of the grid, or add it to the
# allowlist with the reason it holds no grid data (a process-wide setting, an output writer, ...).
# Architect: Hans Bihs
import os, re, sys

here = os.path.dirname(os.path.abspath(__file__))
src = os.path.normpath(os.path.join(here, '..', '..', 'src'))
allow = os.path.join(here, 'static_allowlist.txt')

decl = re.compile(r'^\s*(?:inline\s+)?(?:static|thread_local)\b(?:\s+(?:inline|thread_local|static))*\s+(?!const\b|constexpr\b)([^;]*);\s*$')
func = re.compile(r'\)\s*(?:const)?\s*(?:override|final|noexcept)?\s*$')

def strip_comments(text):
    text = re.sub(r'/\*.*?\*/', lambda m: '\n'*m.group(0).count('\n'), text, flags=re.S)
    return re.sub(r'//[^\n]*', '', text)

def scan():
    found = []
    files = sorted(os.listdir(src))
    # REEF3D::CORAL lives in folders below src/coral (named by their path from src, e.g. coral/core/x.cpp)
    for dirpath, _, names in sorted(os.walk(os.path.join(src, 'coral'))):
        files += sorted(os.path.relpath(os.path.join(dirpath, n), src) for n in names)
    for fn in files:
        if not fn.endswith(('.h', '.cpp')):
            continue
        text = strip_comments(open(os.path.join(src, fn), errors='replace').read())
        for line in text.split('\n'):
            m = decl.match(line)
            if not m:
                continue
            body = ' '.join(m.group(1).split())
            if 'const ' in body.split('=')[0] and not '*' in body.split('=')[0]:
                continue
            name = body.split('=')[0].strip()
            # functions (static member functions, static free functions): a parameter list
            if '(' in name and not re.search(r'\w\s*\[', name):
                continue
            found.append((fn, body))
    return found

def norm(s):
    return ' '.join(s.split())

def read_allow():
    entries = []
    if os.path.exists(allow):
        for line in open(allow):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = [x.strip() for x in line.split('|')]
            if len(parts) >= 2:
                entries.append((parts[0], norm(parts[1]), parts[2] if len(parts) > 2 else ''))
    return entries

found = scan()
if '--list' in sys.argv:
    for fn, body in found:
        print('%s | %s | ' % (fn, body))
    sys.exit(0)

allowed = {(f, d) for f, d, _ in read_allow()}
new = [(f, d) for f, d in found if (f, norm(d)) not in allowed]
present = {(f, norm(d)) for f, d in found}
stale = [(f, d) for f, d, _ in read_allow() if (f, d) not in present]
for f, d in new:
    print('NEW static: %s | %s' % (f, d))
for f, d in stale:
    print('stale allowlist entry (no longer in the sources): %s | %s' % (f, d))
print('%d statics, %d allowed, %d new, %d stale' % (len(found), len(found) - len(new), len(new), len(stale)))
sys.exit(1 if new else 0)
