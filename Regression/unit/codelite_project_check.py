#!/usr/bin/env python3
# CodeLite projects against the sources on disk, no build needed.
#   The CodeLite build compiles only the files listed in src/REEF3D.project (DIVEMesh: src/DIVEMesh.project);
#   Makefile and CMakeLists.txt glob src/*.cpp. A new .cpp missing from the project links with undefined symbols
#   in CodeLite only (e.g. nhflow_momentum_SI_elm.cpp after SI P13). This check finds
#     * .cpp and .h files under src/ that the project does not list (src/coral/** included, coral/app and
#       reef3d_pch.h excluded: the CORAL app is its own executable, the PCH is set in the build settings)
#     * project entries whose file is not on disk (deleted or renamed files)
#   python3 codelite_project_check.py                      REEF3D (../../src) and DIVEMesh (../../../DIVEMesh/src, if there)
#   python3 codelite_project_check.py --reef3d R --divemesh D   other trees (e.g. a scratch clone before git format-patch)
#   python3 codelite_project_check.py --fix                adds the missing files next to the entry with the closest name
#                                                          (same folder and extension), removes the entries of missing files
#   python3 codelite_project_check.py --git                only the files git tracks (ignores untracked work in progress)
#   exit 1 with the differences (after --fix: 0)
# Architect: Hans Bihs
import argparse, os, re, subprocess, sys

here = os.path.dirname(os.path.abspath(__file__))

SKIP_DIRS = {'obj', 'Release', 'Debug', 'All', 'build', '.git', '__pycache__'}
PROJECTS = {
    # name: (project file, files excluded, folders excluded (relative to src))
    'REEF3D':   ('REEF3D.project',   {'reef3d_pch.h'}, {'coral/app'}),
    'DIVEMesh': ('DIVEMesh.project', set(),            set()),
}
ENTRY = re.compile(r'^(\s*)<File Name="([^"]+)"/>\s*$')


def sources(src, excl_files, excl_dirs, use_git):
    out = set()
    if use_git:
        r = subprocess.run(['git', '--no-optional-locks', '-C', src, 'ls-files', '--', '.'],
                           capture_output=True, text=True, check=True)
        names = r.stdout.split()
    else:
        names = []
        for d, dirs, files in os.walk(src):
            rel = os.path.relpath(d, src)
            rel = '' if rel == '.' else rel.replace(os.sep, '/')
            dirs[:] = [x for x in dirs if x not in SKIP_DIRS and (rel + '/' + x).lstrip('/') not in excl_dirs]
            names += [(rel + '/' + f).lstrip('/') for f in files]
    for n in names:
        if not n.endswith(('.cpp', '.h')) or n in excl_files:
            continue
        if any(n == e or n.startswith(e + '/') for e in excl_dirs):
            continue
        if any(p in SKIP_DIRS for p in n.split('/')[:-1]):
            continue
        out.add(n)
    return out


def entries(lines):
    ent = {}
    for i, l in enumerate(lines):
        m = ENTRY.match(l)
        if m:
            ent.setdefault(m.group(2), i)
    return ent


def sibling(name, ent):
    # the listed file of the same folder and extension with the longest common prefix, then the nearest name
    folder, ext = os.path.dirname(name), os.path.splitext(name)[1]
    best = None
    for e in ent:
        if os.path.splitext(e)[1] != ext:
            continue
        same = os.path.dirname(e) == folder
        k = 0
        while k < min(len(e), len(name)) and e[k] == name[k]:
            k += 1
        key = (same, k, e < name, e if e < name else '')
        if best is None or key > best[0]:
            best = (key, e)
    return best[1] if best else None


def check(label, src, use_git, fix):
    pfile, excl_files, excl_dirs = PROJECTS[label]
    path = os.path.join(src, pfile)
    if not os.path.isfile(path):
        print(f'{label}: no {path}')
        return False
    with open(path, encoding='utf-8') as f:
        lines = f.read().split('\n')
    ent = entries(lines)
    disk = sources(src, excl_files, excl_dirs, use_git)
    listed = {e for e in ent if e.endswith(('.cpp', '.h'))}
    missing = sorted(disk - listed)
    stale = sorted(e for e in listed if not os.path.isfile(os.path.join(src, e)))
    for n in missing:
        print(f'{label}: not in {pfile}: {n}')
    for n in stale:
        print(f'{label}: in {pfile}, not on disk: {n}')
    if not (missing or stale):
        print(f'{label}: {pfile} lists all {len(disk)} .cpp/.h files')
        return True
    if not fix:
        return False
    drop = {ent[n] for n in stale}
    after = {}
    for n in missing:
        s = sibling(n, ent)
        if s is None:
            print(f'{label}: no place for {n} in {pfile}, add it in CodeLite')
            continue
        after.setdefault(ent[s], []).append(n)
    new = []
    for i, l in enumerate(lines):
        if i in drop:
            continue
        new.append(l)
        if i in after:
            indent = ENTRY.match(l).group(1)
            new += [f'{indent}<File Name="{n}"/>' for n in sorted(after[i])]
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(new))
    print(f'{label}: {pfile} fixed: {sum(len(v) for v in after.values())} added, {len(drop)} removed')
    return all(n in sum(after.values(), []) for n in missing)


def main():
    a = argparse.ArgumentParser(description='CodeLite project files against the sources on disk')
    a.add_argument('--reef3d', default=os.path.join(here, '..', '..'), help='REEF3D repository root')
    a.add_argument('--divemesh', default=os.path.join(here, '..', '..', '..', 'DIVEMesh'),
                   help='DIVEMesh repository root (skipped if not there)')
    a.add_argument('--git', action='store_true', help='only git-tracked files')
    a.add_argument('--fix', action='store_true', help='add missing files, remove entries of missing files')
    o = a.parse_args()
    ok = check('REEF3D', os.path.normpath(os.path.join(o.reef3d, 'src')), o.git, o.fix)
    dm = os.path.normpath(os.path.join(o.divemesh, 'src'))
    if os.path.isdir(dm):
        ok = check('DIVEMesh', dm, o.git, o.fix) and ok
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
