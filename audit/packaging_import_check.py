#!/usr/bin/env python3
"""Packaging import check for the ADVTIG viral workflow (supplements static_validation.py §2).

static_validation.py §2 resolves each file's imports against that file's own directory. That is
correct for entry-point scripts but not for modules inside pipeline/, which at runtime are
imported with the *entry script's* directory as sys.path[0] (e.g. `from pipeline import
dep_pre_chunker` inside pipeline/mp_demux_reads.py resolves against the repo root when
mp_metagenomic_assessment_v4.py is the entry point).

This script walks the transitive local import graph from each entry point using Python's own
importlib.util.find_spec with sys.path set as it would be at runtime ([entry dir] + stdlib only).
It executes no workflow code: find_spec only imports parent packages (pipeline/__init__.py is
empty; q_control is a namespace package). Third-party packages are reported, not resolved
(they are an environment matter, blocker B6).

Usage:  python3 audit/packaging_import_check.py [repo_root] > audit/packaging_import_check_results.txt
Exit status 0 if every local import from the canonical entry points resolves inside the repo.
"""
import ast
import importlib.util
import os
import sys
import sysconfig

sys.dont_write_bytecode = True  # find_spec imports parent packages; keep the tree free of __pycache__

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
STDLIB = set(sys.stdlib_module_names)
STDLIB_PATHS = [p for p in {sysconfig.get_paths()["stdlib"], sysconfig.get_paths()["platstdlib"],
                            os.path.join(sysconfig.get_paths()["stdlib"], "lib-dynload")} if os.path.isdir(p)]

# (entry point, canonical?)  The root map_genome_to_reads.py is a misplaced duplicate; checked for the record only.
ENTRY_POINTS = [
    ("mp_metagenomic_assessment_v4.py", True),
    ("pipeline/map_genome_to_reads.py", True),
    ("ADVTIG-untargeted/deep_cov.py", True),
    ("map_genome_to_reads.py", False),
]


def rel(p):
    return os.path.relpath(p, ROOT)


def imports_of(path):
    tree = ast.parse(open(path, encoding="utf-8").read())
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                yield a.name, None, node.lineno
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            yield node.module, [a.name for a in node.names], node.lineno


def repo_locations(top):
    hits = []
    for d, dirs, files in os.walk(ROOT):
        dirs[:] = [x for x in dirs if x not in (".git", "audit", "__pycache__")]
        if top + ".py" in files:
            hits.append(rel(os.path.join(d, top + ".py")))
        if top in dirs:
            hits.append(rel(os.path.join(d, top)) + "/")
    return sorted(hits)


def in_repo(spec):
    locs = [spec.origin] if spec.origin and spec.origin not in ("namespace", "built-in", "frozen") else []
    locs += list(spec.submodule_search_locations or [])
    return any(os.path.abspath(l).startswith(ROOT + os.sep) for l in locs if l)


def spec_file(spec):
    return spec.origin if spec.origin and os.path.isfile(spec.origin) else None


def check(entry):
    entry_path = os.path.join(ROOT, entry)
    saved_path, saved_mods = sys.path[:], set(sys.modules)
    sys.path[:] = [os.path.dirname(entry_path)] + STDLIB_PATHS
    importlib.invalidate_caches()
    ok, miss, third = [], [], set()
    seen, todo = set(), [entry_path]
    try:
        while todo:
            f = todo.pop()
            if f in seen:
                continue
            seen.add(f)
            for mod, names, line in imports_of(f):
                top = mod.split(".")[0]
                if top in STDLIB:
                    continue
                try:
                    spec = importlib.util.find_spec(mod)
                except (ModuleNotFoundError, ImportError, ValueError):
                    spec = None
                if spec is None:
                    elsewhere = repo_locations(top)
                    if elsewhere:
                        miss.append(f"{rel(f)}:{line}: import {mod} -- not on runtime sys.path; "
                                    f"exists in repo at {', '.join(elsewhere)}")
                    else:
                        third.add(top)
                    continue
                if not in_repo(spec):
                    third.add(top)
                    continue
                targets = []
                if names and spec.submodule_search_locations is not None:
                    for n in names:
                        sub = importlib.util.find_spec(f"{mod}.{n}")
                        if sub is None:
                            miss.append(f"{rel(f)}:{line}: from {mod} import {n} -- no submodule {mod}.{n}")
                        else:
                            targets.append((f"from {mod} import {n}", sub))
                else:
                    targets.append((f"import {mod}" if names is None else f"from {mod} import {', '.join(names)}", spec))
                for label, s in targets:
                    sf = spec_file(s)
                    ok.append(f"{rel(f)}:{line}: {label} -> {rel(sf) if sf else s.origin}")
                    if sf and sf.endswith(".py"):
                        todo.append(sf)
    finally:
        sys.path[:] = saved_path
        for m in set(sys.modules) - saved_mods:
            del sys.modules[m]
    return ok, miss, sorted(third), sorted(rel(s) for s in seen)


failed = False
for entry, canonical in ENTRY_POINTS:
    print(f"\n## {entry}  (sys.path[0] = {rel(os.path.dirname(os.path.join(ROOT, entry))) or '.'})"
          f"{'' if canonical else '  [NON-CANONICAL: misplaced duplicate, informational only]'}\n")
    if not os.path.isfile(os.path.join(ROOT, entry)):
        print("MISS  entry point not present")
        failed |= canonical
        continue
    ok, miss, third, files = check(entry)
    for l in ok:
        print(f"OK    {l}")
    for l in miss:
        print(f"MISS  {l}")
    print(f"\nLocal files in closure ({len(files)}): {', '.join(files)}")
    print(f"Third-party top-level packages (environment, not packaging): {', '.join(third) or 'none'}")
    print(f"RESULT {'PASS' if not miss else 'FAIL'}: {len(ok)} local imports resolved, {len(miss)} unresolved")
    failed |= canonical and bool(miss)
print(f"\nOVERALL (canonical entry points): {'FAIL' if failed else 'PASS'}")
sys.exit(1 if failed else 0)
