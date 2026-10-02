#!/usr/bin/env python3
"""Static usability checks for the ADVTIG viral workflow repository.

Does NOT run any biological analysis. It only:
  * byte-compiles every .py file (syntax check, into a temp dir)
  * resolves imports statically (AST) as Python would when the script is run directly
  * rebuilds each script's argparse parser in isolation (argparse only, no workflow imports)
  * parses every documented command line from the run-books against those parsers
  * bash -n syntax-checks the run-book code fences
  * greps for cross-stage file-name dependencies

Usage:  python3 audit/static_validation.py [repo_root]  > audit/static_validation_results.txt
"""
import ast
import argparse
import os
import py_compile
import re
import shlex
import subprocess
import sys
import tempfile

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), ".."))
STDLIB = set(sys.stdlib_module_names)


def rel(p):
    return os.path.relpath(p, ROOT)


def py_files():
    out = []
    for d, dirs, files in os.walk(ROOT):
        dirs[:] = [x for x in dirs if x not in (".git", "audit", "__pycache__")]
        out += [os.path.join(d, f) for f in files if f.endswith(".py")]
    return sorted(out)


def section(t):
    print(f"\n## {t}\n")


# 1. syntax ---------------------------------------------------------------
section("1. Syntax (py_compile)")
with tempfile.TemporaryDirectory() as tmp:
    for f in py_files():
        try:
            py_compile.compile(f, cfile=os.path.join(tmp, "x.pyc"), doraise=True)
            print(f"OK    {rel(f)}")
        except py_compile.PyCompileError as e:
            print(f"FAIL  {rel(f)}: {e.msg}")

# 2. imports --------------------------------------------------------------
section("2. Import resolution (script directory on sys.path, as when run directly)")
third_party = {}
for f in py_files():
    tree = ast.parse(open(f, encoding="utf-8").read())
    base = os.path.dirname(f)
    for node in ast.walk(tree):
        mods = []
        if isinstance(node, ast.Import):
            mods = [(a.name, None) for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            mods = [(node.module, [a.name for a in node.names])]
        for mod, names in mods:
            top = mod.split(".")[0]
            if top in STDLIB:
                continue
            local_mod = os.path.join(base, *mod.split(".")) + ".py"
            local_pkg = os.path.join(base, *mod.split("."))
            if os.path.isfile(local_mod):
                print(f"OK    {rel(f)}: import {mod} -> {rel(local_mod)}")
            elif os.path.isdir(local_pkg):
                for n in names or []:
                    sub = os.path.join(local_pkg, n + ".py")
                    print(f"{'OK  ' if os.path.isfile(sub) else 'MISS'}  {rel(f)}: from {mod} import {n} -> {rel(sub)}")
            elif top in {"pandas", "numpy", "pysam", "Bio", "matplotlib", "seaborn", "scipy", "sklearn", "openpyxl"}:
                third_party.setdefault(top, set()).add(rel(f))
            else:
                print(f"MISS  {rel(f)}: import {mod}{' (' + ', '.join(names) + ')' if names else ''} "
                      f"-- no {rel(local_mod)} or {rel(local_pkg)}/ in repo")
print("\nThird-party packages required:")
for k, v in sorted(third_party.items()):
    print(f"  {k}: {', '.join(sorted(v))}")
env_files = [p for p in ("environment.yml", "environment.yaml", "requirements.txt", "Pipfile", "Pipfile.lock",
                         "pyproject.toml", "setup.py", "conda-lock.yml", "Dockerfile", "Singularity", "Apptainer")
             if os.path.exists(os.path.join(ROOT, p))]
print(f"\nEnvironment definition files present in repo: {env_files or 'NONE'}")


# 3. parsers --------------------------------------------------------------
section("3. argparse construction (rebuilt in isolation from the script source)")


class _Stub(dict):
    """Stand-in for database_config.blastN_databases() etc. when building mp_metagenomic's parser."""


def build_parser(path):
    src = open(path, encoding="utf-8").read()
    tree = ast.parse(src)
    fn = next((n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "parse_arguments"), None)
    ns = {"argparse": argparse}
    if fn is not None:
        body = [s for s in fn.body if not (isinstance(s, ast.Return))]
        body.append(ast.Return(value=ast.Name(id="parser", ctx=ast.Load())))
        fn.body = body
        mod = ast.Module(body=[fn], type_ignores=[])
        exec(compile(ast.fix_missing_locations(mod), path, "exec"), ns)
        return ns["parse_arguments"]()
    # module-level or __main__-level parser: collect statements that touch `parser` or the stub tables
    stmts = []

    def collect(body):
        for s in body:
            if isinstance(s, ast.If):
                collect(s.body)
                continue
            seg = ast.get_source_segment(src, s) or ""
            if re.match(r"\s*(parser\s*=|parser\.add_argument|\w+_\s*=\s*\\?\s*$|random_spp_\s*=)", seg) or \
               re.match(r"random_spp_\s*=", seg) or "_libraries = database_config" in seg:
                stmts.append(seg)
    collect(tree.body)
    code = "\n".join(stmts)
    code = re.sub(r"database_config\.\w+\(\)", "_Stub({'uviral25-2':1,'uviral25':1,'c-viral29':1,'16S_16S':1})", code)
    ns["_Stub"] = _Stub
    exec(code, ns)
    return ns["parser"]


SCRIPTS = {
    "run_coverage.py": "ADVTIG-untargeted/run_coverage.py",
    "deep_cov.py": "ADVTIG-untargeted/deep_cov.py",
    "filt_low_complexity.py": "ADVTIG-untargeted/filt_low_complexity.py",
    "counter_screen.py": "ADVTIG-untargeted/counter_screen.py",
    "count_the_screen.py": "ADVTIG-untargeted/count_the_screen.py",
    "mp_metagenomic_assessment_v4.py": "mp_metagenomic_assessment_v4.py",
    "map_genome_to_reads.py": "map_genome_to_reads.py",
}
parsers = {}
for name, p in SCRIPTS.items():
    full = os.path.join(ROOT, p)
    try:
        parsers[name] = build_parser(full)
        print(f"OK    {p}: parser constructed")
    except argparse.ArgumentError as e:
        print(f"FAIL  {p}: argparse.ArgumentError: {e}")
    except Exception as e:  # noqa
        print(f"FAIL  {p}: {type(e).__name__}: {e}")
print("NOTE  collect_taxid_primary_counts.py, split_fa.py, summarise_urvdb.py have no CLI (hard-coded module-level code).")

# 4. documented commands --------------------------------------------------
section("4. Documented command lines parsed against the parsers (no execution)")


class _NoExit(argparse.ArgumentParser):
    pass


def try_parse(parser, argv):
    import io, contextlib
    err = io.StringIO()
    try:
        with contextlib.redirect_stderr(err):
            parser.parse_args(argv)
        return "OK", ""
    except SystemExit:
        return "FAIL", err.getvalue().strip().splitlines()[-1] if err.getvalue() else "SystemExit"


# Classification (added 2026-10-02 with the B3 fix). CURRENT = uncommented command lines in ADVTIG-protocol.md.
# HISTORICAL = commented-out lines in ADVTIG-protocol.md and every line of ADVTIG-untargeted/URVDB.md (a dated
# build/run log). Historical lines are still parsed against the CURRENT parsers and reported; a failure is printed
# as HFAIL (expected: they use the CLI of their date) and counted separately. Current failures stay FAIL.
counts = {("CURRENT", "OK"): 0, ("CURRENT", "FAIL"): 0, ("HISTORICAL", "OK"): 0, ("HISTORICAL", "FAIL"): 0}
cmd_re = re.compile(r"^\s*(#\s*)?(time\s+)?(\S*?)(run_coverage|deep_cov|filt_low_complexity|counter_screen|count_the_screen|"
                    r"mp_metagenomic_assessment_v4|map_genome_to_reads|collect_taxid_primary_counts|split_fa|summarise_urvdb)\.py(.*)$")
for doc in ("ADVTIG-protocol.md", "ADVTIG-untargeted/URVDB.md"):
    print(f"### {doc}  (shell variables substituted from the run-book's own preceding assignments)")
    shell_vars = {}
    for i, line in enumerate(open(os.path.join(ROOT, doc), encoding="utf-8"), 1):
        a = re.match(r'^\s*([A-Z_][A-Z0-9_]*)="([^"]*)"\s*(#.*)?$', line)
        if a:
            shell_vars[a.group(1)] = a.group(2)
            continue
        m = cmd_re.match(line)
        if not m:
            continue
        commented, prefix, script, rest = bool(m.group(1)), m.group(3), m.group(4) + ".py", m.group(5)
        rest = rest.split(" #")[0]
        rest = re.sub(r"\$\{(\w+)\}", lambda mm: shell_vars.get(mm.group(1), mm.group(0)), rest)
        try:
            argv = shlex.split(rest)
        except ValueError:
            argv = rest.split()
        expected = SCRIPTS.get(script)
        loc = (prefix + script)
        if script in ("collect_taxid_primary_counts.py", "split_fa.py", "summarise_urvdb.py"):
            status, msg = ("OK", "no CLI") if not argv else ("WARN", "args given to script without CLI")
        elif script not in parsers:
            status, msg = "FAIL", "parser cannot be constructed"
        else:
            status, msg = try_parse(parsers[script], argv)
        path_ok = ""
        if loc.startswith("./"):
            path_ok = "path OK" if os.path.isfile(os.path.join(ROOT, loc[2:])) else f"PATH MISSING in repo ({loc})"
        else:
            path_ok = f"external path ({loc})"
        cls = "HISTORICAL" if (commented or doc.endswith("URVDB.md")) else "CURRENT"
        counts[(cls, "FAIL" if status == "FAIL" else "OK")] += 1
        shown = "HFAIL" if (status == "FAIL" and cls == "HISTORICAL") else status
        print(f"{shown:5} L{i:<4} [{cls}]{' [commented]' if commented else ''} {script} {' '.join(argv)[:90]} | {path_ok}"
              + (f" | {msg}" if msg else ""))
print(f"\nSUMMARY §4: CURRENT commands: {counts[('CURRENT', 'OK')]} OK, {counts[('CURRENT', 'FAIL')]} FAIL; "
      f"HISTORICAL commands: {counts[('HISTORICAL', 'OK')]} parse under current CLI, "
      f"{counts[('HISTORICAL', 'FAIL')]} HFAIL (intentionally not valid under the current CLI)")

# 5. bash -n --------------------------------------------------------------
section("5. bash -n on run-book code fences")
for doc in ("ADVTIG-protocol.md", "ADVTIG-untargeted/URVDB.md"):
    text = open(os.path.join(ROOT, doc), encoding="utf-8").read()
    blocks = re.findall(r"```bash\n(.*?)```", text, flags=re.S)
    for bi, b in enumerate(blocks):
        r = subprocess.run(["bash", "-n"], input=b, text=True, capture_output=True)
        print(f"{'OK  ' if r.returncode == 0 else 'FAIL'}  {doc} block {bi}: {r.stderr.strip()[:300] or 'syntax OK'}")

# 6. cross-stage names ----------------------------------------------------
section("6. Cross-stage file-name dependencies")
writes = subprocess.run(["grep", "-rn", "--exclude-dir=audit", "MetaFilt-Fourth-Pass\\|MetaFilt-Fifth-Pass", "--include=*.py", ROOT],
                        text=True, capture_output=True).stdout
for l in writes.splitlines():
    print("  " + l.replace(ROOT + "/", ""))
# Stage 6 (count_the_screen.py, --skip) merged output vs stage 7 (collect_taxid_primary_counts.py) FILTER_CSV input.
_s6 = open(os.path.join(ROOT, "ADVTIG-untargeted/count_the_screen.py"), encoding="utf-8").read()
_s7 = open(os.path.join(ROOT, "ADVTIG-untargeted/collect_taxid_primary_counts.py"), encoding="utf-8").read()
s6_out = set(re.findall(r'out_save_merge\s*=\s*f?"[^"]*/(MetaFilt-[^"/]+\.csv)"', _s6))
s7_in = set(re.findall(r'FILTER_CSV\s*=\s*f?"[^"]*/(MetaFilt-[^"/]+\.csv)"', _s7))
print(f"\nstage 6 writes: {sorted(s6_out)}; stage 7 reads: {sorted(s7_in)}")
print(f"{'OK  ' if s6_out and s7_in and s7_in <= s6_out else 'FAIL'}  stage 6 output filename "
      f"{'matches' if s6_out and s7_in and s7_in <= s6_out else 'does NOT match'} stage 7 input filename")

# 7. inventory of supporting material ------------------------------------
section("7. Supporting material present in repo")
for label, pats in {
    "sample sheets / configs": (".csv", ".txt", ".tsv", ".yml", ".yaml", ".json", ".ini", ".cfg"),
}.items():
    found = [rel(os.path.join(d, f)) for d, ds, fs in os.walk(ROOT) if ".git" not in d and "audit" not in d
             for f in fs if f.endswith(pats) and not f.startswith("WORKFLOW_")]
    print(f"{label}: {found or 'NONE'}")
for p in ("pipeline", "q_control", "configs", "docs", "tests", "envs", "references", "Makefile", "Snakefile"):
    print(f"{p}: {'present' if os.path.exists(os.path.join(ROOT, p)) else 'ABSENT'}")
print(f"README.md bytes: {os.path.getsize(os.path.join(ROOT, 'README.md'))}")


# 8. targeted B3 / B4 checks (added 2026-10-02) ----------------------------
section("8. Targeted checks for B3 (run_coverage.py CLI) and B4 (stage 6 -> 7 file name)")
# Current CLI decision (2026-10-02): extraction = -x/--extract-fastq; edit distance = --edit-distance FLOAT
# (default 0.15); run_coverage.py has no -e option at all. Historical -e forms are classified in §4.


def _check(label, ok, detail=""):
    print(f"{'OK  ' if ok else 'FAIL'}  {label}{' -- ' + detail if detail else ''}")


rc = parsers.get("run_coverage.py")
_check("run_coverage.py parser constructs (rebuilt from parse_arguments source)", rc is not None)
if rc is not None:
    try:
        h = rc.format_help()
        _check("run_coverage.py help text generates", "-x, --extract-fastq" in h and "--edit-distance EDIT_DISTANCE" in h
               and not re.search(r"(^|\s)-e[\s,]", h), "help lists -x/--extract-fastq and --edit-distance; no -e option")
    except Exception as e:  # noqa
        _check("run_coverage.py help text generates", False, f"{type(e).__name__}: {e}")
    base = ["-s", "s.csv", "-d", "D", "-o", "D", "--database", "uviral25-2"]

    def _ns(extra):
        import io, contextlib
        with contextlib.redirect_stderr(io.StringIO()):
            try:
                return rc.parse_args(base + extra)
            except SystemExit:
                return None
    n = _ns([])
    _check("default edit-distance threshold is 0.15 (omitting --edit-distance == --edit-distance 0.15)",
           n is not None and n.edit_distance_threshold == 0.15 and vars(n) == vars(_ns(["--edit-distance", "0.15"])))
    n = _ns(["--edit-distance", "0.15"])
    _check("'--edit-distance 0.15' parses", n is not None and n.edit_distance_threshold == 0.15 and not n.extract_fastq)
    n = _ns(["--edit-distance", "0.10"])
    _check("explicit other value '--edit-distance 0.10' still accepted", n is not None and n.edit_distance_threshold == 0.10)
    n = _ns(["-x"])
    _check("'-x' accepted as extract-FASTQ mode", n is not None and n.extract_fastq and n.edit_distance_threshold == 0.15)
    n = _ns(["--extract-fastq"])
    _check("long form '--extract-fastq' accepted", n is not None and n.extract_fastq)
    for combo in (["--edit-distance", "0.15", "-x"], ["-x", "--edit-distance", "0.15"]):
        n = _ns(combo)
        _check(f"'{' '.join(combo)}' coexist", n is not None and n.extract_fastq and n.edit_distance_threshold == 0.15)
    for combo, attr in ((["--edit-distance", "0.15", "-r"], "run_alignment"), (["--edit-distance", "0.15", "-l"], "label_species"),
                        (["--edit-distance", "0.15", "-a", "-x"], "asm5")):
        n = _ns(combo)
        _check(f"'{' '.join(combo)}' parses", n is not None and getattr(n, attr) and n.edit_distance_threshold == 0.15)
    _check("'-e 0.15' rejected (no -e option in current CLI)", _ns(["-e", "0.15"]) is None)
    _check("historical bare '-e' (old extract form) rejected, as expected under the current CLI", _ns(["-e"]) is None)

    # The real file, end to end: `run_coverage.py --help` through its own __main__ block. pandas/pysam/numpy/Bio are
    # not installed here (blocker B6), so inert stand-in modules are injected for those imports only; --help exits
    # inside parse_arguments() before any analysis code runs.
    stub = r"""
import sys, types, runpy
class _M(types.ModuleType):
    def __getattr__(self, k):
        return _M(k)
    def __call__(self, *a, **k):
        return _M('x')
for m in ['pysam', 'pandas', 'numpy', 'Bio', 'Bio.Align', 'Bio.pairwise2', 'Bio.SeqIO']:
    sys.modules.setdefault(m, _M(m))
sys.argv = [sys.argv[1], '--help']
runpy.run_path(sys.argv[0], run_name='__main__')
"""
    rp = os.path.join(ROOT, "ADVTIG-untargeted/run_coverage.py")
    r = subprocess.run([sys.executable, "-B", "-c", stub, rp], text=True, capture_output=True)
    _check("real `run_coverage.py --help` exits 0 without argparse exception (third-party imports stubbed)",
           r.returncode == 0 and "--extract-fastq" in r.stdout and "--edit-distance" in r.stdout and "argparse" not in r.stderr,
           f"rc={r.returncode}{'; ' + r.stderr.strip().splitlines()[-1] if r.stderr.strip() else ''}")

_cur_rc = []
for _i, _line in enumerate(open(os.path.join(ROOT, "ADVTIG-protocol.md"), encoding="utf-8"), 1):
    _m = cmd_re.match(_line)
    if _m and not _m.group(1) and _m.group(4) == "run_coverage":
        _cur_rc.append((_i, shlex.split(_m.group(5).split(" #")[0])))
_bad = [i for i, a in _cur_rc if "-e" in a]
_check(f"no CURRENT run_coverage.py command uses -e ({len(_cur_rc)} current commands scanned)", bool(_cur_rc) and not _bad,
       f"offending lines: {_bad}" if _bad else "")
_check(f"every CURRENT run_coverage.py command passes --edit-distance 0.15",
       bool(_cur_rc) and all("--edit-distance" in a and a[a.index("--edit-distance") + 1] == "0.15" for _, a in _cur_rc))
_check(f"CURRENT extraction commands use '--edit-distance 0.15 -x'",
       any("-x" in a for _, a in _cur_rc) and all(" ".join(a).endswith("--edit-distance 0.15 -x") for _, a in _cur_rc if "-x" in a))
_check("B4: stage 6 output filename equals stage 7 input filename", bool(s6_out) and bool(s7_in) and s7_in <= s6_out,
       f"stage 6 {sorted(s6_out)} / stage 7 {sorted(s7_in)}")

