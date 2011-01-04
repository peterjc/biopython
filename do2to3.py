"""Helper script for building and installing Biopython on Python 3.

Note that we can't just use distutils.command.build_py function build_py_2to3
in setup.py since (as far as I can see) that does not allow us to alter the
2to3 options. In particular, we need to turn off the long fixer for some of
our files.

This code is intended to be called from setup.py automatically under Python 3,
and is not intended for end users. The basic idea follows the approach taken
by NumPy with their setup.py file calling tools/py3tool.py to do the 2to3
conversion automatically.

This calls the lib2to3 library functions to convert the Biopython source code
from Python 2 to Python 3, tracking changes to files so that unchanged files
need not be reconverted making development much easier (i.e. if you edit one
source file, doing 'python setup.py install' will only reconvert the one file).
This is done by recording the md5 checksum of the Python 2 source, and the
converted Python 3 version. This should be clever enough to reconvert files
as needed when switching branches (common during development).

NOTE - This is intended to be run under Python 3 (not under Python 2), but
care has been taken to make it run under Python 2 enough to give a clear error
message. In particular, this meant avoiding with statements etc.
"""
import sys
if sys.version_info[0] < 3:
    sys.stderr.write("Please run this under Python 3\n")
    sys.exit(1)

import shutil
import os
import subprocess
import lib2to3
import lib2to3.main
from io import StringIO
from hashlib import md5


def run2to3(filename):
    stderr = sys.stderr
    handle = StringIO()
    try:
        #Want to capture stderr (otherwise too noisy)
        sys.stderr = handle
        #TODO - Configurable options per file?
        args = ["--nofix=long", "--no-diffs", "-n", "-w"]
        e = lib2to3.main.main("lib2to3.fixes", args + [filename])
        if e != 0:
            sys.stderr = stderr
            sys.stderr.write(handle.get_value())
            raise RuntimeError("Error %i from 2to3")
        #And again for any doctests,
        e = lib2to3.main.main("lib2to3.fixes", args + ["-d", filename])
        if e != 0:
            sys.stderr = stderr
            sys.stderr.write(handle.get_value())
            raise RuntimeError("Error %i from 2to3 (doctests)")
    finally:
        #Restore stderr
        sys.stderr = stderr

def file_md5(filename):
    handle = open(filename, "rb")
    answer = md5(handle.read()).hexdigest()
    handle.close()
    return answer

def load_md5_cache(filename):
    old = dict()
    new = dict()
    handle = open(filename)
    for line in handle:
        line = line.rstrip()
        if not line or line.startswith("#"): continue
        parts = line.rstrip().split("\t")
        if len(parts) == 2:
            old[parts[0]]=parts[1]
        elif len(parts) == 3:
            old[parts[0]]=parts[1]
            new[parts[0]]=parts[2]
    handle.close()
    return old, new
        
def save_md5_cache(filename, old, new):
    handle = open(filename, "w")
    handle.write("#This is an automatically generated file.\n")
    handle.write("#Each line is tab separated, starting with a filename,\n")
    handle.write("#followed by the original file's md5 and (for Python\n")
    handle.write("#files only) the converted file's md5 checksum.\n")
    handle.write("#These checksums are used to avoid re-converting any\n")
    handle.write("#unchanged Python files with 2to3 via setup.py\n")
    for f, m in old.items():
        if f in new:
            handle.write("%s\t%s\t%s\n" % (f,m,new[f]))
        else:
            handle.write("%s\t%s\n" % (f,m))
    handle.close()

def do_update(py2folder, py3folder, md5cache, verbose=False):
    if os.path.isfile(md5cache):
        old_md5, new_md5 = load_md5_cache(md5cache)
        #Check and remove any stale values
        for f in old_md5:
            if not os.path.isfile(os.path.join(py2folder, f))\
            or not os.path.isfile(os.path.join(py3folder, f)):
                del old_md5[f]
        for f in new_md5:
            if not os.path.isfile(os.path.join(py2folder, f))\
            or not os.path.isfile(os.path.join(py3folder, f)):
                del new_md5[f]
        save_md5_cache(md5cache, old_md5, new_md5)
    else:
        old_md5 = dict()
        new_md5 = dict()
    if not os.path.isdir(py2folder):
        raise ValueError("Python 2 folder %r does not exist" % py2folder)
    if not os.path.isdir(py3folder):
        os.mkdir(py3folder)
    for dirpath, dirnames, filenames in os.walk(py2folder):
        if verbose: print("Processing %s" % dirpath)
        relpath = os.path.relpath(dirpath, py2folder)
        for d in dirnames:
            new = os.path.join(py3folder, relpath, d)
            if not os.path.isdir(new):
                os.mkdir(new)
        for f in filenames:
            if f.startswith("."):
                #Ignore hidden files
                continue
            elif f.endswith("~") or f.endswith(".bak") or f.endswith(".swp"):
                #Ignore backup files
                continue
            elif f.endswith(".pyc") or f.endswith("$py.class"):
                #Ignore compiled python
                continue
            rel = os.path.join(relpath, f)
            if rel[:2] == "./":
                #This is just to give cleaner filenames
                rel = rel[2:]
            old = os.path.join(py2folder, rel)
            m = file_md5(old)
            new = os.path.join(py3folder, rel)
            if rel in old_md5 and m == old_md5[rel] \
            and os.path.isfile(new) and file_md5(new) == new_md5.get(rel, m):
                if verbose: print("Unchanged: %s" % new)
                continue
            #Python, C code, data files, etc - copy
            shutil.copy2(old, new)
            if f.endswith(".py"):
                #Also run 2to3 on it
                print("Converting %s" % new)
                run2to3(new)
                new_md5[rel] = file_md5(new)
            else:
                if verbose: print("Updated %s" % new)
                assert rel not in new_md5
            old_md5[rel] = m
            #Bit heavy handed to resave after each conversion/copy,
            #but useful if the task is interrupted and resumed later.
            save_md5_cache(md5cache, old_md5, new_md5)
            
def main(children=["Bio", "BioSQL", "Tests", "Scripts", "Doc"]):
    #Want to use different folders for Python 3.1, 3.2, etc
    #since the 2to3 libraries have changed so the conversion
    #may differ slightly.
    python3_source = "build/py%i.%i" % sys.version_info[:2]
    md5_cache_template = python3_source + "_%s.md5"

    if not os.path.isdir("build"):
        os.mkdir("build")
    if not os.path.isdir(python3_source):
        os.mkdir(python3_source)

    for child in children:
        print("Processing %s" % child)
        do_update(child,
                  os.path.join(python3_source, child),
                  md5_cache_template % child)
    print("Python 2to3 processing done.")
              
if __name__ == "__main__":
    main()
