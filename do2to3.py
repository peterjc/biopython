"""Helper script for building and installing Biopython on Python 3.

Note that we can't just use distutils.command.build_py function build_py_2to3
in setup.py since (as far as I can see) that does not allow us to alter the
2to3 options. In particular, we need to turn off the long fixer for some of
our files.

Furthermore, we have also decided to make some custom changes as part of the
Python 2 to Python 3 migration - adopting lower case PEP8 conformant module
names. This means we also need to rename files and folders to lower case,
and adjust the import lines. e.g. 'import Bio' becomes 'import bio'.

This code is intended to be called from setup.py automatically under Python 3,
and is not intended for end users. The basic idea follows the approach taken
by NumPy with their setup.py file calling tools/py3tool.py to do the 2to3
conversion automatically.

This calls the lib2to3 library functions to convert the Biopython source code
from Python 2 to Python 3, tracking changes to files so that unchanged files
need not be reconverted making development much easier (i.e. if you edit one
source file, doing 'python setup.py install' will only reconvert the one file).
This is done by the last modified date stamps (which will be updated by git if
you switch branches).

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
import re
import lib2to3.main
from io import StringIO

#Initial idea was to define a new fixer, but how do I invoke it?
#from lib2to3.fixes.fix_imports import FixImports
#MAPPING = {'Bio' : 'bio'}
#class LowerCaseNamespace(FixImports):
#    mapping = MAPPING

#Second idea was to monkey patch an existing import fixer.
#This works, but sadly only seems intended to handle the top
#level renames (in our case, only simple imports from Bio or
#BioSQL were changed to lower case).
#from lib2to3.fixes import fix_imports2 
#for name in ["Bio", "Bio.Alphabet", "Bio.Seq", "Bio.SeqRecord", "Bio.SeqFeature", "BioSQL"]:
#    fix_imports2.MAPPING[name] = name.lower()

#Third idea - brute force hackery, see OLD_NAMES list
OLD_NAMES = [] #Updated later!

def strip_comment(text):
    if "#" in text:
        return text[:text.find("#")].strip()
    else:
        return text.strip()

def hack_file_import_lines(f):
    """Evil hack to change our import lines to use lower case..."""
    global OLD_NAMES

    if f.endswith("/__init__.py"):
        m = f[:-12].split(os.path.sep)
    elif f.endswith(".py"):
        m = f[:-3].split(os.path.sep)
    else:
        assert False, f
    if f.endswith("/__init__.py"):
        b = m = ".".join(m[2:])
    else:
        b = ".".join(m[2:-1])
        m = ".".join(m[2:])
    #assert m.startswith("bio"), ("%r from %s" % (m, f))

    TEMPLATES = [",%s.", "=%s.", "[%s.", "(%s ", "(%s.", " %s.", " %s ", " %s)", " %s,"]

    #Top level imports:
    NAMES = list(OLD_NAMES)
    #Relative imports:
    for name in OLD_NAMES:
        if name.lower().startswith(b + "."):
            #print("Adding %s as a local import" % name)
            NAMES.append(name[len(b)+1:])
    #Now turn them into regular expressions:
    re_plain_import = re.compile("^import (%s)$" % "|".join(NAMES))
    re_from_import = re.compile("^from (%s) import (.+)$" % "|".join(NAMES))

    file_mapping = set()
    doc_mapping = set()

    h = open(f, "r")
    lines = list(h)
    h.close()

    h = open(f, "w")
    in_triple_quote = False
    for line in lines:
        assert line.count('"""') <= 2, line
        core = strip_comment(line)

        if line.count('"""') == 1:
            in_triple_quote = not in_triple_quote
            h.write(line)
            #Start/end/one line docstring => clear doc_mapping,
            doc_mapping = set()
        elif core.startswith(">>> "):
            assert in_triple_quote
            core = core[4:]
            #Will need to update rest of this docstring post import...
            if re_plain_import.match(core):
                if " as " in core:
                    #No need to change use of the " as " name later
                    i = line.find(" as ")
                    h.write(line[:i].lower() + line[i:])
                else:
                    name = core[7:].strip()
                    assert name in NAMES, \
                        "Hmm. %r from file %s, Line:\n%s" % (name, f, line)
                    doc_mapping.add(name)
                    h.write(line.lower()) #TODO - Preserve case of any comment?
            elif re_from_import.match(core):
                base = line.split("from ",1)[1].split(" import ",1)[0].strip()
                assert base in NAMES
                if " as " in core:
                    name = line.split(" import ",1)[1].split(" as ",1)[0]
                    if base + "." + name in NAMES:
                        i = line.find(" as ")
                    else:
                        i = line.find(" import ")
                    h.write(line[:i].lower() + line[i:])
                    #Don't need to change later use of the 'as' name
                else:
                    names = [x.strip() for x in line.split(" import ",1)[1].split(",")]
                    new_names = []
                    for name in names:
                        if base + "." + name in NAMES:
                            doc_mapping.add(name)
                            new_names.append(name.lower())
                        else:
                            new_names.append(name)
                    h.write(line.split(" import ",1)[0].lower() \
                                + " import " + ", ".join(new_names) + "\n")
            else:
                #Boring line; do we need to apply any import replacements?
                for name in doc_mapping:
                    for template in TEMPLATES:
                        x = template % name
                        line = line.replace(x, x.lower())
                    if line.startswith(name + "."):
                        line = name.lower() + line[len(name):]
                h.write(line)
        elif re_plain_import.match(core):
            if " as " in core:
                name = line.split("import ",1)[1].split(" as ",1)[0]
                assert name in NAMES
                i = line.find(" as ")
                h.write(line[:i].lower() + line[i:])
                #Don't need to change later use of the 'as' name
            else:
                name = core[7:].strip()
                assert name in NAMES, "Hmm. %r from file %s, Line:\n%s" % (name, f, line)
                file_mapping.add(name)
                h.write(line.lower()) #TODO - Preserve case of any comment?
        elif core.startswith("import "):
            #This could be a local import, e.g. 'import FastaIO'
            #in Bio/SeqIO/__init__.py should become 'import fastaio'
            if " as " in core:
                name = core.split("import ",1)[1].split(" as ")[0].strip()
                if (m + "." + name).lower() in [x.lower() for x in NAMES]:
                    i = line.find(" as ")
                    h.write(line[:i].lower() + line[i:])
                    #Don't need to change later use of the 'as' name
                else:
                    h.write(line)
            else:
                name = core.split("import ",1)[1].strip()
                if (m + "." + name).lower() in [x.lower() for x in NAMES]:
                    h.write(line.lower())
                    file_mapping.add(name)
                else:
                    #Nope, not one of our imports
                    h.write(line)
        elif re_from_import.match(core):
            base = line.split("from ",1)[1].split(" import ",1)[0].strip()
            assert base in NAMES
            if " as " in core:
                name = line.split(" import ",1)[1].split(" as ",1)[0]
                if base + "." + name in NAMES:
                    i = line.find(" as ")
                else:
                    i = line.find(" import ")
                h.write(line[:i].lower() + line[i:])
                #Don't need to change later use of the 'as' name
            else:
                names = [x.strip() for x in line.split(" import ",1)[1].split(",")]
                new_names = []
                for name in names:
                    #Don't collect object/function names
                    if base + "." + name in NAMES:
                        file_mapping.add(name)
                        new_names.append(name.lower())
                    else:
                        new_names.append(name)
                h.write(line.split(" import ",1)[0].lower() \
                                + " import " + ", ".join(new_names) + "\n")
        else:
            #Boring line; do we need to apply any import replacements?
            for name in file_mapping:
                for template in TEMPLATES:
                    x = template % name
                    line = line.replace(x, x.lower())
                if line.startswith(name + "."):
                    line = name.lower() + line[len(name):]
            h.write(line)
    h.close()
    #TODO - Work out why this fails (commented out to allow
    #all later files to be processed, to spot other issues):
    if in_triple_quote:
        print("Triple quote confused in %s" % f)
        #raise ValueError("Triple quote confused in %s" % f)

def run2to3(filenames):
    stderr = sys.stderr
    for filename in filenames:
        sys.stderr = stderr
        print("Converting %s" % filename)
        #First, our evil hackery of the import lines:
        hack_file_import_lines(filename)
        #TODO - Configurable options per file?
        try:
            #Want to capture stderr (otherwise too noisy)
            handle = StringIO()
            sys.stderr = handle
            args = ["--nofix=long", "--no-diffs", "-n", "-w"]
            e = lib2to3.main.main("lib2to3.fixes", args + [filename])
            if e != 0:
                sys.stderr = stderr
                sys.stderr.write(handle.getvalue())
                os.remove(filename) #Don't want a half edited file!
                raise RuntimeError("Error %i from 2to3 on %s" \
                                   % (e, filename))
            #And again for any doctests,
            e = lib2to3.main.main("lib2to3.fixes", args + ["-d", filename])
            if e != 0:
                sys.stderr = stderr
                sys.stderr.write(handle.getvalue())
                os.remove(filename) #Don't want a half edited file!
                raise RuntimeError("Error %i from 2to3 (doctests) on %s" \
                                   % (e, filename))
        except KeyboardInterrupt:
            sys.stderr = stderr
            sys.stderr.write("Interrupted during %s\n" % filename)
            os.remove(filename) #Don't want a half edited file!
            for filename in filenames:
                if os.path.isfile(filename):
                    #Don't want uncoverted files left behind:
                    os.remove(filename)
            sys.exit(1)
        finally:
            #Restore stderr
            sys.stderr = stderr

def get_module_names(py2folder):
    for top in ["Bio", "BioSQL"]:
        yield top
        for dirpath, dirnames, filenames in os.walk(os.path.join(py2folder, top)):
            for f in filenames:
                if dirpath.startswith("./"):
                    dirpath = dirpath[2:]
                parts = dirpath.split(os.path.sep)
                if f == "__init__.py":
                    pass
                elif f.endswith(".py"):
                    parts.append(f[:-3])
                else:
                    continue
                yield ".".join(parts)
    #Special cases like C modules,
    yield "Bio.PDB.mmCIF.MMCIFlex"

def do_update(py2folder, py3folder, verbose=False):
    if not os.path.isdir(py2folder):
        raise ValueError("Python 2 folder %r does not exist" % py2folder)
    if not os.path.isdir(py3folder):
        os.mkdir(py3folder)
    #First remove any files from the 3to2 conversion which no
    #longer exist in the Python 2 origin (only expected to happen
    #on a development machine).
    for dirpath, dirnames, filenames in os.walk(py3folder):
        relpath = os.path.relpath(dirpath, py3folder)
        for d in dirnames:
            #TODO - Handle case changes... tricky...
            new = os.path.join(py3folder, relpath, d)
            old = os.path.join(py2folder, relpath, d)
            if not os.path.isdir(old):
                print("Removing %s" % new)
                shutil.rmtree(new)
        for f in filenames:
            new = os.path.join(py3folder, relpath, f)
            old = os.path.join(py2folder, relpath, f)
            if not os.path.isfile(old):
                print("Removing %s" % new)
                os.remove(new)
    #Check all the Python 2 original files have been copied/converted
    #Note we need to do all the conversions *after* copying the files
    #so that 2to3 can detect local imports successfully.
    to_convert = []
    for dirpath, dirnames, filenames in os.walk(py2folder):
        if verbose: print("Processing %s" % dirpath)
        relpath = os.path.relpath(dirpath, py2folder)
        #This is just to give cleaner filenames
        if relpath[:2] == "/.":
            relpath = relpath[2:]
        elif relpath == ".":
            relpath = ""
        for d in dirnames:
            #Note use of lower to change module names:
            new = os.path.join(py3folder, relpath.lower(), d.lower())
            if not os.path.isdir(new):
                if verbose:
                    print ("Creating directory %s" % new)
                os.mkdir(new)
        for f in filenames:
            if f.startswith("."):
                #Ignore hidden files
                continue
            elif f.endswith("~") or f.endswith(".bak") \
            or f.endswith(".swp"):
                #Ignore backup files
                continue
            elif f.endswith(".pyc") or f.endswith("$py.class"):
                #Ignore compiled python
                continue
            if f.endswith(".py"):
                f_new = f.lower()
            old = os.path.join(py2folder, relpath, f)
            if f.endswith(".py"):
                new = os.path.join(py3folder, relpath.lower(), f.lower())
            else:
                #Do not make non-python filenames lowercase (e.g. DTD files)
                new = os.path.join(py3folder, relpath.lower(), f)
            #The filesystem can (in Linux) record nanoseconds, but
            #when copying only microsecond accuracy is used.
            #See http://bugs.python.org/issue10148
            #Compare modified times down to milliseconds only. In theory
            #might able to use times down to microseconds (10^-6), but
            #that doesn't work on the Windows machine I'm testing on.
            if os.path.isfile(new) \
            and round(os.stat(new).st_mtime*1000) >= \
                round(os.stat(old).st_mtime*1000):
                if verbose: print("Current: %s" % new)
                continue
            #Python, C code, data files, etc - copy with date stamp etc
            shutil.copy2(old, new)
            assert abs(os.stat(old).st_mtime-os.stat(new).st_mtime)<0.0001, \
                   "Modified time not copied! %0.8f vs %0.8f, diff %f" \
                   % (os.stat(old).st_mtime, os.stat(new).st_mtime,
                      abs(os.stat(old).st_mtime-os.stat(new).st_mtime))
            if f.endswith(".py"):
                #Also run 2to3 on it
                to_convert.append(new)
                if verbose: print("Will convert %s" % new)
            else:
                if verbose: print("Updated %s" % new)
    if to_convert:
        print("Have %i python files to convert" % len(to_convert))
        run2to3(to_convert)

            
def main(python2_source, python3_source,
         children=["Bio", "BioSQL", "Tests", "Scripts", "Doc"]):
    #Note want to use different folders for Python 3.1, 3.2, etc
    #since the 2to3 libraries have changed so the conversion
    #may differ slightly.
    print("The 2to3 library will be called automatically now,")
    print("and the converted files cached under %s" % python3_source)
    if not os.path.isdir(python3_source):
        os.mkdir(python3_source)
    global OLD_NAMES
    OLD_NAMES = list(get_module_names(python2_source))

    for child in children:
        print("Processing %s" % child)
        if child in ["Bio", "BioSQL"]:
            #Note we make the path name lower case on Python 3
            do_update(os.path.join(python2_source, child),
                      os.path.join(python3_source, child.lower()))
        else:
            #Don't change the case of the Tests folder etc
            do_update(os.path.join(python2_source, child),
                      os.path.join(python3_source, child))
    print("Python 2to3 processing done.")
              
if __name__ == "__main__":
    python2_source = "."
    python3_source = "build/py%i.%i" % sys.version_info[:2]
    main(python2_source, python3_source)
