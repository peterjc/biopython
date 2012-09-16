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

#Third idea - brute force hackery?

#TODO - get list of module names either from setup.py or
#walking the file system.

NAMES = [] #Updated later!

def hack_simple_import(line, core):
    global NAMES
    assert core.startswith("import "), core
    old = core[7:].strip()
    if old in NAMES:
        #print("%r --> %r" % (line, line.replace(old, old.lower())))
        return line.replace(old, old.lower())
    else:
        return line

def hack_from_import(line, core):
    global NAMES
    old = line
    assert core.startswith("from ") and " import " in core, core
    if "#" in line:
        i = line.find("#")
        comment = " " + line[i:].lstrip()
        line = line[:i]
    else:
        comment = ""
    from_part = core[5:].split(" import ", 1)[0]
    if from_part not in NAMES:
        return line
    if " as " in line:
        i = line.find(" as ")
        as_part = line[i:]
        line = line[:i]
    else:
        as_part = ""
    bits = []
    for x in line.split(" import ", 1)[1].split(","):
        x = x.strip()
        if from_part + "." + x in NAMES:
            x = x.lower()
        bits.append(x)
    new = line.split(" import ", 1)[0].lower() + " import " + ", ".join(bits) + as_part + comment
    new = new.rstrip() + "\n"
    #print("OLD: %r\nNEW: %r\n" % (old, new))
    return new

def hack_import(module, line):
    core = line.strip()
    if "#" in core:
        core = core[:core.find("#")].strip()
    if core.startswith("from "):
        if core.endswith("\\") or core.endswith(","):
            raise NotImplementedError("Continued line %r" % line)
        return hack_from_import(line, core)
    elif line.lstrip().startswith("import "):
        if core.endswith("\\") or core.endswith(","):
            raise NotImplementedError("Continued line %r" % line)
        return hack_simple_import(line, core)
    else:
        return line

def hack_file_import_lines(f):
    #Evil hack to change our import lines to use lower case...
    if f.endswith("/__init__.py"):
        m = f[:-12].split(os.path.sep)
    elif f.endswith(".py"):
        m = f[-3:].split(os.path.sep)
    else:
        assert False, f
    h = open(f, "r")
    lines = list(h)
    h.close()
    h = open(f, "w")
    in_triple_quote = False
    for line in lines:
        assert line.count('"""') <= 2, line
        if line.count('"""') == 1:
            in_triple_quote = not in_triple_quote
            h.write(line)
        elif in_triple_quote:
            #Probably a docstring...
            if line.lstrip().startswith(">>> "):
                i = line.find(">>> ")
                h.write(line[:i] + ">>> " + hack_import(m, line[i+4:]))
            else:
                h.write(line)
        else:
            #Should be normal code...
            h.write(hack_import(m, line))
    h.close()
    #TODO - Work out why this fails (commented out to allow
    #all later files to be processed, to spot other issues):
    #if in_triple_quote:
    #    raise ValueError("Triple quote confused in %s" % f)

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
    global NAMES #Used in our import hack
    NAMES = list(get_module_names(python2_source))
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
