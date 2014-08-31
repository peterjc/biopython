# Copyright 2014 by Peter Cock.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Provides code to access the REST-style KEGG online API.

This module aims to make the KEGG online REST-style API easier to use. See:
http://www.kegg.jp/kegg/rest/keggapi.html

The KEGG REST-style API provides simple access to a range of KEGG databases.
This works using simple URLs (which this module will construct for you),
with any errors indicated via HTTP error levels.

The functionality is somewhat similar to Biopython's Bio.TogoWS and Bio.Entrez
modules.

Currently KEGG does not provide any usage guidelines (unlike the NCBI whose
requirements are reasonably clear). To avoid risking overloading the service,
Biopython will only allow three calls per second.
"""

from __future__ import print_function

import time
from Bio._py3k import _binary_to_string_handle, _as_bytes

#Importing these functions with leading underscore as not intended for reuse
from Bio._py3k import urlopen as _urlopen
from Bio._py3k import quote as _quote

#Constant
_BASE_URL = "http://rest.kegg.jp"

#Caches:
_db_names = set([
    'pathway', 'path',
    'brite', 'br',
    'module', 'md',
    'orthology', 'ko',
    'genome',
    'genomes', 'gn',
    'genes',
    'ligand', 'ligand',
    'compound', 'cpd',
    'compound_ja', 'cpd_ja',
    'glycan', 'gl',
    'reaction', 'rn',
    'rpair', 'rp',
    'rclass', 'rc',
    'enzyme', 'ec',
    'disease', 'ds',
    'disease_j', 'ds_ja',
    'drug', 'dr',
    'drug_ja', 'dr_ja',
    'dgroup', 'dg',
    'dgroup_ja', 'dg_ja',
    'environ', 'ev',
    'environ_ja', 'ev_ja',
    #Auxillirary databases:
    'dgenes', 'egenes', 'mgenes', 'egenome', 'mgenome',
])


def kegg_info(db):
    """KEGG info - Displays the current statistics of a given database.

    db - database or organism (string)

    The argument db can be a KEGG database name (e.g. 'pathway' or its
    official abbreviation, 'path'), or a KEGG organism code or T number
    (e.g. 'hsa' or 'T01001' for human).

    A valid list of organism codes and their T numbers can be obtained
    via kegg_info('organism') or http://rest.kegg.jp/list/organism

    Returns a string.
    """
    #TODO - chache and validate the organism code / T numbers?
    #TODO - can we parse the somewhat formatted output?
    pass

def kegg_list_database(db, org=None):
    """KEGG list - Entry list for database.

    db - database or organism (string)
    org - optional organism (string), see below.

    For the pathway and module databases the optional organism can be
    used to restrict the results.

    Returns a list of tuples.
    """
    pass

def kegg_list_dbentries(dbentries):
    """KEGG list - Entry list for specified database entries.

    dbentries - Up to 100 identifies (list of strings).

    Returns a list of tuples.

    Example
    -------
    Using kegg_list_dbentries(['hsa:10458', 'ece:Z5100']) would give
    two lines via http://rest.kegg.jp/list/hsa:10458+ece:Z5100
    """
    pass

def kegg_find(db, query, option=None):
    """KEGG find - Data search.

    Finds entries with matching query keywords or other query data in
    a given database.

    db - database or organism (string)
    query - search terms (string)
    option - search option (string), see below.

    For the compound and drug database, set option to the string 'formula',
    'exact_mass' or 'mol_weight' to search on that field only. The
    chemical formula search is a partial match irrespective of the order
    of atoms given. The exact mass (or molecular weight) is checked by
    rounding off to the same decimal place as the query data. A range of
    values may also be specified with the minus(-) sign.

    Returns a list of tuples.
    """
    #Quote encoding!
    pass

def kegg_get(dbentries, option=None):
    """KEGG get - Data retrieval.

    dbentries - Identifiers (single string, or list of strings), see below.
    option - One of "aaseq", "ntseq", "mol", "kcf", "image", "kgml" (string)

    The input is limited up to 10 entries.
    The input is limited to one pathway entry with the image or kgml option.
    The input is limited to one compound/glycan/drug entry with the image option.

    Returns a handle.
    """
    if isinstance(dbentries, list):
        dbentries = "+".join(dbentries)        
    url = _BASE_URL + "/get/%s" % _quote(dbentries)
    binary = False
    if option:
        if option not in ("aaseq", "ntseq", "mol", "kcf", "image", "kgml"):
            raise ValueError("The option argument should be 'aaseq', 'ntseq', "
                             "'mol', 'kcf', 'image', or 'kgml', "
                             "and not %r" % option)
        url += "/" + option
        if option in ("image", ):
            binary = True
    return _open(url, binary=binary)


def kegg_convert(data, target_db, source_db, option=None):
    """KEGG conv – convert KEGG identifiers to/from outside identifiers

    target_db - 
    source_db_or_dbentries - 
    option - Can be "turtle" or "n-triple" (string).
    """
    #TODO functions? 
    pass

def _open(url, post=None, binary=False):
    """Helper function to build the URL and open a handle to it (PRIVATE).

    Open a handle to KEGG, will raise an IOError if it encounters an error.

    In the absense of clear guidelines, this function enforces a limit of
    "up to three queries per second" to avoid abusing the KEGG servers.
    """
    delay = 0.333333333  # one third of a second
    current = time.time()
    wait = _open.previous + delay - current
    if wait > 0:
        time.sleep(wait)
        _open.previous = current + wait
    else:
        _open.previous = current

    #print(url)
    if post:
        handle = _urlopen(url, _as_bytes(post))
    else:
        handle = _urlopen(url)

    #We now trust KEGG to have set an HTTP error code.

    if binary:
        return handle
    else:
        return _binary_to_string_handle(handle)

_open.previous = 0
