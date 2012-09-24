# Copyright 2001 by Katharine Lindner.  All rights reserved.
# Copyright 2006 by PeterC.  All rights reserved.
# Copyright 2007 by Michiel de Hoon.  All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
"""Parser for files from NCBI's Gene Expression Omnibus (GEO).

http://www.ncbi.nlm.nih.gov/geo/
"""

import Record


def _read_key_value(line):
    words = line[1:].split("=", 1)
    try:
        key, value = words
        value = value.strip()
    except ValueError:
        key = words[0]
        value = ""
    key = key.strip()
    return key, value


def parse(handle):
    rec = None
    for line in handle:
        line = line.strip('\n').strip('\r')
        if not line: continue # Ignore empty lines
        c = line[0]
        if c=='^':
            if rec:
                yield rec
            rec = Record.Record()
            rec.entity_type, rec.entity_id = _read_key_value(line)
        elif c=='!':
            if line in ('!Sample_table_begin',
                        '!Sample_table_end',
                        '!Platform_table_begin',
                        '!Platform_table_end'):
                continue
            key, value = _read_key_value(line)
            if key in rec.entity_attributes:
                if type(rec.entity_attributes[key])==list:
                    rec.entity_attributes[key].append(value)
                else:
                    existing = rec.entity_attributes[key]
                    rec.entity_attributes[key] = [existing, value]
            else:
                rec.entity_attributes[key] = value
        elif c=='#':
            key, value = _read_key_value(line)
            assert key not in rec.col_defs
            rec.col_defs[key] = value
        else:
            row = line.split("\t")
            rec.table_rows.append(row)
    if rec:
        yield rec
