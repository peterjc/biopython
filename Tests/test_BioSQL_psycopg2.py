#!/usr/bin/env python
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

"""Run BioSQL tests using PostgreSQL"""
from Bio import MissingExternalDependencyError
from BioSQL import BioSeqDatabase

from common_BioSQL import *

##################################
# Start of user-editable section #
##################################

# I ran this by hand on a buildbot machine. Note it is an unfortunate
# limitation of _do_db_create() that currently we need permission to
# create a database (the code uses DROP DATABASE; CREATE DATABASE)
# and I don't see an easy portable way to do "DROP TABLE *;" instead.
#
# $ psql -U postgres -c "CREATE USER biosql WITH PASSWORD 'testing' CREATEDB"
# CREATE ROLE
# $ psql -U postgres -c "CREATE DATABASE biosql_test"
# CREATE DATABASE
# $ psql -U postgres -c "GRANT ALL PRIVILEGES ON DATABASE biosql_test TO biosql"
# GRANT
# $ psql -U postgres -c ""
#
# And added this to the pg_hba.conf file:
#
# # Access from local host for running BioSQL tests
# local   "biosql_testing" biosql                            md5
# host    "biosql_testing" biosql      127.0.0.1/32          md5
# host    "biosql_testing" biosql      ::1/128               md5

# Constants for the database driver
DBHOST = 'localhost'
DBUSER = 'biosql'  # not just 'postgres'
DBPASSWD = 'testing'  # not just ''
TESTDB = 'biosql_test'

################################
# End of user-editable section #
################################

DBDRIVER = 'psycopg2'
DBTYPE = 'pg'

# This will abort if driver not installed etc:
check_config(DBDRIVER, DBTYPE, DBHOST, DBUSER, DBPASSWD, TESTDB)

# Some of the unit tests don't create their own database,
# so just in case there is no database already:
create_database()

if __name__ == "__main__":
    # Run the test cases
    runner = unittest.TextTestRunner(verbosity=2)
    unittest.main(testRunner=runner)
