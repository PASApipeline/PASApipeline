# PASApipeline
PASA software

See wiki tab for documentation.


remaining to include in documentation:

1.  Defaults to SQLite. Users can override this by setting the DBI_DRIVER environment variable (see the DBI perldoc) to 'mysql' (if defaulting to SQLite is too bold for the next release, the default could be set the other way around, and SQLite chosen with DBI_DRIVER=SQLite).

2. The user instead sets a DATABASE parameter to either the name of the MySQL database, or the absolute pathname of the SQLite database.

3. When using SQLite, the $PASAHOME/pasa_conf/conf.txt file is optional (defaulting to ${PASAHOME}/pasa_conf/pasa.TEMPLATE for hooks) when SQLite is used. If it exists, it still will be used, and users can still override the pathname with PASACONF environment variable.
