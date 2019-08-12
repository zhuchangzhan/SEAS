#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
"""
Background and base for making database queries.
"""

import os
import time
import urllib

import logging
logger = logging.getLogger(__name__)


class Database(object):
    """
    Database mixin provides a standard process for connecting to a database
    """
    CREDENTIALS_FILENAME = '~/.config/tsig/tsig-dbinfo'
    CREDENTIALS_URL = 'http://tessellate.mit.edu/tsig/tsig-dbinfo'

    def __init__(self, dbinfo_loc=None,
                 dbhost=None, dbname=None, dbport=None,
                 dbuser=None, dbpass=None, dbtable=None):
        # start with database parameters from file or url
        self.dbinfo = self.read_dbinfo(dbinfo_loc)
        # override with credentials and other connection parameters
        if dbhost is not None:
            self.dbinfo['dbhost'] = dbhost
        if dbname is not None:
            self.dbinfo['dbname'] = dbname
        if dbport is not None:
            self.dbinfo['dbport'] = dbport
        if dbuser is not None:
            self.dbinfo['dbuser'] = dbuser
        if dbpass is not None:
            self.dbinfo['dbpass'] = dbpass
        if dbtable is not None:
            self.dbinfo['dbtable'] = dbtable

        if self.dbinfo.get('dbport') is None:
            self.dbinfo['dbport'] = 5432

        for k in self.dbinfo:
            v = self.dbinfo[k]
            if k == 'dbpass':
                v = "XXXX"
            logger.debug("%s: %s" % (k, v))

    @classmethod
    def read_dbinfo(cls, path=None):
        """
        Try to read from the specified location, which could be a file or a
        ULR.  If the specified location does not work, fail fast.

        If nothing is specified, then fallback to first trying the default
        file location, then try the web server and cache result.
        """
        default_filename = cls.CREDENTIALS_FILENAME
        default_url = cls.CREDENTIALS_URL
        if not path or '://' in path:
            url = path or default_url
            filename = os.path.abspath(os.path.expanduser(default_filename))
            try:
                if not os.path.isfile(filename):
                    try:
                        os.makedirs(os.path.dirname(filename))
                    except os.error:
                        pass
                    testfile = urllib.URLopener()
                    testfile.retrieve(url, filename)
            except Exception as err:
                logger.error("Download of credentials failed: "
                             "%s (%s)" % (str(err), url))
        else:
            filename = os.path.abspath(os.path.expanduser(path))

        # now read the parameters
        dbinfo = dict()
        try:
            with open(filename, 'r') as f:
                for line in f:
                    name, value = line.split('=')
                    dbinfo[name.strip()] = value.strip()
        except IOError:
            logger.error("No %s database credentials at %s" %
                         (cls.__name__, filename))
        return dbinfo

    @staticmethod
    def check_dbinfo(dbinfo):
        missing = []
        for k in ['dbhost', 'dbname', 'dbuser', 'dbpass', 'dbtable']:
            if k not in dbinfo:
                missing.append(k)
        if missing:
            raise Exception("Missing database parameters: %s" %
                            ','.join(missing))

    def select(self, *args, **where):
        """
        Basic SQL select function. fields is a list of fields to select
        if non provided, then all fields are selected.

        where is a dictionary of field names or selection conditionals
        """
        limit = where.pop('limit', '')
        if limit:
            limit = "LIMIT %s" % str(limit)

        order = where.pop('order', '')
        if order:
            order = "ORDER BY %s" % str(order)

        fields = "%s" % ", ".join(args or ['*'])
        vals = []
        whr = self._where(vals, **where)
        sql = "SELECT %s FROM %s %s %s %s" % (fields, self.table, whr, order, limit)
        return self._do_query(sql, vals)

    def insert(self, **values):
        """Insert one row into this table"""
        fields = []
        vals = []
        valt = []

        for (field, value) in values.items():
            fields.append(field)
            if isinstance(value, tuple) and len(value) == 2:
                valt.append("POINT(%s, %s)")
                vals += list(value)
            else:
                valt.append("%s")
                vals.append(value)

        sql = "INSERT INTO %s (%s) VALUES (%s)" % \
            (self.table, ", ".join(fields), ", ".join(valt))
        return self._do_query(sql, vals)

    def create_table(self):
        fields = []
        for (name, kind) in self.schema:
            fields.append("%s %s" % (name, kind.upper()))
        sql = "CREATE TABLE %s (%s)" % (self.table, ",\n".join(fields))
        return self._do_query(sql)

    def drop_table(self):
        sql = "DROP TABLE %s" % (self.table)
        return self._do_query(sql)

    def _where(self, vals, **where):
        """Generates a where clause for sql queries"""
        if not where:
            return ""
        ret = []
        for (name, value) in where.items():
            op = '='
            if '__' in name:
                (name, op) = name.rsplit('__', 1)

            if '%' in name:
                ret.append(name)
                vals.append(value)
                continue

            # Special within clause for positional selections
            special = "where_" + op
            if hasattr(self, special):
                ret.append(getattr(self, special)(name, value, vals))
                continue

            # Add more operators here as needed.
            op = {'lt': '<', 'gt': '>'}.get(op, op)

            ret.append("%s %s %%s" % (name, op))
            vals.append(value)
        return "WHERE " + " AND ".join(ret)

    @staticmethod
    def where_within(name, value, vals):
        """Select within a geometric circle radius"""
        if len(value) != 3:
            raise IOError("When selecting within a radius, three values must be provided")
        vals += list(value)
        return "%s <@ circle '((%%s,%%s), %%s)'" % name

    @property
    def connection(self):
        """Return the database connect"""
        try:
            import psycopg2
        except ImportError as e:
            logger.error("PostgreSQL python module psycopg2 is required")
            logger.debug("psycopg2 import failed: %s" % e)
            return None

        # Decide if we want to keep the connection live of not.
        return psycopg2.connect(
            host=self.dbinfo['dbhost'], database=self.dbinfo['dbname'],
            port=self.dbinfo['dbport'],
            user=self.dbinfo['dbuser'], password=self.dbinfo['dbpass'],
            connect_timeout=5)

    def _do_query(self, sqlcmd, *args, **kw):
        logger.debug("query: %s" % sqlcmd)
        t0 = time.time()

        conn = self.connection
        cur = conn.cursor()
        t1 = time.time()
        logger.debug("query start")
        cur.execute(sqlcmd, *args, **kw)
        t2 = time.time()
        logger.debug("query complete")
        if 'select' in sqlcmd.lower():
            result = QueryResults(cur)
        else:
            result = conn.commit()

        t3 = time.time()
        logger.debug("data retrieved")
        t4 = time.time()
        logger.debug("query: setup=%.3f execute=%.3f fetch=%.3f close=%.3f"
                     % (t1 - t0, t2 - t1, t3 - t2, t4 - t3))
        return result

    def _do_query_raw(self, sqlcmd):
        logger.debug("query: %s" % sqlcmd)
        t0 = time.time()
        conn = self.connection
        cur = conn.cursor()
        t1 = time.time()
        logger.debug("query start")
        cur.execute(sqlcmd)
        t2 = time.time()
        logger.debug("query complete")
        result = cur.fetchall()
        column_names = [desc[0] for desc in cur.description]
        t3 = time.time()
        logger.debug("data retrieved")
        cur.close()
        cur = None
        conn.close()
        conn = None
        t4 = time.time()
        logger.info("query: setup=%.3f execute=%.3f fetch=%.3f close=%.3f"
                    % (t1 - t0, t2 - t1, t3 - t2, t4 - t3))
        return result, column_names



class QueryResults(object):
    """Returned when we make a query"""
    def __init__(self, cur):
        self.cur = cur

    def __iter__(self):
        for row in self.cur.fetchall():
            yield row
        self.close()

    def as_dict(self):
        """Return the query result as a dictionary per row (generator)"""
        cols = self.columns
        for row in self:
            yield dict(zip(cols, row))

    @property
    def columns(self):
        return [desc[0] for desc in self.cur.description]

    def close(self):
        self.cur.close()
        self.cur.connection.close()

    def __len__(self):
        return self.cur.rowcount
