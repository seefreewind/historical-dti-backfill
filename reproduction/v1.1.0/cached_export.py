"""I/O adapter only: run unchanged exporter with a bounded read-only SQLite cache."""
import sqlite3,runpy,sys
original=sqlite3.connect
def connect(*a,**k):
 c=original(*a,**k);c.execute('PRAGMA cache_size=-1048576');c.execute('PRAGMA query_only=ON');return c
sqlite3.connect=connect
script=sys.argv.pop(1)
runpy.run_path(script,run_name='__main__')
