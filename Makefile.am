
# Removed mpi for now
SUBDIRS = src util tests

EXTRA_DIST = VERSION LICENSE NOTICE config.site

CLEANFILES = so_locations

ACLOCAL_AMFLAGS = -I config

install-data-local:
	$(mkinstalldirs) $(DESTDIR) $(GMAPDB)

dist-hook:
#	svn log -v --xml | ./svncl.pl > ChangeLog

