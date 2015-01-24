#!/bin/sh

# Use this script to regenerate all of the auto* files, e.g. when
# you check out a fresh repository.  (Thereafter, as long as
# you configure with the --enable-maintainer-mode during development,
# the files are regenerated automatically as needed whenever you
# change configure.ac, Makefile.am, or similar.)

# automake requires a ChangeLog file to be GNU-ly correct; this will be
# overwritten later by the output of "git log"
touch ChangeLog

# paranoia: sometimes autoreconf doesn't get things right the first time
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force
autoreconf --verbose --install --symlink --force

./configure --enable-maintainer-mode "$@"
