#!/bin/bash

# 
# 20160830 switching over to a new compilation scheme
# to eliminate the lex/yacc dependencies for installing
# SCUFF-EM.
#
# The lex/yacc input files (scanner.l and parser.y) will
# continue to be distributed with the code, but the
# derived files generated from them (parser.c, parser.h, scanner.c)
# will also be distributed, and Makefile.am will refer
# only to these derived files, so flex/yacc will not be 
# needed to build SCUFF-EM from its repository.
#
# For developers: Whenever you modify the files 
# scanner.l or parser.y, run this script to update
# the resulting output files.
# 
#

YACC=bison
LEX=flex

${YACC} -d parser.y -o parser.c
${LEX} -o scanner.c scanner.l
