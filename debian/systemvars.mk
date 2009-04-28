
# for SHELL, do not change the type of shell - only use Bourne or BASH
SHELL = /bin/sh

# System dependent commands (NB: the first two are the most platform dependent)

INSTALL = install -p
RANLIB = ranlib

RM = /bin/rm
CP = /bin/cp
CHMOD = /bin/chmod
MKDIR = /bin/mkdir
TCLSH = /usr/bin/tclsh

# Compiler dependent variables

CC = gcc
CXX = g++
CSTATICFLAGS =
CXXSTATICFLAGS =

DEPENDFLAGS = -MM
