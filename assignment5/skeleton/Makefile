#      _________   _____________________  ____  ______
#     / ____/   | / ___/_  __/ ____/ __ \/ __ \/ ____/
#    / /_  / /| | \__ \ / / / /   / / / / / / / __/
#   / __/ / ___ |___/ // / / /___/ /_/ / /_/ / /___
#  /_/   /_/  |_/____//_/  \____/\____/_____/_____/
#
#  http://www.inf.ethz.ch/personal/markusp/teaching/
#  How to Write Fast Numerical Code 263-2300 - ETH Zurich
#  Copyright (C) 2015  Alen Stojanov      (astojanov@inf.ethz.ch)
#                      Daniele Spampinato (daniele.spampinato@inf.ethz.ch)
#                      Singh Gagandeep    (gsingh@inf.ethz.ch)
#	                   Markus Pueschel    (pueschel@inf.ethz.ch)
#  Copyright (C) 2013  Georg Ofenbeck     (ofenbeck@inf.ethz.ch)
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program. If not, see http://www.gnu.org/licenses/.

# ==================================================================================== #
# = Folder structure / tests
# ==================================================================================== #

EXEC   = compute
INCDIR = include
LIBDIR = lib
BINDIR = bin
SRCDIR = src

# ==================================================================================== #
# = Compiler settings
# ==================================================================================== #

AFLAGS = -O3 -m64 -march=corei7 # remove avx to compare to sse4 compiler optimizations
DFLAGS = -g -Wall
LIBS   = -lm

# ==================================================================================== #
# = Default targets
# ==================================================================================== #

TARGET = $(BINDIR)/$(EXEC)

MACROS = 
		  
CFLAGS += -I$(INCDIR)/ $(MACROS)

LIBOBJ = $(BINDIR)/complex.o \
		 $(BINDIR)/validate.o
		 
SRCOBJ = $(BINDIR)/comp.o \
		 $(BINDIR)/main.o  
 
$(BINDIR)/%.o : $(LIBDIR)/%.c $(INCDIR)/%.h
	$(CC) -c $(CFLAGS) $< -o $@

$(BINDIR)/%.o : $(SRCDIR)/%.c $(SRCDIR)/%.h
	$(CC) -c $(AFLAGS) $(DFLAGS) $(CFLAGS) $< -o $@

$(BINDIR)/%.o : $(SRCDIR)/%.c $(INCDIR)/%.h
	$(CC) -c $(AFLAGS) $(DFLAGS) $(CFLAGS) $< -o $@


all: $(LIBOBJ) $(SRCOBJ)
	$(CC) $(LIBOBJ) $(SRCOBJ) $(LIBS) -o $(TARGET)
	
clean:
	rm -rf $(LIBOBJ)
	rm -rf $(SRCOBJ)
	rm -rf $(TARGET)
	
