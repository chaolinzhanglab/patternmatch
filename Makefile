# Makefile for the 'STORM' program

NAME = PatternMatch
BASE = $(CREAD)
LIBDIR = $(BASE)/lib
INCLUDEDIR = $(BASE)/include
BINDIR = $(BASE)/bin
LIBS =  -lpopt -lm  #-lefence
CC = gcc
CPP = g++
DEBUGFLAGS = -Wall -g
CFLAGS = -static #-O2
MAIN = main.cpp

$(NAME):	$(MAIN)
	$(CPP) -static $(DEBUGFLAGS) $(CFLAGS) -o $(NAME) $(MAIN) $(LIBS) -L$(LIBDIR) -I$(INCLUDEDIR)

install:	$(NAME)
	cp $(NAME) $(BINDIR)
.PHONY: install

clean:
	-rm -f $(NAME)
.PHONY: clean
