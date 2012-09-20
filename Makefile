# Makefile for the 'STORM' program

NAME = PatternMatch
LIBS =  -lpopt -lm  #-lefence
CC = gcc
CPP = g++
DEBUGFLAGS = -Wall -g
CFLAGS = -O2 #-static #-O2
MAIN = main.cpp

$(NAME):	$(MAIN)
	$(CPP) $(DEBUGFLAGS) $(CFLAGS) -o $(NAME) $(MAIN) $(LIBS)

