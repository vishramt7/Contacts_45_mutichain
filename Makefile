P=Contacts_45_multichain Contacts_45_multichain_parallel
OBJECTS =
CFLAGS =-g -Wall -std=gnu11 -O3 -ggdb3
LDLIBS = -lm -fopenmp
CC=gcc
$(P):$(OBJECTS)
