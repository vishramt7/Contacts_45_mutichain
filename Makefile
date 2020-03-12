P=Contacts_45_multichain Contacts_45_multichain_parallel temp temp2 temp_parallel Contacts_CA_test Contacts_AA_test
OBJECTS =
CFLAGS =-g -Wall -std=gnu11 -O3 -ggdb3
LDFLAGS = -L/home/user/XDR_FILES/lib
LDLIBS = -lxdrfile -lm -fopenmp
CC=gcc
$(P):$(OBJECTS)

Contacts_CA : Contacts_CA_modxdr.c
	gcc Contacts_CA_modxdr.c -L/home/user/XDR_FILES/lib -lxdrfile -lm -o Contacts_CA

Contacts_AA : Contacts_AA.c
	gcc Contacts_AA.c -L/home/user/XDR_FILES/lib -lxdrfile -lm -o Contacts_AA

Contacts_45 : Contacts_45.c
	gcc Contacts_45.c -lm -o Contacts_45

Contacts_45_xtc : Contacts_45_xtc.c
	gcc Contacts_45_xtc.c -L/home/user/XDR_FILES/lib -lxdrfile -lm -o Contacts_45_xtc

AA_CA_weights : AA_CA_weights.c
	gcc AA_CA_weights.c -lm -o AA_CA_weights

Contact_map : Contact_map.c
	gcc Contact_map.c -L/home/user/XDR_FILES/lib -lxdrfile -lm -o Contact_map

