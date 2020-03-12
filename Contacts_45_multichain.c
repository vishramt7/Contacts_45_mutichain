#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#define CUTOFF 4.5
#define DIFF 4

void DIST (double *x2, double *x1, double *y2, double *y1, double *z2, double *z1, int *dist) {

  if ((sqrt(((*x2 - *x1) * (*x2 - *x1)) + ((*y2 - *y1) * (*y2 - *y1)) + ((*z2 - *z1) * (*z2 - *z1)))) <= CUTOFF) {
    *dist = 1;
  }
  else {
  	*dist = 0;
  }
}

void UPDATE (char *str_point, long int *ATOM_NO, long int *RESIDUE_NO, double *X, double *Y, double *Z){
  char atom_no_array[6] = {}, residue_no_array[5] = {}, x_array[9] = {}, y_array[9] = {}, z_array[9] = {};
  //printf("%s\n", str_point);    

  for (int a = 0; a < strlen (str_point); a++) {

    if ((a >= 6) && (a <= 10)) {
      atom_no_array[a-6] = str_point[a];
    }

    if ((a >= 22) && (a <= 25)) {
      residue_no_array[a-22] = str_point[a];
    }

    if ((a >= 30) && (a <= 37)) {
      x_array[a-30] = str_point[a];
    }

    if ((a >= 38) && (a <= 45)) {
      y_array[a-38] = str_point[a];
    }

    if ((a >= 46) && (a <= 53)) {
      z_array[a-46] = str_point[a];
    }

    if (a >53) {
      break;
    }
  }

  char *end;
  *ATOM_NO = strtol (atom_no_array, &end, 10);
  *RESIDUE_NO = strtol (residue_no_array, &end, 10);
  *X = strtod (x_array, &end);
  *Y = strtod (y_array, &end);
  *Z = strtod (z_array, &end);
  //printf("%ld %ld %f %f %f\n",*ATOM_NO, *RESIDUE_NO, *X, *Y, *Z);
}

void contact_calculation (FILE *pdb_point, int *a, int *b, long int *chA_start, long int *chA_end, long int *chB_start, long int *chB_end) {
	
	long int chA = (*chA_end - *chA_start) + 50; // Here we set the size of arrays for chain A and B, 50 is the buffer for memory
	long int chB = (*chB_end - *chB_start) + 50;

	char strg2[300] = {};
	long int line_number2 = 0, atom_no = 0, residue_no = 0;
	double atom_X = 0.0, atom_Y = 0.0, atom_Z = 0.0;
	
	int length; 
	if (*a == *b) {
		length = chA;
  }
	else {
		length = chA + chB;
  }

  //printf("%d %d %ld %ld %ld %ld %ld\n", *a, *b, *chA_start, *chA_end, *chB_start, *chB_end, length);

  int *atomno_arr = calloc (length, sizeof (int));
  int *residueno_arr = calloc (length, sizeof ( int));
  int *MAP = calloc (length, sizeof ( int)); 
  int *MAP2 = calloc (length, sizeof ( int));
  double *X_arr = calloc (length, sizeof (double));
  double *Y_arr = calloc (length, sizeof (double));
  double *Z_arr = calloc (length, sizeof (double));

  if ((!residueno_arr) || (!atomno_arr) || (!MAP) || (!MAP2) || (!X_arr) || (!Y_arr) || (!Z_arr)) {
    printf ("calloc error\n");
  }

  rewind (pdb_point);
  long int index = 0, k = 0;
  MAP[k] = index;
  long int x_min = 0, x_max = 0, y_min = 0, y_max = 0;
  //printf("%d %d %d %d %d %d %d\n", *a, *b, *chA_start, *chA_end, *chB_start, *chB_end, length);
  while (fscanf(pdb_point, "%299[^\n]\n", strg2) == 1) {

    line_number2 ++;

    if (((line_number2 >= *chA_start) && (line_number2 <= *chA_end)) || ((line_number2 >= *chB_start) && (line_number2 <= *chB_end))) {

      //printf("%s\n", strg2);
      if (line_number2 > *chB_end) {
        break;
      }

      if (line_number2 == *chA_start) {
        x_min = index;
      }

      if (line_number2 == *chA_end) {
        x_max = index;
      }

      if (line_number2 == *chB_start) {
        y_min = index;
      }

      if (line_number2 == *chB_end) {
        y_max = index;
      }

      UPDATE (strg2, &atom_no, &residue_no, &atom_X, &atom_Y, &atom_Z);

      atomno_arr[index] = atom_no;
      residueno_arr[index] = residue_no;
      X_arr[index] = atom_X;
      Y_arr[index] = atom_Y;
      Z_arr[index] = atom_Z;

      if (residueno_arr[MAP[k]] != residueno_arr[index]) {
        k++;
        MAP[k] = index;
      }
      MAP2[index] = k;
      index ++;
    }
  }
	
  // This part is allocating memory to all atom contacts
  int **AA_array = calloc ( length, sizeof ( int *));

  if (!AA_array) {
    printf ("mem not allocated to rows\n");
  }

  for ( long int z = 0; z < length; z++) {
    AA_array[z] = calloc (length, sizeof ( int));

    if (!AA_array[z]) {
      printf ("mem not allocated to columns\n");        
    }
  }

  // This part allocates memory to c-alpha contacts
  int **CA_array = calloc ( k + 1, sizeof ( int *));

  if (!CA_array) {
    printf ("mem not allocated to rows\n");
  }

  for (long int z1 = 0; z1 < (k + 1); z1++) {
    CA_array[z1] = calloc ( k + 1, sizeof ( int));

    if (!CA_array[z1]) {
      printf ("mem not allocated to columns\n");        
    }
  }
  // Memory allocation ends here 

  //printf("%ld %ld %ld %ld %ld %ld\n", x_min, x_max, y_min, y_max, chA, chB);
  for (long int x = x_min; x <=  x_max; x++) {
    for (long int y = y_min; y <= y_max; y++) {
      //AA_array[x][y] = 0;
      DIST (&X_arr[y], &X_arr[x], &Y_arr[y], &Y_arr[x], &Z_arr[y], &Z_arr[x], &AA_array[x][y]);
      //printf("%ld\n",AA_array[x][y] );
    }
  }
//  printf("%ld %ld %ld %ld %ld %ld %ld\n", x_min, x_max, y_min, y_max, chA, chB, AA_array[x_max][y_max]);

  char *AAname;
  char *CAname;
  asprintf (&AAname,"AA_ch%d_ch%d.cont",*a + 1, *b + 1);
  asprintf (&CAname,"CA_ch%d_ch%d.cont",*a + 1, *b + 1);
  FILE *AA_file = fopen (AAname,"w");
  FILE *CA_file = fopen (CAname,"w");

  for (long int x1 = x_min; x1 <= x_max; x1++) {
    for (long int y1 = y_min; y1 <= y_max; y1++) {
      
      if (AA_array[x1][y1] == 1) {

        CA_array[MAP2[x1]][MAP2[y1]] ++;

        if (*a == *b) {
          if ( residueno_arr[y1] >= residueno_arr[x1] + DIFF) { 
            fprintf(AA_file,"%d %d %d %d \n", *a + 1, atomno_arr[x1], *b + 1, atomno_arr[y1]); // This can be used to print all atom contacts
          } 
        }
        else {
          fprintf(AA_file,"%d %d %d %d \n", *a + 1, atomno_arr[x1], *b + 1 , atomno_arr[y1]); // This can be used to print all atom contacts
        }
      }
    }
  }

  
  for (long int x2 = 0; x2 <= k ; x2++) {
    for (long int y2 = x2 ; y2 <= k; y2++) {
      if (CA_array[x2][y2] > 0) {
      	if (*a == *b) {
      		if (residueno_arr[MAP[y2]] >= residueno_arr[MAP[x2]] + DIFF) {
      			fprintf(CA_file,"%d %d %d %d %d \n", *a + 1, residueno_arr[MAP[x2]], *b + 1, residueno_arr[MAP[y2]], CA_array[x2][y2]); // This can be used to print C-alpha contacts		
      		}
      	}
      	else {
      	 fprintf(CA_file,"%d %d %d %d %d \n", *a + 1, residueno_arr[MAP[x2]], *b + 1, residueno_arr[MAP[y2]], CA_array[x2][y2]); // This can be used to print C-alpha contacts
      	} 
      } 
    }
  }

  for (long int z2 = 0; z2 < length; z2++) { // This part is to free the 2d array
    free (AA_array[z2]);
  }

  for (long int z3 = 0; z3 < (k + 1); z3++) { // This part is to free the 2d array
    free (CA_array[z3]);
  }  

	free (AA_array);
	free (CA_array);
	free (atomno_arr);
	free (residueno_arr);	
	free (MAP);
  free (MAP2);
	free (X_arr);
	free (Y_arr);
	free (Z_arr);
  fclose (AA_file);
  fclose (CA_file);
  free (AAname);
  free (CAname);
}

int main(int argc, char const *argv[]){

  if (argc < 2) {
  	fprintf (stdout,"Minimum 2 arguments needed\n");
  }

	FILE *pdb = fopen (argv[1],"r");

	if (pdb == NULL) {
  	fprintf (stdout,"CANNOT OPEN %s\n",argv[1]);
	}

  char column1[7] = {}, strg[300] = {};
  long int line_number = 0;
  long int *chainid_start_array = calloc (50, sizeof (long int)); // Initial size of 50 for these arrays, should be enough
  long int *chainid_array = 	calloc (50, sizeof (long int));
  int chainid_start = 0;
  int chainid = 0;

  if ((!chainid_array) || (!chainid_array)) {
  	printf("calloc error\n");
  	exit (1);
  }

  while (fscanf(pdb, "%299[^\n]\n", strg) == 1){

  	line_number ++;

  	for (int k = 0; k < 5; k++) //This loop copies the first 6 columns of each line 
  		column1[k] = strg[k];	

  	if ((chainid_start == chainid) && (strstr (column1,"ATOM"))) {
  		//printf("%s This is ATOM , line_number %ld\n", column1, line_number);
  		chainid_start_array[chainid_start] = line_number;
  		chainid_start ++;
  	}

  	if (strstr (column1,"TER")) {
  		chainid_array[chainid] = line_number - 1;
  		chainid ++;
  	}
    
    if ((chainid > 0) && (chainid % 40 == 0)) {
  		chainid_array = realloc (chainid_array, (chainid + 40) * sizeof(long int));
  		chainid_start_array = realloc (chainid_start_array, (chainid_start + 40) * sizeof(long int));
  		//printf(" More than %d chains in this file\n", chainid);
  	}
    
  	if ((feof(pdb)) || (strstr(column1,"END"))) {

  		for (int i = 0; i < chainid; i++) {
  			for (int j = i; j < chainid; j++) {
  				contact_calculation (pdb, &i, &j, &chainid_start_array[i], &chainid_array[i], &chainid_start_array[j], &chainid_array[j]);
  				//printf("%d %d %d %ld %ld %ld %ld\n", chainid, i, j, chainid_start_array[i], chainid_array[i], chainid_start_array[j], chainid_array[j]);
  			}
  		}
  	}
  }
  printf("Total no. of chains are %d, the CUTOFF is %0.5lf A and diff is %d\n", chainid, CUTOFF, DIFF);
  free (chainid_array);
  free (chainid_start_array);
	fclose (pdb);
}
