/* Monte-Carlo simulation of a CO/Pt(111) surface system to replicate observed dynamics.	*/
/* Written by Arthur Gordon-Wright, ajgw20@bath.ac.uk.						*/
/* 												*/
/* Version number 0.0.2, date 20/11/09								*/
/* Update notes: Added initialisation, data structures.						*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include <ctype.h>*/
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <unistd.h>
#include <sys/stat.h>

#define SQRT2 1.414213562373095 /* Square root of 2, to save repeated calculation */
#define SQRT3 1.732050807568877 /* Square root of 3, " */
#define UTT 210			/* Top-top diffusion barrier */
#define USAGE "Usage: %s dim n_co outname\n" /*[a A b B c C]\n"*/
/*#define OPTIONS "a, b, c.\n"*/

/* Variable declarations */
char name[30], outname[40], outmode='w';/* a[10]="A", b[10]="B", c[10]="C";*/
int dim, n_sites, n_co, *occupations, *adatoms;

/* Function declarations */
int args ( int , char ** );
void error_close ( int );
void init_random ( );
int init_adatoms ( );

int main ( int argc, char *argv[] )
{
	int err_code,i,j;

	if( argc==1 || strcmp(argv[1],"help")==0 ) { /* Test if help has been selected or no arguments set. */
		printf(USAGE, argv[0]);
		return 3;
	}

	init_random();	/* Set random seed. */

	if( (err_code=args( argc, argv )) ) { /* Get and test command-line arguments, exit on failure. */
		error_close(err_code);
	}

	if( ( err_code=init_adatoms() ) ) {	/* Initialise adatom positions, exit on failure. */
		error_close(err_code);
	}

	printf("name=%s, dim=%d, n_co=%d, outname=%s\n", name, dim, n_co, outname );
	/*printf("a=%s, b=%s, c=%s\n", a, b, c );*/

	printf("\nOccupations.\n");
	for(i=0;i<n_sites;i++) {
		printf("%d ",occupations[i]);
		if( !((i+1)&(dim-1)) ) printf("\n");
	}

	printf("\nCoordinates.\n");
	for(i=0;i<n_co;i++)
		printf("%d ",adatoms[i]);
	printf("\n\n");

	return 0;
}

int args ( int argc, char *argv[] )	/* Processes command-line arguments and sets global variables accordingly. */
{					/* Returns 0 if successful. */
	int i;
	FILE *testfile;
	struct stat buf;

	/* Extract program name. */
	sprintf(name,"%s",argv[0]);

	/* Test for minimum valid number of arguments. */
	if( argc<4 ) {
		printf("Incorrect number of arguments.\n"USAGE, argv[0]);
		return 2;
	}
	
	/* Extract lattice dimension, number of sites and number of CO adatoms. */
	dim=atoi(argv[1]);
	n_sites=dim*dim;
	n_co=atoi(argv[2]);

	/* Test for successful extraction. */
	if( !dim || !n_co ) {
		printf("Invalid dimension or number of adatoms.  Please specify integer values.\n");
		return 2;
	}

	/* Check that specified lattice dimension is a power of 2, to allow bitwise modulus function. */
	if ( (dim-1) & dim ) {
		printf("For program efficiency, please specify a lattice dimension that is an integer\npower of two.\n");
		return 2;
	}

	/* Check for valid specified coverage (i.e. number of sites is >= number of adatoms. */
	if( (dim*dim)<n_co ) {
		printf("Specified coverage is greater than 100%!\n");
		return 2;
	}

	/* Extract output file name. */
	sprintf(name,"%s",argv[0]);
	sprintf(outname,"%s",argv[3]);

	/* Check if output file is a directory. */
	if( !stat(outname,&buf) )
		if( S_ISDIR(buf.st_mode) ) {
			printf("Specified output file is a directory.\n");
			return 2;
		}
	
	/* Test if output file already exists. */
	if( access(outname,F_OK) ) {	/* File does not exist.  Test file opening. */

		testfile=fopen(outname,&outmode);

		if( testfile==NULL ) {	/* File could not be created. */
			printf("Specified output file cannot be created.\n");
			return 2;
		}

		fclose(testfile);	/* Close and remove file for now if creation is possible. */
		remove(outname);
		printf("Output file ok.\n");
	}
	else {	/* File exists.  Test if it is writable. */

		if( access(outname,W_OK) ) {	/* Not writable. */
			printf("You do not have write access to the specified output file.\n");
			return 2;
		return 0;
		}
		else {	/* File exists and is writable.  Ask for overwrite or append mode. */
			printf("Specified output file already exists.\nover(w)rite (a)ppend (q)uit? ");
			scanf("%c",&outmode);

			/* Test if user selected valid output mode. */
			while( outmode!='w' && outmode!='a' ) {
				if( outmode=='q' ) {	/* User selected quit. */
					return 1;
				}
				printf("Invalid output mode selected.\nover(w)rite (a)ppend (q)uit? ");
				scanf("%c",&outmode);
			}
		}
	}

	/*if( argc==4 ) return 0;	Exit successfully if no further arguments. */

	/* Process further arguments, with flags and values alternating. */
	/* for(i=4;i<argc;i+=2) {
		if( !strcmp(argv[i],"a") && argv[i+1] ) {
			sprintf(a,"%s",argv[i+1]);
		}
		else if( !strcmp(argv[i],"b") && argv[i+1] ) {
			sprintf(b,"%s",argv[i+1]);
		}
		else if( !strcmp(argv[i],"c") && argv[i+1] ) {
			sprintf(c,"%s",argv[i+1]);
		}
		else {  If argument was not a valid flag, exit
			printf("Optional arguments must be of the form 'name VALUE'.\n");
			printf("Valid names: "OPTIONS);
			return 1;
		}

	} */

	return 0;
}

void init_random ( )	/* Initialise random number seed */
{
	struct timeb tp;

	ftime(&tp);
	srandom( (tp.time*1000)+tp.millitm ); /* Use time in milliseconds to ensure no repeats */
}

int init_adatoms ( )	/* Initialise adatom and occupation arrays, placing atoms at random positions */
{
	int i, unallocated=n_sites-1, index;
	int free_sites[n_sites];

	occupations=(int *) calloc(n_sites,sizeof(int));	/* Allocate occupation and adatom coordinate arrays */
	adatoms=(int *) calloc(n_co,sizeof(int));

	if( occupations==NULL || adatoms==NULL ) return 3;	/* Test for successful allocation */

	/*for(i=0;i<n_co;i++) {
		do {
			adatoms[i]=(n_sites-1)*random()/RAND_MAX; Place adatoms at randomly chosen coordinates 
		} while ( occupations[ adatoms[i] ] );	 Test if site is already occupied
	
		occupations[ adatoms[i] ]=1;	 Set occupation state for chosen coordinate 
	}*/

	for(i=0;i<n_sites;i++) free_sites[i]=i;	/* Set free sites matrix: initially all sites are included */

	for(i=0;i<n_co;i++) {
		index=(unallocated)*random()/RAND_MAX;	/* Choose random free site */
		adatoms[i]=free_sites[index];	/* Place adatom at site */
		occupations[ adatoms[i] ]=1;	/* Set occupation state for chosen coordinate */
		free_sites[index]=free_sites[unallocated];/* Replace free site coordinate with highest unalloc. one */
		unallocated--;	/* Decrement the number of remaining free sites */
	}

	return 0;
}

void error_close ( int code )	/* Terminates the program on error with appropriate message and return value. */
{
	if( code==1 )
		printf("(q)uit selected.\nExiting %s...\n",name);
	else if( code==2 )
		printf("Invalid command-line arguments.\nExiting %s...\n",name);
	else if( code==3 )
		printf("Memory allocation failed when initialising arrays.\nExiting %s...\n",name);		
	else {
		printf("Program failed with unknown error.\nExiting %s...\n",name);
	}

	exit(code);
}
