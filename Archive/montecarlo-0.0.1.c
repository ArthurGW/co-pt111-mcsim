/* Monte-carlo simulation of a cO/Pt(111) surface system to replicate observed dynamics.	*/
/* Written by arthur Gordon-Wright, ajgw20@bath.ac.uk.						*/
/* 												*/
/* Version number 0.0.1, date 18/11/09								*/
/* Update notes: added command-line argument processing. All tests passed.			*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include <ctype.h>*/
/*#include <math.h>*/
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>

#define SQRT2 1.414213562373095 /* Square root of 2, to save repeated calculation */
#define SQRT3 1.732050807568877 /* Square root of 3, " */
#define UTT 210
#define USAGE "Usage: %s dim n_co outname [a A b B c C].\n"
/*#define OPTIONS "a, b, c.\n"*/

char name[30], outname[40], outmode='w', a[10]="A", b[10]="B", c[10]="C";
int dim, n_co;

/* Function declarations */
int args ( int , char ** );
void error_close( int );

int main ( int argc, char *argv[] )
{
	int arg_status;
	srandom( (unsigned int)time(NULL) );

	if( argc==1 || strcmp(argv[1],"help")==0 ) {
		printf(USAGE, argv[0]);
		return 3;
	}

	if( (arg_status=args( argc, argv )) ) { 
		error_close(arg_status);
	}

	printf("name=%s, dim=%d, n_co=%d, outname=%s\n", name, dim, n_co, outname );
	printf("a=%s, b=%s, c=%s\n", a, b, c );

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
	
	/* Extract lattice dimension and number of CO adatoms. */
	dim=atoi(argv[1]);
	n_co=atoi(argv[2]);

	/* Test for successful extraction. */
	if( !dim || !n_co ) {
		printf("Invalid dimension or number of adatoms.  Please specify integer values.\n");
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

void error_close ( int status )
{
	if( status==1 )
		printf("(q)uit selected.\nExiting %s...\n",name);
	else if( status==2 )
		printf("Invalid command-line arguments.\nExiting %s...\n",name);	
	else {
		printf("Program failed with unknown error.\nExiting %s...\n",name);
	}

	exit(status);
}
