/* Monte-Carlo simulation of a CO/Pt(111) surface system to replicate observed dynamics.	*/
/* Written by Arthur Gordon-Wright, ajgw20@bath.ac.uk.						*/
/* 												*/
/* Version number 0.2.0, date 29/11/09								*/
/* Update notes: Added repulsion								*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include <ctype.h>*/
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <unistd.h>
#include <sys/stat.h>

#define ONEOVSQRT2 0.707106781186548	/* One over the square root of 2, to save repeated calculation */
/*#define SQRT3 1.732050807568877		/* Square root of 3 */
/*#define SQRT3OV8 0.612372435695795	/* Square root of 3 over 8 */
#define SQRT1OV8 0.353553390593274	/* Square root of 1 over 8 */
#define UTT 210		/* Top-top diffusion barrier */
#define USAGE "Usage: %s dim n_co int_mode out_name\n"
#define OPTIONS "g nsweeps NSWEEPS\n"
#define LINEWIDTH 22	/* Number of '-' separators in output file dividers */
#define OUT_NAME_LEN 100	/* Maximum length of output file name (and specified path) */
#define MAX_OUTSIZE 1048576	/* Maximum output file size (1Mb) */
#define GRES 100

/* Variable declarations */
char name[30], out_name[OUT_NAME_LEN], outmode='w';
int dim, n_sites, n_co, int_mode, skip_sweeps=10000, n_sweeps=10000, g=0, *occupations, *adatoms, *G, bins;
time_t start_time;
double bin_res;

/* Function declarations */
int args ( int , char ** );
void init_random ( );
int init_adatoms ( );
void error_close ( int );
int output ( int );
int test_file ( );
double interaction ( int, int );
char * int_name ( );
int sweep ( );
int select_move ( int , int * );
double separation ( int, int );
void pair ( );

int main ( int argc, char *argv[] )
{
	int err_code, i, j, k, lineindex, bin;

	if( argc==1 || strcmp(argv[1],"help")==0 ) { /* Test if help has been selected or no arguments set */
		printf(USAGE, argv[0]);
		return 3;
	}

	time( &start_time );

	init_random();	/* Set random seed */

	if( (err_code=args( argc, argv )) ) { /* Get and test command-line arguments, exit on failure */
		error_close(err_code);
	}

	output(1);

	bins=GRES*dim;
	bin_res=ONEOVSQRT2/GRES;
	G=(int *)calloc(bins,sizeof(int));

	for(i=0;i<bins;i+=2) G[i]=0;

	if( ( err_code=init_adatoms() ) ) {	/* Initialise adatom positions, exit on failure */
		error_close(err_code);
	}

	output(2);

	printf("name=%s, dim=%d, n_co=%d, int_mode=%d, n_sweeps=%d, out_name=%s\n", name, dim, n_co, int_mode, n_sweeps, out_name );
	/*printf("a=%s, b=%s, c=%s\n", a, b, c );*/

	printf("\nOccupations.\n");
	printf("+ ");
	for(j=0;j<dim;j++) printf("- ");
	printf("+\n");
	for(j=0;j<n_sites;j++) {
		if( !((j)&(dim-1)) ) {
			lineindex=j/dim;
			for(k=0;k<lineindex+1;k++) printf(" ");
			printf("\\ ");
		}
		if( occupations[j] ) printf("%d ",occupations[j]);
		else printf("  ");
		if( !((j+1)&(dim-1)) ) printf("\\\n");
	}
	for(j=0;j<dim+1;j++) printf(" ");
	printf("+ ");
	for(j=0;j<dim;j++) printf("- ");
	printf("+\n");

	printf("\nCoordinates.\n");
	for(i=0;i<n_co;i++)
		printf("%d ",adatoms[i]);
	printf("\n\n");

	for(i=0;i<skip_sweeps;i++) sweep();

	for(i=0;i<n_sweeps;i++) {
		sweep();

		if( g ) {
			system("clear");
			printf("Sweep #%d:\n",i+1);
			printf("+ ");
			for(j=0;j<dim;j++) printf("- ");
			printf("+\n");
			for(j=0;j<n_sites;j++) {
				if( !((j)&(dim-1)) ) {
					lineindex=j/dim;
					for(k=0;k<lineindex+1;k++) printf(" ");
					printf("\\ ");
				}
				if( occupations[j] ) printf("%d ",occupations[j]);
				else printf("  ");
				if( !((j+1)&(dim-1)) ) printf("\\\n");
			}
			for(j=0;j<dim+1;j++) printf(" ");
			printf("+ ");
			for(j=0;j<dim;j++) printf("- ");
			printf("+\n");
			sleep(1);
		}

		pair();
	}

	/*printf("\aNo. of bins: %d Bin resolution: %1.6f\n",bins,bin_res);
	printf("Final separations.\n");
	for(i=0;i<n_co-1;i++) {
		printf("Coord1: %d Coord2: %d\n", adatoms[i],adatoms[i+1]); 
		printf("Bin: %d\n",separation(adatoms[i],adatoms[i+1]),
				(int)floor(.5+separation(adatoms[i],adatoms[i+1])*bin_res) );
	}
	printf("\n\n");*/

	output(8);

	if( ( err_code=output(9) ) ) {	 /* Write simulation results to file, exit on error */
		error_close(err_code);
	}

	return 0;
}

int args ( int argc, char *argv[] )	/* Processes command-line arguments and sets global variables accordingly */
{					/* Returns 0 if successful */
	int i, error_code;

	/* Extract program name */
	sprintf(name,"%s",argv[0]);

	/* Test for minimum valid number of arguments */
	if( argc<4 ) {
		printf("Incorrect number of arguments.\n"USAGE, argv[0]);
		return 2;
	}
	
	/* Extract lattice dimension, number of sites and number of CO adatoms */
	dim=atoi(argv[1]);
	n_sites=dim*dim;
	n_co=atoi(argv[2]);

	/* Test for successful extraction */
	if( !dim || !n_co ) {
		printf("Invalid dimension or number of adatoms.  Please specify integer values.\n");
		return 2;
	}

	/* Check that specified lattice dimension is a power of 2, to allow bitwise modulus function */
	if ( (dim-1) & dim ) {
		printf("For program efficiency, please specify a lattice dimension that is an integer\npower of two.\n");
		return 2;
	}

	/* Check for valid specified coverage (i.e. number of sites is >= number of adatoms */
	if( (dim*dim)<n_co ) {
		printf("Specified coverage is greater than 100%!\n");
		return 2;
	}

	/* Extract required interaction, test for valid option */
	int_mode=atoi(argv[3]);
	if( int_mode<0 || int_mode>3 ) {
		printf("Please enter an interaction mode between 0 and 3.\n");
		return 2;	
	}

	/* Extract output file name */
	if( strlen(argv[4])>100 ) {
		printf("Specified output file name must be less than 100 characters.\n");
		return 2;
	}
	sprintf(out_name,"%s",argv[4]);

	if( error_code=test_file(1) ) return error_code;

	if( argc==5 ) return 0;	/* Exit successfully if no further arguments */

	/* Process further arguments */
	for(i=5;i<argc;i++) {
		if( !strcmp(argv[i],"g") ) {
			g=1;
		}
		else if( !strcmp(argv[i],"nsweeps") && argv[i+1] ) {
			n_sweeps=atoi(argv[i+1]);
			if( !n_sweeps ) {
				printf("Please enter an integer value for optional parameter nsweeps.\n");
				return 2;
			}
			i++;
		}
		else if( !strcmp(argv[i],"ssweeps") && argv[i+1] ) {
			skip_sweeps=atoi(argv[i+1]);
			if( !skip_sweeps ) {
				printf("Please enter an integer value for optional parameter ssweeps.\n");
				return 2;
			}
			i++;
		}
		else {  /* If argument was not a valid flag, exit */
			printf("Invalid optional argument.\n");
			printf("Valid arguments: "OPTIONS);
			return 2;
		}

	}

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
	int i, unallocated=n_sites, index;
	int free_sites[n_sites];

	printf("Allocating occupation matrix...");
	occupations=(int *) calloc(n_sites,sizeof(int));	/* Allocate occupation and adatom coordinate arrays */
	printf("done!\nAllocating adatom coordinate matrix...");
	adatoms=(int *) calloc(n_co,sizeof(int));
	printf("done!\n");

	if( occupations==NULL || adatoms==NULL ) return 3;	/* Test for successful allocation */

	for(i=0;i<n_sites;i++) free_sites[i]=i;	/* Set free sites matrix: initially all sites are included */

	for(i=0;i<n_co;i++) {
		index=(int)floor(unallocated*random()/RAND_MAX);	/* Choose random free site */
		adatoms[i]=free_sites[index];	/* Place adatom at site */
		occupations[ adatoms[i] ]=1;	/* Set occupation state for chosen coordinate */
		free_sites[index]=free_sites[unallocated];/* Replace free site coordinate with highest unalloc. one */
		unallocated--;	/* Decrement the number of remaining free sites */
	}

	return 0;
}

int output ( int out_num )	/* Writes data to output file.  Attempts writing three times, then asks user to */	
{				/* try a different file, then another, then fails 				*/
	int i, n, attempts=1, success=0, file_no=1, normalisation=(n_co)*(n_co-1)*n_sweeps;
	FILE *tmp, *out;
	char tmp_name[15];
	char * buffer;

	sprintf(tmp_name,"%d.tmp",(int)start_time);

	if( out_num==1 ) {
		printf("tmp_name=%s\n",tmp_name);
	}

	tmp=fopen(tmp_name,"a");	/* Open temp file */

	/* Write output to temp file. out_num flag sets which stage output to write */	
	if( out_num==1 ) {	/* Header section */

		fprintf(tmp,"CO/Pt111 System Monte-Carlo Simulation.\t\tStarted: %s\n\n",ctime( &start_time ) );

		fprintf(tmp,"name=%s, dim=%d, n_co=%d, n_sweeps=%d, out_name=%s\n", name, dim, n_co, n_sweeps, out_name );
		fprintf(tmp,"Interaction mode: %s\n", int_name() );

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}
	else if( out_num==2 ) {	/* Initial positions section */
		/*fprintf(tmp,"Initial coordinates:\n");

		for(i=0;i<n_sites;i++) {
			fprintf(tmp,"%d ",occupations[i]);
			if( !((i+1)&(dim-1)) )
				fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");*/
	}
	else if( out_num==8 ) {	/* Closing section */
		/*fprintf(tmp,"Final coordinates:\n");

		for(i=0;i<n_sites;i++) {
			fprintf(tmp,"%d ",occupations[i]);
			if( !((i+1)&(dim-1)) )
				fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");*/

		fprintf(tmp,"Pair correlation function.\n");
		fprintf(tmp,"Bin resolution: %1.9f, normalisation factor: %d\n",bin_res,normalisation);
		fprintf(tmp,"Bin centre\tProbability\tNumber in bin\n",ONEOVSQRT2/GRES);
		for(i=0;i<bins;i++)
			if( G[i] )
				fprintf(tmp,"%2.7f\t%1.9f\t%d\n",i/GRES,(float)G[i]/normalisation,G[i]);

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}


	fclose(tmp);

	if( out_num==9 ) {	/* Transfer temporary output to specified file */ 
		while( !success ) {
			success=1; /* Assume run is successful unless later set otherwise */

			if( !(tmp=fopen(tmp_name,"r")) ) success=0;	/* Open temp file for reading */

			if( !(out=fopen(out_name,&outmode)) ) success=0;	/* Open output file */

			buffer=(char *) malloc(MAX_OUTSIZE);	/* Read temp to buffer then write to output */
			n=fread(buffer,1,MAX_OUTSIZE,tmp);
			printf("n=%d\n",n);
			if( ferror(tmp) ) success=0;
			else {
				fwrite(buffer,1,n,out);
				fputc('\n',out);
				if( ferror(out) ) success=0;
			}

			fclose(tmp);	/* Close files */
			fclose(out);

			if( !success ) {
				printf("Attempt #%d to write output file failed.\n",attempts);

				attempts++;

				if( attempts==4 ) {	/* After three attempts at using a file, try another, 	*/
					file_no++;	/* assuming less than three have been tried 		*/
					if( file_no<4 ) {
						attempts=1;
						do {	/* Keep asking for a new file until a correct one is specified */
							printf("Cannot write to output file.  Please specify another:\n");
							scanf("%s",&out_name);
						} while( test_file() );
					}
					else {	/* If next file would be the fourth, exit */
						return 4;
					}
				}
			}
		}
		remove(tmp_name);
	}

	printf("Output #%d successfully written.\n",out_num);
	return 0;
}

void error_close ( int code )	/* Terminates the program on error with appropriate message and return value */
{
	if( code==1 )
		printf("(q)uit selected.\nExiting %s...\n",name);
	else if( code==2 )
		printf("Invalid command-line arguments.\nExiting %s...\n",name);
	else if( code==3 )
		printf("Memory allocation failed when initialising arrays.\nExiting %s...\n",name);
	else if( code==4 )
		printf("Output file writing failed.  Temporary file not deleted.\nExiting %s...\n",name);		
	else {
		printf("Program failed with unknown error.\nExiting %s...\n",name);
	}

	exit(code);
}

int test_file ( )	/* Tests output file is valid. Returns 0 on success or various errors on failure */
{
	FILE *testfile;
	struct stat buf;

	/* Check if output file is a directory */
	if( !stat(out_name,&buf) )
		if( S_ISDIR(buf.st_mode) ) {
			printf("Specified output file is a directory.\n");
			return 2;
		}
	
	/* Test if output file already exists */
	if( access(out_name,F_OK)==-1 ) {	/* File does not exist.  Test file opening */

		testfile=fopen(out_name,&outmode);

		if( testfile==NULL ) {	/* File could not be created */
			printf("Specified output file cannot be created.\n");
			return 2;
		}

		fclose(testfile);	/* Close and remove file for now if creation is possible */
		remove(out_name);
	}
	else {	/* File exists.  Test if it is writable */

		if( access(out_name,W_OK) ) {	/* Not writable */
			printf("You do not have write access to the specified output file.\n");
			return 2;

		}
		else {	/* File exists and is writable.  Ask for overwrite or append mode */
			printf("Specified output file already exists.\nover(w)rite (a)ppend (q)uit? ");
			scanf("%c",&outmode);

			/* Test if user selected valid output mode */
			while( outmode!='w' && outmode!='a' ) {
				if( outmode=='q' ) return 1;	/* User selected quit */
				printf("Invalid output mode selected.\nover(w)rite (a)ppend (q)uit? ");
				scanf("%c",&outmode);
			}
			if( outmode=='w' ) printf("Overwrite mode selected.  ");
			else printf("Append mode selected.  ");
		}
	}

	printf("Output file ok.\n");	
	return 0;
}

double interaction ( int coord1, int coord2 )
{
	return 0;	
}

double separation ( int coord1, int coord2 )
{
	unsigned int i, num_mirrors, sep, new_sep; 
	signed int twodx, twodx_m, dy, dy_m;

	static signed int mirrors[8];
	signed int current_mirrors[6]={0,0,0,0,0,0};
	static int flag=0;

	if( !flag ) {
		mirrors[0]=dim;
		mirrors[1]=2*dim;
		mirrors[2]=3*dim-1;
		mirrors[3]=2*dim-2;
		mirrors[4]=2*dim-1;
		mirrors[5]=-2;
		mirrors[6]=dim-1;
		mirrors[7]=-2*dim-2;
		flag++;
	}

	dy=(int)floor(coord2/dim) - (int)floor(coord1/dim);
	twodx=2 * ( (coord2&(dim-1)) - (coord1&(dim-1)) ) + dy;

	sep=(twodx*twodx)+(3*dy*dy);

	if( twodx-dy>0 ) {
		if( dy>0 ) {
			num_mirrors=3;
			current_mirrors[0]=-1*mirrors[0];
			current_mirrors[1]=-1*mirrors[1];
			current_mirrors[2]=-1*mirrors[2];
			current_mirrors[3]=-1*mirrors[3];
			current_mirrors[4]=-1*mirrors[4];
			current_mirrors[5]=-1*mirrors[5];
		}
		else if( dy<-1 ) {
			num_mirrors=3;
			current_mirrors[0]=mirrors[0];
			current_mirrors[1]=mirrors[1];
			current_mirrors[2]=-1*mirrors[4];
			current_mirrors[3]=-1*mirrors[5];
			current_mirrors[4]=-1*mirrors[6];
			current_mirrors[5]=-1*mirrors[7];
		}
		else {
			num_mirrors=1;
			current_mirrors[0]=-1*mirrors[4];
			current_mirrors[1]=-1*mirrors[5];
		}
	}
	else if( twodx-dy<0 ) {
		if( dy>1 ) {
			num_mirrors=3;
			current_mirrors[0]=-1*mirrors[0];
			current_mirrors[1]=-1*mirrors[1];
			current_mirrors[2]=mirrors[4];
			current_mirrors[3]=mirrors[5];
			current_mirrors[4]=mirrors[6];
			current_mirrors[5]=mirrors[7];
		}
		else if( dy<0 ) {
			num_mirrors=3;
			current_mirrors[0]=mirrors[0];
			current_mirrors[1]=mirrors[1];
			current_mirrors[2]=mirrors[2];
			current_mirrors[3]=mirrors[3];
			current_mirrors[4]=mirrors[4];
			current_mirrors[5]=mirrors[5];
		}
		else {
			num_mirrors=1;
			current_mirrors[0]=mirrors[4];
			current_mirrors[1]=mirrors[5];
		}
	}
	else {
		num_mirrors=1;
		if( dy>0 ) {
			current_mirrors[0]=-1*mirrors[0];
			current_mirrors[1]=-1*mirrors[1];
		}
		else {
			current_mirrors[0]=mirrors[0];
			current_mirrors[1]=mirrors[1];
		}
	}


	for(i=0;i<2*num_mirrors;i+=2) {
		twodx_m=twodx+current_mirrors[i];
		dy_m=dy+current_mirrors[i+1]/2;

		new_sep=(twodx_m*twodx_m)+(3*dy_m*dy_m);
		sep = ( sep <= new_sep ) ? sep : new_sep;
	}

	return SQRT1OV8*sqrt(sep);

}

int select_move ( int coord, int * nn )
{
	int index;

	/* Select random move */
	index=(int)floor(6*random()/RAND_MAX);

	if( occupations[ nn[index] ] ) return coord;
	else return nn[index];
}

int sweep ( )
{
	int i,j,k,lineindex,atom,coord,move_coord,nn[6],nns,mnns;

	for(i=0;i<n_co;i++) {

		/* Select random adatom */
		atom=(int)floor(n_co*random()/RAND_MAX);
		coord=adatoms[atom];

		nns=get_nns(coord,nn);	/* Populate nearest neighbour coordinate matrix */

		move_coord=select_move(coord,nn);	/* Select random move */
		if( move_coord==coord ) continue;	/* If no move selected, go to next atom */

		mnns=get_nns(move_coord,nn)-1;

		/* Simple test interaction model */
		if( !int_mode || mnns<nns || (random()/RAND_MAX)< (5+nns-mnns)/6 ) {
			occupations[adatoms[atom]]=0;
			occupations[move_coord]=1;
			adatoms[atom]=move_coord;
		}
	}

	return 0;
}

char * int_name ( )
{
	if( int_mode==0 ) return "None";
	if( int_mode==1 ) return "McEwen";
	if( int_mode==2 ) return "Petrova";
	if( int_mode==3 ) return "Persson";
}

void pair ( )
{
	int i,j,bin;
	double ijsep;

	for(i=0;i<n_co;i++) {
		j=0;
		while( j<n_co ) {
			if( i==j ) {
				j++;
				continue;
			}
				
			ijsep=separation(adatoms[i],adatoms[j]);
			bin=(int)floor(.5+ijsep*GRES/ONEOVSQRT2);

			G[bin]++;
			
			j++;
		}
	}
}

int get_nns ( int coord, int * nn )
{
	int i,sum=0;

	nn[0]=(coord+1)&(n_sites-1);
	nn[1]=(coord-1+n_sites)&(n_sites-1);
	nn[2]=(coord+dim)&(n_sites-1);
	nn[3]=(coord-dim+n_sites)&(n_sites-1);
	nn[4]=(coord+(dim-1)+n_sites)&(n_sites-1);
	nn[5]=(coord-(dim-1)+n_sites)&(n_sites-1);

	for(i=0;i<6;i++) sum+=occupations[ nn[i] ];
}
