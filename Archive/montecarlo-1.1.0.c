/* Monte-Carlo simulation of a CO/Pt(111) surface system to replicate observed dynamics.	*/
/* Written by Arthur Gordon-Wright, ajgw20@bath.ac.uk.						*/
/* 												*/
/* Version number 0.4.0, date 08/12/09								*/
/* Update notes: Reduced area of k-space used							*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fftw3.h>

#define ONEOVSQRT2 0.707106781186548	/* One over the square root of 2, to save repeated calculation */
#define SQRT1OV8 0.353553390593274	/* Square root of 1 over 8 */
#define UTT 210		/* Top-top diffusion barrier */
#define USAGE "Usage: %s dim n_co int_mode out_name\n"
#define OPTIONS "g nsweeps NSWEEPS tfact TFACT tsteps TSTEPS repeats REPEATS\n"
#define LINEWIDTH 22		/* Number of '-' separators in output file dividers */
#define OUT_NAME_LEN 100	/* Maximum length of output file name (and specified path) */
#define MAX_OUTSIZE 231072	/* Maximum output file size (~0.1Mb) */
/*#define GRES 100	/* Resolution for pair correlation function - number of bins per lattice vector */

/* Variable declarations */
char name[30], out_name[OUT_NAME_LEN], outmode='w';
int dim, n_sites, n_co, int_mode, skip_sweeps=100000, n_sweeps=100000, g=0, bins, time_steps=20, time_factor=2, repeats=10;
time_t start_time;
double bin_res;
fftw_complex *fftout, *fft_steps, *ISF_av, *temp_ISF;
double *temp_occs;
int *occupations, *adatoms;

/* Function declarations */
int args ( int , char ** );
void init_random ( );
int init_adatoms ( );
void error_close ( int );
int output ( int );
int test_file ( );
double interaction ( int, int );
char * int_name ( );
void sweep ( );
int select_move ( int , int * );
void pair ( );
void fft_run ( int, signed int );
void fft_mult ( int, int );
void display_data ( int );

int main ( int argc, char *argv[] )
{
	int err_code, i, j, k, lineindex, bin, time_index, repeat=0, sweeps, next_step;

	if( argc==1 || strcmp(argv[1],"help")==0 ) { /* Test if help has been selected or no arguments set */
		printf(USAGE, argv[0]);
		return 3;
	}

	time( &start_time );	/* Get current time */

	init_random();	/* Set random seed */

	if( (err_code=args( argc, argv )) ) { /* Get and test command-line arguments, exit on failure */
		error_close(err_code);
	}

	if( ( err_code=init_adatoms() ) ) {	/* Initialise adatom positions, exit on failure */
		error_close(err_code);
	}

	output(1); output(2);	/* Output headers and initial positions to temp. file */

	printf("name=%s, dim=%d, n_co=%d, int_mode=%d, n_sweeps=%d, skip_sweeps=%d, out_name=%s\n", name, dim, n_co, int_mode, n_sweeps, skip_sweeps, out_name );	/* Display some parameters */

	printf("\nOccupations.\n");	/* Display current occupations matrix */
	display_data(0);

	for(i=0;i<skip_sweeps;i++) sweep();	/* Skip some sweeps */

	/* FFT - initialise input and output arrays, then run (with planning) */
	temp_occs = (double*) fftw_malloc( sizeof(double) * n_sites );	/* Copy of occupations matrix */
	fftout=(fftw_complex *) fftw_malloc( 2 * sizeof(double) * n_sites );	/* FFT output array */
	fft_steps=(fftw_complex *) malloc( (time_steps+1) * 2 * sizeof(double) * n_sites );
								/* Array for storing FFTs for each time step */
	fft_run(1,time_steps);	/* Plan and run FFT */

	printf("\nTransformed data.\n");	/* Display FFT result */
	display_data(1);

	output(4);	/* Output FFT result */

	ISF_av=(fftw_complex *) malloc( 2 * sizeof(double) * (dim+2) * (time_steps+1) );
								/* Array for storing ISF means and variances */
	temp_ISF=(fftw_complex *) malloc( 2 * sizeof(double) * (dim+2) );	/* Array for storing each ISF result */

	for(i=0;i<(dim+2)*(time_steps+1);i++) {	/* Set ISF_av to zero intially */
		ISF_av[i][0]=0.;
		ISF_av[i][1]=0.;
	}

	do {	/* Main loop (over repeats) */
		time_index=1;	/* Reset time index, sweeps index */
		sweeps=0;
		next_step=pow(time_factor,1);	/* Calculate next sweep index to record at */

		for(i=0;i<n_sites/2+dim;i++) {	/* Set t=0 FFT record to current FFT result */
			fft_steps[i][0]=fft_steps[time_steps*(n_sites/2+dim)+i][0];
			fft_steps[i][1]=fft_steps[time_steps*(n_sites/2+dim)+i][1];
		}

		do {	/* Loop over time steps */
			sweep();	/* Run sweep */
			sweeps++;

			if( g ) {	/* If in graphical mode, display occupations matrix */
				system("clear");
				printf("Sweep #%d:\n",sweeps);
				display_data(0);
				sleep(1);
			}

			if( sweeps==next_step ) {	/* If sweep is to be recorded */
				fft_run(0,time_index);	/* Calculate FFT, store in fft_steps */
				time_index++;		/* Calculate next sweep to be recorded */
				next_step=pow(time_factor,time_index);
			}

		} while(time_index<time_steps+1);
		
		for(time_index=0;time_index<time_steps+1;time_index++) {
			fft_mult(0,time_index);	/* Multiply FFTs for t=time_index and t=0 */
			for(i=0;i<dim/2+1;i++) {	/* Transfer results from temp. array to ISF average array */
				ISF_av[time_index*(dim+2)+2*i][0]+=temp_ISF[i][0]/repeats;
				ISF_av[time_index*(dim+2)+2*i][1]+=temp_ISF[i][1]/repeats;
				ISF_av[time_index*(dim+2)+2*i+1][0]+=(temp_ISF[i][0]*temp_ISF[i][0])/repeats;
				ISF_av[time_index*(dim+2)+2*i+1][1]+=(temp_ISF[i][1]*temp_ISF[i][1])/repeats;
			}		
		}

		repeat++;
		if( !(repeat&8191) ) printf("Repeat no.: %d\n",repeat);	/* Display repeat number occasionally */

	} while(repeat<repeats);

	free(temp_ISF);	free(fft_steps); free(adatoms);	fftw_free(temp_occs); /* Free arrays that are no longer needed */

	for(time_index=0;time_index<time_steps+1;time_index++) {
		for(i=0;i<dim+2;i+=2) {	/* Calculate variance by subtracting mean^2 from mean-squared value */
			ISF_av[time_index*(dim+2)+i+1][0]-=
					ISF_av[time_index*(dim+2)+i][0]*ISF_av[time_index*(dim+2)+i][0];
			ISF_av[time_index*(dim+2)+i+1][1]-=
					ISF_av[time_index*(dim+2)+i][1]*ISF_av[time_index*(dim+2)+i][1];
		}		
	}

	output(6);	/* Output ISF data then free array */
	free(ISF_av);

	printf("\nFinal occupations.\n");	/* Display final positions */
	display_data(0);

	printf("\nTransformed data.\n");	/* Display final FFT */
	display_data(1);

	output(2); output(4);	/* Output final positions and FFT, then free arrays */
	free(occupations); fftw_free(fftout);

	if( ( err_code=output(11) ) ) {	 /* Write simulation results to file, exit on error */
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
		else if( !strcmp(argv[i],"tfact") && argv[i+1] ) {
			time_factor=atoi(argv[i+1]);
			if( !time_factor ) {
				printf("Please enter an integer value for optional parameter tfact.\n");
				return 2;
			}
			i++;
		}
		else if( !strcmp(argv[i],"tsteps") && argv[i+1] ) {
			time_steps=atoi(argv[i+1]);
			if( !time_steps ) {
				printf("Please enter an integer value for optional parameter tsteps.\n");
				return 2;
			}
			i++;
		}
		else if( !strcmp(argv[i],"repeats") && argv[i+1] ) {
			repeats=atoi(argv[i+1]);
			if( !repeats ) {
				printf("Please enter an integer value for optional parameter repeats.\n");
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

	for(i=0;i<n_sites;i++) {
		occupations[i]=0;	/* Set occupations matrix to all unoccupied */
		free_sites[i]=i;	/* Set free sites matrix: initially all sites are included */
	}

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
	int i, j, n, attempts=1, success=0, file_no=1, normalisation;
	FILE *tmp, *out;
	static char tmp_name[15];
	char * buffer;

	if( out_num==1 ) {	/* Create and display temp. file name on first run */
		sprintf(tmp_name,"%d.tmp",(int)start_time);
		printf("tmp_name=%s\n",tmp_name);
	}

	tmp=fopen(tmp_name,"a");	/* Open temp file */

	/* Write output to temp file. out_num flag sets which stage output to write */	
	if( out_num==1 ) {	/* Header section */

		fprintf(tmp,"CO/Pt111 System Monte-Carlo Simulation.\t\tStarted: %s\n\n",ctime( &start_time ) );

		fprintf(tmp,"name=%s, dim=%d, n_co=%d, n_sweeps=%d, skip_sweeps=%d, out_name=%s\n", name, dim, n_co, n_sweeps, skip_sweeps, out_name );
		fprintf(tmp,"Interaction mode: %s\n", int_name() );

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}
	else if( out_num==2 ) {	/* Initial positions section */
		fprintf(tmp,"Coordinates:\n");

		for(i=0;i<n_sites;i++) {
			fprintf(tmp,"%d ",occupations[i]);
			if( !((i+1)&(dim-1)) )
				fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}
	else if( out_num==3 ) {	/* Temp. occs. section */
		fprintf(tmp,"Temp. occs.:\n");

		for(i=0;i<n_sites;i++) {
			fprintf(tmp,"%f\t",temp_occs[i]);
			if( !((i+1)&(dim-1)) )
				fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}
	else if( out_num==4 ) {	/* Most recent Fourier transform */
		fprintf(tmp,"Fourier transform:\n");
		for(i=0;i<n_sites/2+dim;i++) {
			fprintf(tmp,"%.3f\t%.3f\t\t",fftout[i][0],fftout[i][1]);
			if( !((i+1)%(dim/2+1)) ) fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}	
	else if( out_num==5 ) {	/* Fourier transform over all time steps */
		fprintf(tmp,"Fourier timestep matrix:\n");

		for(i=0;i<(n_sites/2+dim)*(time_steps+1);i++) {
			fprintf(tmp,"%.3f\t%.3f\t\t",fft_steps[i][0],fft_steps[i][1]);

			if( !((i+1)%(dim/2+1)) ) fprintf(tmp,"\n");
			if( !((i+1)%(n_sites/2+dim)) ) fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}
	else if( out_num==6 ) {	/* ISF real and imaginary means and variances */
		fprintf(tmp,"ISF:\n");
		fprintf(tmp,"time_factor: %d, time_steps: %d, repeats: %d\n\n",time_factor,time_steps,repeats);
		fprintf(tmp,"Real\nTimestep 1\tTimestep 2\t");
		for(i=0;i<dim/2+1;i++) fprintf(tmp,"k%d\t",i);
		fprintf(tmp,"\n");
		for(j=0;j<(time_steps+1);j++) {
			fprintf(tmp,"%d\t%d\t",0,j);
			for(i=0;i<dim+2;i+=2)
				fprintf(tmp,"%.5f\t",ISF_av[j*(dim+2)+i][0]);
			fprintf(tmp,"\n");
		}

		fprintf(tmp,"\nImaginary\nTimestep 1\tTimestep 2\t");
		for(i=0;i<dim/2+1;i++) fprintf(tmp,"k%d\t",i);
		fprintf(tmp,"\n");
		for(j=0;j<(time_steps+1);j++) {
			fprintf(tmp,"%d\t%d\t",0,j);
			for(i=0;i<dim+2;i+=2)
				fprintf(tmp,"%.5f\t",ISF_av[j*(dim+2)+i][1]);
			fprintf(tmp,"\n");
		}

		fprintf(tmp,"\nReal variance\nTimestep 1\tTimestep 2\t");
		for(i=0;i<dim/2+1;i++) fprintf(tmp,"k%d\t",i);
		fprintf(tmp,"\n");
		for(j=0;j<(time_steps+1);j++) {
			fprintf(tmp,"%d\t%d\t",0,j);
			for(i=0;i<dim+2;i+=2)
				fprintf(tmp,"%.5f\t",ISF_av[j*(dim+2)+i+1][0]);
			fprintf(tmp,"\n");
		}
		fprintf(tmp,"\nImaginary variance\nTimestep 1\tTimestep 2\t");
		for(i=0;i<dim/2+1;i++) fprintf(tmp,"k%d\t",i);
		fprintf(tmp,"\n");
		for(j=0;j<(time_steps+1);j++) {
			fprintf(tmp,"%d\t%d\t",0,j);
			for(i=0;i<dim+2;i+=2)
				fprintf(tmp,"%.5f\t",ISF_av[j*(dim+2)+i+1][1]);
			fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}
	else if( out_num==10 ) {	/* Closing section */
		/*fprintf(tmp,"Final coordinates:\n");

		for(i=0;i<n_sites;i++) {
			fprintf(tmp,"%d ",occupations[i]);
			if( !((i+1)&(dim-1)) )
				fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");*/

		/*fprintf(tmp,"Pair correlation function.\n");
		normalisation=(n_co)*(n_co-1)*n_sweeps;
		fprintf(tmp,"Bin resolution: %1.9f, normalisation factor: %d\n\n",bin_res,normalisation);
		fprintf(tmp,"Bin centre\tProbability\tNumber in bin\n");
		for(i=0;i<bins;i++)
			if( G[i] )
				fprintf(tmp,"%2.4f\t%1.6f\t%d\n",(float)i/GRES,(float)G[i]/normalisation,G[i]);

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");*/
	}


	fclose(tmp);

	if( out_num==11 ) {	/* Transfer temporary output to specified file */ 
		while( !success ) {
			success=1; /* Assume run is successful unless later set otherwise */

			if( !(tmp=fopen(tmp_name,"r")) ) success=0;	/* Open temp file for reading */

			if( !(out=fopen(out_name,&outmode)) ) success=0;	/* Open output file */

			buffer=(char *) malloc(MAX_OUTSIZE);/* Read temp to buffer then write to output */
			n=fread(buffer,1,MAX_OUTSIZE,tmp);
			printf("n=%d\n",n);
			if( ferror(tmp) ) success=0;
			else {
				fwrite(buffer,1,n,out);
				fputc('\n',out);
				if( ferror(out) ) success=0;
			}

			free(buffer);	/* Free buffer */
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

		remove(tmp_name);	/* Remove temp. file */
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

void sweep ( )
{
	int i,atom,coord,move_coord,nn[6],nns,mnns;

	for(i=0;i<n_co;i++) {

		/* Select random adatom */
		atom=(int)floor(n_co*random()/RAND_MAX);
		coord=adatoms[atom];

		nns=get_nns(coord,nn);	/* Populate nearest neighbour coordinate matrix */

		move_coord=select_move(coord,nn);	/* Select random move */
		if( move_coord==coord ) continue;	/* If no move selected, go to next atom */

		mnns=get_nns(move_coord,nn)-1;	/* Calculate potential new number of nns, -1 as adatom	*/
						/*  will have moved from an occcupied neighbour 	*/

		/* Simple test interaction model */
		if( !int_mode || mnns<=nns ||(float)(random()/RAND_MAX) < (float)((5.+nns-mnns)/6.) ) {
			occupations[adatoms[atom]]=0;	/* Move adatom */
			occupations[move_coord]=1;
			adatoms[atom]=move_coord;
		}
	}
}

int select_move ( int coord, int * nn )
{
	int index;

	/* Select random move */
	index=(int)floor(6*random()/RAND_MAX);

	/* Test if selected site is occupied, if it is return unmoved coordinate */
	if( occupations[ nn[index] ] ) return coord;
	else return nn[index];	/* Otherwise return new coordinate */
}

char * int_name ( )
{
	if( int_mode==0 ) return "None";
	if( int_mode==1 ) return "Nearest neighbour";
	if( int_mode==2 ) return "Petrova";
	if( int_mode==3 ) return "Persson";
}

int get_nns ( int coord, int * nn )
{
	int j,sum=0;

	/* For given coord, calculate nearest neighbour coords using helical BCs */
	nn[0]=(coord+1)&(n_sites-1);
	nn[1]=(coord-1+n_sites)&(n_sites-1);
	nn[2]=(coord+dim)&(n_sites-1);
	nn[3]=(coord-dim+n_sites)&(n_sites-1);
	nn[4]=(coord+(dim-1)+n_sites)&(n_sites-1);
	nn[5]=(coord-(dim-1)+n_sites)&(n_sites-1);

	/* Calculate number of occupied neighbours */
	for(j=0;j<6;j++) sum+=occupations[ nn[j] ];

	return sum;
}

void fft_run ( int mode, int time_index )
{
	int i;
	static fftw_plan fft_occs, ifft_occs;

	if( mode ) {	/* Plan mode */
		for(i=0;i<n_sites;i++) temp_occs[i]=(double)occupations[i]; /* Set up temp. occupations matrix, as */
									    /* matrix is modified during planning  */

		printf("Planning forward FFT...");
		fft_occs=fftw_plan_dft_r2c_2d(dim,dim,temp_occs,fftout,FFTW_MEASURE); /* Plan forward transform */
		printf("done!\n");
	}

	/* Normal mode - writes output to appropriate time step position in storage matrix */
	for(i=0;i<n_sites;i++) {
		temp_occs[i]=(double)occupations[i];	/* Populate temp. occs. array */
	}
	fftw_execute(fft_occs);		/* Run */

	if( time_index>=0 )	/* If storage index set, transfer result to fft_steps matrix */
		for(i=0;i<n_sites/2+dim;i++) {
			fft_steps[time_index*(n_sites/2+dim)+i][0]=fftout[i][0];
			fft_steps[time_index*(n_sites/2+dim)+i][1]=fftout[i][1];
		}
}

void fft_mult ( int index1, int index2 ) /* Multiplies two FFTs together */
{
	int i,offset1,offset2;
	double step1_re,step1_im,step2_re,step2_im;
		
	offset1=index1*(n_sites/2+dim); /* Set starting offsets in fft_steps matrix for data for the timesteps */
	offset2=index2*(n_sites/2+dim);
	
	for(i=0;i<dim/2+1;i++) {	/* Perform multiplication */
			step1_re=(double)fft_steps[offset1+i][0];	/* Take elements from lower timestep */
			step1_im=(double)fft_steps[offset1+i][1];

			step2_re=(double)fft_steps[offset2+i][0];	/* Take elements from higher timestep */
			step2_im=(double)(-1*fft_steps[offset2+i][1]);	/* Take complex conjugate */
				
			/* Multiply, calculate the real and imaginary parts separately */
			temp_ISF[i][0]=step1_re*step2_re-step1_im*step2_im;
			temp_ISF[i][1]=step1_im*step2_re+step1_re*step2_im;
	}
}

void display_data ( int mode )	/* Display occupations matrix or FFT */
{
	int lineindex,i,k;

	if( !mode ) {
		printf("+ ");
		for(i=0;i<dim;i++) printf("- ");
		printf("+\n");
		for(i=0;i<n_sites;i++) {
			if( !((i)&(dim-1)) ) {
				lineindex=i/dim;
				for(k=0;k<lineindex+1;k++) printf(" ");
				printf("\\ ");
			}
			if( occupations[i] ) printf("1 ");
			else printf("  ");
			if( !((i+1)&(dim-1)) ) printf("\\\n");
		}
		for(i=0;i<dim+1;i++) printf(" ");
		printf("+ ");
		for(i=0;i<dim;i++) printf("- ");
		printf("+\n");
	}
	else	for(i=0;i<n_sites/2+dim;i++) {
				printf("%+.3f %+2.3fi\t",fftout[i][0],fftout[i][1]);
				if( !((i+1)%(dim/2+1)) ) printf("\n");
		}
}
