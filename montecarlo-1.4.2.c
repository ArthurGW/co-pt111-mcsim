/*
Monte-Carlo simulation of a CO/Pt(111) surface system to replicate observed dynamics.
Written by Arthur Gordon-Wright, https://github.com/ArthurGW/						
 																							
Version number 1.4.2, date 04/01/19																
Update notes: Various fixes to run on a different compiler, especially random()->rand() issues
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/timeb.h>
#include <unistd.h>
#include <sys/stat.h>
#include "fftw3.h"

#define ONEOVSQRT2 0.707106781		/* One over the square root of 2, to save repeated calculation */
#define SQRT1OV8 0.353553391		/* Square root of 1 over 8 */
#define UTT 210		/* Top-top diffusion barrier in meV */
#define USAGE "Usage: %s dim n_co int_mode out_name\n"
#define OPTIONS "g p s tfact TFACT tsteps TSTEPS repeats REPEATS\n"
#define LINEWIDTH 22		/* Number of '-' separators in output file dividers */
#define OUT_NAME_LEN 100	/* Maximum length of output file name (and specified path) */
#define MAX_OUTSIZE 1048576	/* Maximum output file size (1Mb) */
#define GRES 100	/* Resolution for pair correlation function - number of bins per lattice vector */
#define BETA 0.0341308965	/* 1/kT in meV^-1 at 340K */
#define ISF_VECT_SCL 1.960985419	/* Factor used to calculate length of <11-2> inverse azimuth */
#define TIMESCALE 0.16693	/* Time of each step in ns */

/* Variable declarations */
char name[30], out_name[OUT_NAME_LEN], outmode='w';
int dim, n_sites, n_co, int_mode, skip_sweeps=100000, g=0, s=0, p=0, e=0, bins, time_steps=8, time_factor=2, repeats=1000000, total_sweeps=0, count=0, rep=10, energy;
time_t start_time;
double bin_res;
fftw_complex *fftout, *fft_steps, *ISF_av, *temp_ISF;
double *temp_occs;
int *occupations, *adatoms, *G;

/* Function declarations */
int args ( int , char ** );
void init_rand ( );
int init_adatoms ( );
void error_close ( int );
int output ( int );
int test_file ( );
double interaction ( int );
char * int_name ( );
void sweep ( );
void pair ( );
void fft_run ( int, signed int );
void fft_mult ( int, int );
void display_data ( int );
void pair ( );
double separation ( int, int );
int total_energy ( );
int randint ( int );

int randint ( int nmax ) {
	// Generate a random int with 0 <= result < nmax
	int rn;
	double fraction;

	// Ensure rn != RAND_MAX so result != nmax
	do {
		rn = rand();
	} while (rn == RAND_MAX);

	fraction = (double)rn / RAND_MAX;

	return (int) floor(nmax * fraction);
}

int main ( int argc, char *argv[] )
{
	int err_code, i, j, k, lineindex, bin, time_index, repeat=0, sweeps, energy_sweeps=0, next_step;
	double *pair_av;
	char en_out_name[20];
	FILE *en_out;

	if( argc==1 || strcmp(argv[1],"help")==0 ) { /* Test if help has been selected or no arguments set */
		printf(USAGE, argv[0]);
		return 0;
	}

	time( &start_time );	/* Get current time */

	init_rand();	/* Set random seed */

	if( (err_code=args( argc, argv )) ) { /* Get and test command-line arguments, exit on failure */
		error_close(err_code);
	}

	if( ( err_code=init_adatoms() ) ) {	/* Initialise adatom positions, exit on failure */
		error_close(err_code);
	}

	interaction(7);	/* Set up interaction energy lookup table */
	if( e ) {
		energy=total_energy();
		sprintf(en_out_name,"%d-%.3f-%d.en",dim,(float)n_co/n_sites,int_mode?rep:0);
		en_out=fopen(en_out_name,"w");
		fprintf(en_out,"0\t%f\n",energy*rep/1000.);
		fclose(en_out);		
	}
	
	output(1); output(2);	/* Output headers and initial positions to temp. file */

	printf("name=%s\ndim=%d, n_co=%d, int_mode=%d, rep=%d\nskip_sweeps=%d, time_factor=%d, time_steps=%d, repeats=%d\nout_name=%s\n", name, dim, n_co, int_mode, int_mode?rep:0, skip_sweeps, time_factor, time_steps, repeats, out_name );	/* Display some parameters */
	printf("options: %s%s%s\n",g?"g":" ",p?"p":" ",s?"s":" ");	/* Display options */

	printf("\nOccupations.\n");	/* Display current occupations matrix */
	display_data(0);
	
	if( e ) printf("Initial energy: %f\n",energy*rep/1000.);

	for(i=0;i<skip_sweeps;i++) {	/* Skip some sweeps */
		sweep();
		if( e && i<((10000/rep)<100?100:(10000/rep)) ) {
			en_out=fopen(en_out_name,"a");
			fprintf(en_out,"%d\t%f\n",++energy_sweeps,energy*rep/1000.);
			fclose(en_out);
		}
	}

	/* FFT - initialise input and output arrays, then run (with planning) */
	temp_occs = (double*) fftw_malloc( sizeof(double) * n_sites );	/* Copy of occupations matrix */
	fftout=(fftw_complex *) fftw_malloc( 2 * sizeof(double) * n_sites );	/* FFT output array */
	fft_steps=(fftw_complex *) fftw_malloc( (time_steps+1) * 2 * sizeof(double) * (n_sites/2+dim) );
							/* Array for storing needed bits of FFTs for each time step */
	fft_run(1,time_steps);	/* Plan and run FFT */

	output(4);	/* Output FFT result */

	ISF_av=(fftw_complex *) fftw_malloc( 2 * sizeof(double) * dim * (time_steps+1) );
								/* Array for storing ISF means and variances */
	pair_av=(double *) malloc( sizeof(double) * dim/2 * (time_steps+1) );
								/* Array for storing ISF means and variances */
	
	for(i=0;i<dim*(time_steps+1);i++) {	/* Set ISF_av to zero intially */
		ISF_av[i][0]=0.;
		ISF_av[i][1]=0.;
	}
	for(i=0;i<(dim/2)*(time_steps+1);i++) {	/* Set ISF_av to zero intially */
		pair_av[i]=0.;
	}
	
	temp_ISF=(fftw_complex *) fftw_malloc( 2 * sizeof(double) * dim/2 );	/* Array for storing each ISF result */

	if( p ) {
		bins=GRES*dim;	/* Set up pair correlation function: number of bins, bin width, counting array */
		bin_res=ONEOVSQRT2/GRES;
		G=(int *)malloc(bins*sizeof(int));

		for(i=0;i<bins;i++) G[i]=0;	/* Set all separation counts to zero */
	}

	do {	/* Main loop (over repeats) */
		time_index=1;	/* Reset time index, sweeps index */
		sweeps=0;
		next_step=time_factor;	/* Calculate next sweep index to record at */

		for(i=0;i<dim/2;i++) {	/* Set t=0 FFT record to current FFT result */
			fft_steps[i][0]=fft_steps[time_steps*dim/2+i][0];
			fft_steps[i][1]=fft_steps[time_steps*dim/2+i][1];
		}

		do {	/* Loop over time steps */
			sweep();	/* Run sweep */
			sweeps++;
			total_sweeps++;
			
			if( 0 && !((++energy_sweeps)&1023) ) {	/* Output total energy if in energy mode every 1024 sweeps */
				en_out=fopen(en_out_name,"a");
				fprintf(en_out,"%d\t%f\n",energy_sweeps,energy*rep/1000.);
				fclose(en_out);
			}

			if( p ) pair();	/* Calculate pair correlation function */

			if( g ) {	/* If in graphical mode, display occupations matrix */
				system("clear");
				printf("Sweep #%d:\n",sweeps);
				display_data(0);
				sleep(1);
			}

			if( sweeps==next_step ) {	/* If sweep is to be recorded */
				fft_run(0,time_index);	/* Calculate FFT, store in fft_steps */
				time_index++;		/* Calculate next sweep to be recorded */
				next_step=s?time_factor*time_index:pow(time_factor,time_index);
			}

		} while(time_index<time_steps+1);
		
		for(time_index=0;time_index<time_steps+1;time_index++) {
			fft_mult(0,time_index);	/* Multiply FFTs for t=time_index and t=0 */
			for(i=0;i<dim/2;i++) {	/* Transfer results from temp. array to ISF average array */
				ISF_av[time_index*dim+2*i][0]+=temp_ISF[i][0];
				ISF_av[time_index*dim+2*i][1]+=temp_ISF[i][1];
				ISF_av[time_index*dim+2*i+1][0]+=(temp_ISF[i][0]*temp_ISF[i][0]);
				ISF_av[time_index*dim+2*i+1][1]+=(temp_ISF[i][1]*temp_ISF[i][1]);
			}		
		}

		repeat++;
		if( !(repeat&131071) ) {
			printf("Repeat no.: %d\n",repeat);	/* Display repeat number occasionally */
		}

	} while(repeat<repeats);

	fftw_free(temp_ISF);
	fftw_free(fft_steps);
	free(adatoms);
	fftw_free(temp_occs); /* Free arrays that are no longer needed */

	// Calculate means and variances
	for(time_index=0;time_index<time_steps+1;time_index++) {
		for(i=0;i<dim;i+=2) {
			// Normalise mean and mean-squared by number of repeats
			ISF_av[time_index*dim+i][0] /= repeats;
			ISF_av[time_index*dim+i][1] /= repeats;
			ISF_av[time_index*dim+i+1][0] /= repeats;
			ISF_av[time_index*dim+i+1][1] /= repeats;

			// Subract (mean)^2 from mean-squared
			ISF_av[time_index*dim+i+1][0]-=
					ISF_av[time_index*dim+i][0]*ISF_av[time_index*dim+i][0];
			ISF_av[time_index*dim+i+1][1]-=
					ISF_av[time_index*dim+i][1]*ISF_av[time_index*dim+i][1];
		}		
	}
	
	if( e ) printf("Final energy: %f (Total energy: %f)\n",energy*rep/1000.,total_energy()*rep/1000.);

	/*printf("\nFinal occupations.\n");	/* Display final positions */
	/*display_data(0);*/

	/*printf("\nTransformed data.\n");	/* Display final FFT */
	/*display_data(1);*/

	output(2); output(4);	/* Output final positions and FFT*/

	output(6);	/* Output ISF data*/
	
	if( p ) {
		output(10);	/* Output pair correlation function*/
		free(G);
	}

	if( ( err_code=output(11) ) ) {	 /* Write simulation results to file, exit on error */
		error_close(err_code);
	}

	fftw_free(ISF_av);
	fftw_free(fftout);
	free(occupations);
	
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
		printf("Specified coverage is greater than 100%%!\n");
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
		else if( !strcmp(argv[i],"s") ) {
			s=1;
		}
		else if( !strcmp(argv[i],"p") ) {
			p=1;
		}
		else if( !strcmp(argv[i],"e") ) {
			e=1;
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
		else if( !strcmp(argv[i],"rep") && argv[i+1] ) {
			rep=atoi(argv[i+1]);
			if( !repeats ) {
				printf("Please enter an integer value for optional parameter rep.\n");
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

void init_rand ( )	/* Initialise random number seed */
{
	struct timeb tp;

	ftime(&tp);
	srand( (tp.time*1000)+tp.millitm ); /* Use time in milliseconds to ensure no repeats */
}

int init_adatoms ( )	/* Initialise adatom and occupation arrays, placing atoms at random positions */
{
	int i, unallocated=n_sites, index;
	int free_sites[n_sites];
	double tmp;

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

	printf("Positioning adatoms:\n");
	for(i=0;i<n_co;i++) {
		index = randint(unallocated); /* Choose random free site */
		printf("Chosen site for adatom %d: free index: %d site: %d\n", i, index, free_sites[index]);
		adatoms[i]=free_sites[index];	/* Place adatom at site */
		occupations[ adatoms[i] ]=1;	/* Set occupation state for chosen coordinate */
		free_sites[index]=free_sites[unallocated - 1];/* Replace free site coordinate with highest unalloc. one */
		unallocated--;	/* Decrement the number of remaining free sites */
	}

	return 0;
}

int output ( int out_num )	/* Writes data to output file.  Attempts writing three times, then asks user to	*/	
{							/* try a different file, then another, then fails 						*/
	int i, j, n, attempts=1, success=0, file_no=1, normalisation;
	FILE *tmp, *out, *ISF_out, *param_out;
	static char tmp_name[15];
	char *  buffer;
	char ISF_out_name[20], param_out_name[20];

	if( out_num==1 ) {	/* Create and display temp. file name on first run */
		sprintf(tmp_name,"%d.tmp",(int)start_time);
		printf("tmp_name=%s\n",tmp_name);
	}

	tmp=fopen(tmp_name,"a");	/* Open temp file */

	/* Write output to temp file. out_num flag sets which stage output to write */	
	if( out_num==1 ) {	/* Header section */

		fprintf(tmp,"CO/Pt111 System Monte-Carlo Simulation.\t\tStarted: %s\n\n",ctime( &start_time ) );

		fprintf(tmp,"name=%s\ndim=%d, n_co=%d, rep=%d\nskip_sweeps=%d, time_factor=%d, time_steps=%d, repeats=%d\n", name, dim, n_co, int_mode?rep:0, skip_sweeps, time_factor, time_steps, repeats );
		fprintf(tmp,"options: %s%s\n",s?"s":" ",p?"p":" ");
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

		for(i=0;i<dim/2*(time_steps+1);i++) {
			fprintf(tmp,"%.3f\t%.3f\t\t",fft_steps[i][0],fft_steps[i][1]);
			if( !((i+1)&(dim/2-1)) ) fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}
	else if( out_num==6 ) {	/* ISF real and imaginary means and variances */
		fclose(tmp);
		double vect_len=ISF_VECT_SCL*(1-2/dim), error_fact=repeats>1?sqrt(401/(repeats-1)):1000000;
		int count[dim/2];
		for(i=0;i<dim/2;i++) count[i]=0;
		
		sprintf(ISF_out_name,"%d-%.3f-%d.ISF",dim,(float)n_co/n_sites,int_mode?rep:0);
		sprintf(param_out_name,"%d-%.3f-%d.prm",dim,(float)n_co/n_sites,int_mode?rep:0);

		ISF_out=fopen(ISF_out_name,"w");
		for(j=0;j<(time_steps+1);j++) {
			fprintf(ISF_out,"%.5f\t",s?TIMESCALE*j*time_factor:TIMESCALE*pow(time_factor,j));
			for(i=0;i<dim;i+=2) {
				if( log(ISF_av[j*dim+i][0])>=(log(ISF_av[i][0])-4) ) {
					fprintf(ISF_out,"%.5f",log(ISF_av[j*dim+i][0]));
					count[i/2]++;
				}
				fprintf(ISF_out,"\t");
			}
			fprintf(tmp,"\n");
		}
		fclose(ISF_out);

		param_out=fopen(param_out_name,"w");
		for(i=0;i<dim/2;i++) {
			fprintf(param_out,"m%d = %.5f\nc%d = %.5f\n",i,
					log(ISF_av[i][0]/ISF_av[(count[i]-1)*dim+i][0])/(s?TIMESCALE*(count[i]-1)*time_factor:TIMESCALE*pow(time_factor,(count[i]-1))),
					i,log(ISF_av[i][0]));
		}
		fclose(param_out);

		tmp=fopen(tmp_name,"a");

		fprintf(tmp,"ISF along <11-2> azimuth:\n");
		fprintf(tmp,"time_factor: %d, time_steps: %d, repeats: %d\n\n",time_factor,time_steps,repeats);
		fprintf(tmp,"Real\nTimestep 1\tTimestep 2\t");
		for(i=0;i<dim/2;i++) fprintf(tmp,"%.5f\t",vect_len*i/(dim/2-1));
		fprintf(tmp,"\n");
		for(j=0;j<(time_steps+1);j++) {
			fprintf(tmp,"%d\t%d\t",0,s?(int)j*time_factor:(int)pow(time_factor,j));
			for(i=0;i<dim;i+=2)
				fprintf(tmp,"%.5f\t",ISF_av[j*dim+i][0]);
			fprintf(tmp,"\n");
		}

		fprintf(tmp,"\nImaginary\nTimestep 1\tTimestep 2\t");
		for(i=0;i<dim/2;i++) fprintf(tmp,"%.5f\t",vect_len*i/(dim/2-1));
		fprintf(tmp,"\n");
		for(j=0;j<(time_steps+1);j++) {
			fprintf(tmp,"%d\t%d\t",0,s?(int)j*time_factor:(int)pow(time_factor,j));
			for(i=0;i<dim;i+=2)
				fprintf(tmp,"%.5f\t",ISF_av[j*dim+i][1]);
			fprintf(tmp,"\n");
		}

		fprintf(tmp,"\nReal variance\nTimestep 1\tTimestep 2\t");
		for(i=0;i<dim/2;i++) fprintf(tmp,"%.5f\t",vect_len*i/(dim/2-1));
		fprintf(tmp,"\n");
		for(j=0;j<(time_steps+1);j++) {
			fprintf(tmp,"%d\t%d\t",0,s?(int)j*time_factor:(int)pow(time_factor,j));
			for(i=0;i<dim;i+=2)
				fprintf(tmp,"%.5f\t",error_fact*sqrt(ISF_av[j*dim+i+1][0]));
			fprintf(tmp,"\n");
		}
		fprintf(tmp,"\nImaginary variance\nTimestep 1\tTimestep 2\t");
		for(i=0;i<dim/2;i++) fprintf(tmp,"%.5f\t",vect_len*i/(dim/2-1));
		fprintf(tmp,"\n");
		for(j=0;j<(time_steps+1);j++) {
			fprintf(tmp,"%d\t%d\t",0,s?(int)j*time_factor:(int)pow(time_factor,j));
			for(i=0;i<dim;i+=2)
				fprintf(tmp,"%.5f\t",error_fact*sqrt(ISF_av[j*dim+i+1][1]));
			fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");
	}
	else if( out_num==10 ) {	/* Closing section */
		fprintf(tmp,"Final coordinates:\n");

		for(i=0;i<n_sites;i++) {
			fprintf(tmp,"%d ",occupations[i]);
			if( !((i+1)&(dim-1)) )
				fprintf(tmp,"\n");
		}

		fprintf(tmp,"\n");
		for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
		fprintf(tmp,"\n\n");

		if( p ) {
			normalisation=total_sweeps*n_co*(n_co-1);
			fprintf(tmp,"Pair correlation function.\n");
			fprintf(tmp,"Bin resolution: %1.9f, normalisation factor: %d\n\n",bin_res,normalisation);
			fprintf(tmp,"Bin centre\tProbability\tNumber in bin\n");
			for(i=0;i<GRES*(dim/2)+1;i++)
				if( G[i] )
					fprintf(tmp,"%2.4f\t%1.6f\t%d\n",(double)i/GRES,(double)G[i]/normalisation,G[i]);
	
			fprintf(tmp,"\n");
			for(i=0;i<LINEWIDTH;i++) fprintf(tmp,"-");
			fprintf(tmp,"\n\n");
		}
	}


	fclose(tmp);

	if( out_num==11 ) {	/* Transfer temporary output to specified file */ 
		while( !success ) {
			/* Open tmp for reading and out for writing */
			if( (tmp=fopen(tmp_name,"r")) && (out=fopen(out_name, &outmode)) ) {
				buffer=(char *) malloc(MAX_OUTSIZE);/* Read temp to buffer then write to output */
				if ( buffer != NULL ) {
					n=fread(buffer,1,MAX_OUTSIZE,tmp);

					printf("Read bytes: %d\n",n);
					fwrite(buffer,1,n,out);
					fputc('\n',out);

					if (!ferror(tmp) && !ferror(out)) success = 1;

					free(buffer);	/* Free buffer */
				}
				fclose(tmp);	/* Close files */
				fclose(out);
			}

			if( !success ) {
				printf("Attempt #%d to write output file failed.\n",attempts);

				attempts++;

				if( attempts==4 ) {	/* After three attempts at using a file, try another, 	*/
					file_no++;	/* assuming less than three have been tried 		*/
					if( file_no<4 ) {
						attempts=1;
						do {	/* Keep asking for a new file until a correct one is specified */
							printf("Cannot write to output file.  Please specify another:\n");
							scanf("%s",(char *)&out_name);
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

double interaction ( int delta_nn )
{
	static double lookup[3];
	int i;

	if( delta_nn>3 ) {	/* Set up lookup table */
		for(i=0;i<3;i++) {
			lookup[i]=exp(-BETA*(i+1)*rep);
		}
		return 0.;
	}
	else if( delta_nn<=0 ) return 1.;
	else return lookup[delta_nn-1];	
}

void sweep ( )
{
	int i,j,atom,coord,nn[6],nns,mnns;
	double rand_factor;

	for(i=0;i<n_co;i++) {	/* For each sweep, iterate over n_co moves */

		/* Select random adatom */
		atom=randint(n_co);
		coord=adatoms[atom];

		/* For initial coord., calculate nearest neighbour coords. using helical BCs */
		nn[0]=(coord+1)&(n_sites-1);
		nn[1]=(coord-1+n_sites)&(n_sites-1);
		nn[2]=(coord+dim)&(n_sites-1);
		nn[3]=(coord-dim+n_sites)&(n_sites-1);
		nn[4]=(coord+(dim+1))&(n_sites-1);
		nn[5]=(coord-(dim+1)+n_sites)&(n_sites-1);

		coord=randint(6);	/* Select random move */

		if( occupations[ nn[coord] ] ) continue;	/* If no move selected, go to next move */
		else coord=nn[coord];

		nns=0; mnns=0;	/* Reset numbers of neighbours */

		/* Calculate current number of occupied neighbours */
		for(j=0;j<6;j++) nns+=occupations[ nn[j] ];

		/* For move coord., calculate nearest neighbour coords using helical BCs */
		nn[0]=(coord+1)&(n_sites-1);
		nn[1]=(coord-1+n_sites)&(n_sites-1);
		nn[2]=(coord+dim)&(n_sites-1);
		nn[3]=(coord-dim+n_sites)&(n_sites-1);
		nn[4]=(coord+(dim+1))&(n_sites-1);
		nn[5]=(coord-(dim+1)+n_sites)&(n_sites-1);

		/* Calculate potential new number of occupied neighbours */
		for(j=0;j<6;j++) mnns+=occupations[ nn[j] ];
		mnns--;	/* Reduce by one due to vacated space */

		rand_factor=(double)rand()/RAND_MAX;	/* Calculate random fraction */

		/* Simple test interaction model */
		if( !int_mode || mnns<=nns || ( rand_factor <= interaction(mnns-nns) ) ) {
			occupations[adatoms[atom]]=0;	/* Move adatom */
			occupations[coord]=1;
			adatoms[atom]=coord;
			energy+=(mnns-nns);	/* Update total energy */
		}
	}
}

char * int_name ( )
{
	if( int_mode==0 ) return "None";
	if( int_mode==1 ) return "Nearest neighbour";
	if( int_mode==2 ) return "Petrova";
	if( int_mode==3 ) return "Persson";
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
		for(i=0;i<dim/2;i++) {
			fft_steps[time_index*dim/2+i][0]=fftout[2*i*(dim/2+1)+i][0];
			fft_steps[time_index*dim/2+i][1]=fftout[2*i*(dim/2+1)+i][1];
		}
}

void fft_mult ( int index1, int index2 ) /* Multiplies two FFTs together */
{
	int i,offset1,offset2;
	double step1_re,step1_im,step2_re,step2_im;
		
	offset1=index1*dim/2; /* Set starting offsets in fft_steps matrix for data for the timesteps */
	offset2=index2*dim/2;
	
	for(i=0;i<dim/2;i++) {	/* Perform multiplication */
			step1_re=(double)fft_steps[offset1+i][0];	/* Take elements from lower timestep */
			step1_im=(double)fft_steps[offset1+i][1];

			step2_re=(double)fft_steps[offset2+i][0];	/* Take elements from higher timestep */
			step2_im=(double)(-1*fft_steps[offset2+i][1]);	/* Take complex conjugate */
				
			/* Multiply, calculate the real and imaginary parts separately */
			temp_ISF[i][0]=(step1_re*step2_re-step1_im*step2_im)/n_co;
			temp_ISF[i][1]=(step1_im*step2_re+step1_re*step2_im)/n_co;
	}
}

void display_data ( int mode )	/* Display occupations matrix or FFT */
{
	int lineindex,i,k;

	if( !mode ) {
		for(i=0;i<dim+1;i++) printf(" ");
		printf("+ ");
		for(i=0;i<dim;i++) printf("- ");
		printf("+\n");
		for(i=0;i<n_sites;i++) {
			if( !((i)&(dim-1)) ) {
				lineindex=i/dim;
				for(k=0;k<dim-lineindex;k++) printf(" ");
				printf("/ ");
			}
			if( occupations[i] ) printf("1 ");
			else printf("  ");
			if( !((i+1)&(dim-1)) ) printf("/\n");
		}
		printf("+ ");
		for(i=0;i<dim;i++) printf("- ");
		printf("+\n");
	}
	else	for(i=0;i<n_sites/2+dim;i++) {
				printf("%+.3f %+2.3fi\t",fftout[i][0],fftout[i][1]);
				if( !((i+1)%(dim/2+1)) ) printf("\n");
		}
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

double separation ( int coord1, int coord2 )
{
	unsigned int i, num_mirrors, sep, new_sep; 
	signed int twodx, twodx_m, dy, dy_m;

	static signed int mirrors[8];
	signed int current_mirrors[6]={0,0,0,0,0,0};
	static int flag=0;

	if( !flag ) {
		mirrors[0]=-1*dim;
		mirrors[1]=2*dim;
		mirrors[2]=dim+1;
		mirrors[3]=2*dim-2;
		mirrors[4]=2*dim+1;
		mirrors[5]=-2;
		mirrors[6]=3*dim+1;
		mirrors[7]=-2*dim-2;
		flag++;
	}

	dy=(int)floor(coord2/dim) - (int)floor(coord1/dim);
	twodx=2 * ( (coord2&(dim-1)) - (coord1&(dim-1)) ) - dy;

	sep=(twodx*twodx)+(3*dy*dy);

	if( twodx+dy>0 ) {
		if( dy>1 ) {
			num_mirrors=3;
			current_mirrors[0]=-1*mirrors[0];
			current_mirrors[1]=-1*mirrors[1];
			current_mirrors[2]=-1*mirrors[2];
			current_mirrors[3]=-1*mirrors[3];
			current_mirrors[4]=-1*mirrors[4];
			current_mirrors[5]=-1*mirrors[5];
		}
		else if( dy<0 ) {
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
	else if( twodx+dy<0 ) {
		if( dy>0 ) {
			num_mirrors=3;
			current_mirrors[0]=-1*mirrors[0];
			current_mirrors[1]=-1*mirrors[1];
			current_mirrors[2]=mirrors[4];
			current_mirrors[3]=mirrors[5];
			current_mirrors[4]=mirrors[6];
			current_mirrors[5]=mirrors[7];
		}
		else if( dy<-1 ) {
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

int total_energy ( )
{
	int i, j, nn[6], bonds=0;

	for(i=0;i<n_sites;i++)
		if(occupations[i]) {
			nn[0]=(i+1)&(n_sites-1);
			nn[1]=(i+dim)&(n_sites-1);
			nn[2]=(i+(dim+1))&(n_sites-1);

			for(j=0;j<3;j++) bonds+=occupations[ nn[j] ];
		}
		
	return bonds;
}
