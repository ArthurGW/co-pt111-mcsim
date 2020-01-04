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
	display_occs();

	/* FFT - initialise input and output arrays, then run (with planning) */
	temp_occs = (double*) fftw_malloc( sizeof(double) * n_sites );	/* Copy of occupations matrix */
	fftout=(fftw_complex *) fftw_malloc( 2 * sizeof(double) * n_sites );	/* FFT output matrix */
	fftout_full=(fftw_complex *) malloc( 2 * sizeof(double) * n_sites );	/* Extended FFT output matrix */
	fft_steps=(fftw_complex *) malloc( (time_steps+1) * 2 * sizeof(double) * n_sites );
								/* Matrix for storing FFTs for each time step */
	ISF_av=(fftw_complex *) malloc( 2 * sizeof(double) * 2 * dim * (time_steps+1) );
								/* Matrix for storing ISF means and variances */

	for(i=0;i<skip_sweeps;i++) sweep();	/* Skip some sweeps */

	fft_run(1,time_steps);	/* Plan and run FFT */

	printf("\nTransformed data.\n");	/* Display FFT result */
	for(i=0;i<n_sites;i++) {
				printf("%+.3f %+2.3fi\t",fftout_full[i][0],fftout_full[i][1]);
				if( !((i+1)&(dim-1)) ) printf("\n");
	}

	output(4);	/* Output FFT result */

	do {	/* Main loop (over repeats) */
		time_index=1;	/* Reset time index, sweeps index */
		sweeps=0;
		next_step=pow(time_factor,1);	/* Calculate next sweep index to record at */

		for(i=0;i<n_sites;i++) {	/* Set t=0 FFT record to current FFT result */
			fft_steps[i][0]=fft_steps[time_steps*n_sites+i][0];
			fft_steps[i][1]=fft_steps[time_steps*n_sites+i][1];
		}

		do {	/* Loop over time steps */
			sweep();	/* Run sweep */
			sweeps++;

			if( g ) {	/* If in graphical mode, display occupations matrix */
				system("clear");
				printf("Sweep #%d:\n",sweeps);
				display_occs();
				sleep(1);
			}

			if( sweeps==next_step ) {	/* If sweep is to be recorded */
				fft_run(0,time_index);	/* Calculate FFT, store in fft_steps */
				time_index++;		/* Calculate next sweep to be recorded */
				next_step=pow(time_factor,time_index);
			}

		} while(time_index<time_steps+1);
		
		temp_ISF=(fftw_complex *) fftw_malloc( 2 * sizeof(double) * n_sites );	/* Create temporary array */
		for(time_index=0;time_index<time_steps+1;time_index++) {
			fft_mult(0,time_index);	/* Multiply FFTs for t=time_index and t=0 */
			for(i=0;i<dim;i++) {	/* Transfer results from temp. array to ISF average array */
				ISF_av[time_index*(2*dim)+2*i][0]+=temp_ISF[i][0]/repeats;
				ISF_av[time_index*(2*dim)+2*i][1]+=temp_ISF[i][1]/repeats;
				ISF_av[time_index*(2*dim)+2*i+1][0]+=temp_ISF[i][0]*temp_ISF[i][0]/repeats;
				ISF_av[time_index*(2*dim)+2*i+1][1]+=temp_ISF[i][1]*temp_ISF[i][1]/repeats;
			}		
		}
		free(temp_ISF);	/* Free temp. array */

		repeat++;
		if( !(repeat&8191) ) printf("Repeat no.: %d\n",repeat);	/* Display repeat number occasionally */

	} while(repeat<repeats);

	for(time_index=0;time_index<time_steps+1;time_index++) {
		for(i=0;i<dim;i++) {	/* Calculate variance by subtracting mean^2 from mean squared value */
			ISF_av[time_index*(2*dim)+2*i+1][0]-=
					ISF_av[time_index*(2*dim)+2*i][0]*ISF_av[time_index*(2*dim)+2*i][0];
			ISF_av[time_index*(2*dim)+2*i+1][1]-=
					ISF_av[time_index*(2*dim)+2*i][1]*ISF_av[time_index*(2*dim)+2*i][1];
		}		
	}

	printf("\nFinal occupations.\n");	/* Display final positions */
	display_occs();

	printf("\nTransformed data.\n");	/* Display final FFT */
	for(i=0;i<n_sites;i++) {
				printf("%+.3f %+2.3fi\t",fftout_full[i][0],fftout_full[i][1]);
				if( !((i+1)&(dim-1)) ) printf("\n");
	}

	output(2); output(4);	/* Output final positions and FFT */

	output(6);	/* Output ISF data */

	free_arrays();	/* Free arrays */

	if( ( err_code=output(11) ) ) {	 /* Write simulation results to file, exit on error */
		error_close(err_code);
	}

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
		if( !int_mode || mnns<=nns || (random()/RAND_MAX) < (5+nns-mnns)/6 ) {
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
	int i,j;
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

	/* Populate full FFT output array for returning */
	for(i=dim-1;i>=0;i--) {
		for(j=dim/2;j>=0;j--) {	/* Position +k points */
			fftout_full[i*dim+(dim/2-1)+j][0]=fftout[i*(dim/2+1)+j][0];
			fftout_full[i*dim+(dim/2-1)+j][1]=fftout[i*(dim/2+1)+j][1];
		}
		for(j=1;j<dim/2;j++) {	/* Set -k points using complex conjugates of +k */
			fftout_full[i*dim+(dim/2-1)-j][0]=fftout_full[i*dim+(dim/2-1)+j][0];
			fftout_full[i*dim+(dim/2-1)-j][1]=-1*fftout_full[i*dim+(dim/2-1)+j][1];
		}
	}

	if( time_index>=0 )	/* If storage index set, transfer result to fft_steps matrix */
		for(i=0;i<n_sites;i++) {
			fft_steps[time_index*n_sites+i][0]=fftout_full[i][0];
			fft_steps[time_index*n_sites+i][1]=fftout_full[i][1];
		}
}

fftw_complex * fft_mult ( int index1, int index2 ) /* Multiplies two FFTs together */
{
	int i,j,k,offset1,offset2,mult_row;
	fftw_complex *step1, *step2;
	double step1_re,step1_im,step2_re,step2_im;
	
	/* Initialise temporary matrices to store the two timesteps */
	step1=(fftw_complex *) calloc(n_sites,2*sizeof(double));
	step2=(fftw_complex *) calloc(n_sites,2*sizeof(double));
	
	offset1=index1*n_sites; /* Set starting offsets in fft_steps matrix for data for the timesteps */
	offset2=index2*n_sites;
	
	for(i=0;i<n_sites;i+=dim)	/* Populate temporary matrices */
		for(j=0;j<dim;j++) {
			step1[i+j][0]=fft_steps[offset1+i+j][0];
			step1[i+j][1]=fft_steps[offset1+i+j][1];
			
			step2[i+j][0]=fft_steps[offset2+i+j][0]; /* Take complex conjugate for higher time */
			step2[i+j][1]=-1*fft_steps[offset2+i+j][1];
		}
	
	for(i=0;i<n_sites;i+=dim)	/* Perform multiplication */
		for(j=0;j<dim;j++) {
			temp_ISF[i+j][0]=0.;	/* Clear results matrix */
			temp_ISF[i+j][1]=0.;
			
			for(k=0;k<dim;k++) {	/* Index to iterate over elements in relevant row and column */
				step1_re=(double)step1[k*dim+j][0];/* For lower timestep, take elements from column j */
				step1_im=(double)step1[k*dim+j][1];

				step2_re=(double)step2[i+k][0];	/* For higher timestep, take elements from row i */
				step2_im=(double)step2[i+k][1];
				
				/* Multiply, and add the real and imaginary parts separately to totals */
				temp_ISF[i+j][0]+=step1_re*step2_re-step1_im*step2_im;
				temp_ISF[i+j][1]+=step1_im*step2_re+step1_re*step2_im;
			}
		}
	
	free(step1);	/* Free temporary arrays */
	free(step2);

	return temp_ISF;
}
