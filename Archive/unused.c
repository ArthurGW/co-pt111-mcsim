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
				
			ijsep=separation1(adatoms[i],adatoms[j]);
			bin=(int)floor(.5+ijsep*GRES/ONEOVSQRT2);

			G[bin]++;
			
			j++;
		}
	}
}

int ** xycoords ( )
{
	int i,**coords;

	coords=(int *) malloc(sizeof(int)*n_co*2);

	for(i=0;i<n_co;i++) {
		coords[i][0]=adatoms[i]%(dim-1);
		coords[i][1]=(int)floor(adatoms[i]/dim);
	}

	return coords;
}

int ** twod_occs ( )
{
	int i,j,**occs;

	occs=(int **) malloc(sizeof(int)*n_sites);

	for(i=0;i<dim;i++)
		for(j=0;j<dim;j++)
			occs[i][j]=occupations[i*dim+j];

	return occs;
}

double separation1 ( int coord1, int coord2 )
{
	unsigned int i, j, k, num_mirrors, sep, new_sep, count=1; 
	signed int twodx, twodx_m, dy, dy_m;

	signed int mirrors[8];
	signed int current_mirrors[6]={0,0,0,0,0,0};

	if( !coord1 & !coord2 ) {
		printf("Creating separation lookup table...");

		sep_lookup=(double *) fftw_malloc( (6*dim-5) * dim * sizeof(double) );

		mirrors[0]=dim;
		mirrors[1]=2*dim;
		mirrors[2]=3*dim-1;
		mirrors[3]=2*dim-2;
		mirrors[4]=2*dim-1;
		mirrors[5]=-2;
		mirrors[6]=dim-1;
		mirrors[7]=-2*dim-2;		

		for(i=0;i<6*dim-5;i++) {
			for(j=0;j<dim;j++){
			
				twodx=2*(-1.5*(dim-1)+0.5*i);
				dy=j;

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

				for(k=0;k<2*num_mirrors;k+=2) {
					twodx_m=twodx+current_mirrors[k];
					dy_m=dy+current_mirrors[k+1]/2;
	
					new_sep=(twodx_m*twodx_m)+(3*dy_m*dy_m);
					sep = ( sep <= new_sep ) ? sep : new_sep;
				}

				sep_lookup[i*(6*dim-5)+j]=sqrt(sep);
			}
		}

		printf("done!\n");
	}

	dy=(int)floor(coord2/dim) - (int)floor(coord1/dim);
	twodx=2 * ( (coord2&(dim-1)) - (coord1&(dim-1)) ) + dy;

	if(dy<0) { dy*=-1; twodx*=-1; }

	return SQRT1OV8*sep_lookup[(twodx+3*(dim-1))*(6*dim-5)+dy];

}

void fft_run ( int mode, int time_index )
{
	int i,j;
	static int last_run;
	static fftw_plan fft_occs, ifft_occs;

	if( mode&1 ) {	/* Plan mode */
		for(i=0;i<n_sites;i++) temp_occs[i]=(double)occupations[i]; /* Set up temp. occupations matrix as */
									    /* matrix is modified during planning */

		printf("Planning forward FFT...");
		fft_occs=fftw_plan_dft_r2c_2d(dim,dim,temp_occs,fftout,FFTW_MEASURE); /* Plan normal transform */
		printf("done!\n");

		for(i=0;i<n_sites;i++) temp_occs[i]=(double)occupations[i];	/* Repopulate occs. matrix */
		fftw_execute(fft_occs);				/* Run normal transform */

		printf("Planning inverse FFT...");	/* Plan the inverse transform */
		ifft_occs=fftw_plan_dft_c2r_2d(dim,dim,fftout,temp_occs,FFTW_MEASURE);
		printf("done!\n");
	}

	if( mode&2 ) {	/* Inverse mode */
		/*printf("Running inverse FFT...\n");*/

		if( last_run&2 ) {	/* Do not run inverse twice in a row - result would be meaningless */
			printf("Inverse already run!\n");
		}
		else if( mode&1 ) {/* If in plan mode, inverse also not valid - forward has not been carried out yet */
			printf("Run forward transform first.\n");
		}
		else {
			fftw_execute(ifft_occs);	/* Run inverse FFT */
		}

		last_run=mode;
	}
	else {	/* Normal mode - writes output to appropriate time step position in storage matrix */
		/*printf("Running forward FFT...\n");*/

		for(i=0;i<n_sites;i++) {
			temp_occs[i]=(double)occupations[i];	/* Populate temp. occs. array */
		}

		fftw_execute(fft_occs);		/* Run */

		/* Populate full FFT output array for returning */
		for(i=dim-1;i>=0;i--) {
			for(j=dim/2;j>=0;j--) {
				fftout_full[i*dim+j+(dim/2-1)][0]=fftout[i*(dim/2+1)+j][0];
				fftout_full[i*dim+j+(dim/2-1)][1]=fftout[i*(dim/2+1)+j][1];
			}
			for(j=1;j<dim/2;j++) {
				fftout_full[i*dim+(dim/2-1)-j][0]=fftout_full[i*dim+(dim/2-1)+j][0];
				fftout_full[i*dim+(dim/2-1)-j][1]=-1*fftout_full[i*dim+(dim/2-1)+j][1];
			}
		}

		if( time_index>=0 )
			for(i=0;i<n_sites;i++) {
				fft_steps[time_index*n_sites+i][0]=fftout_full[i][0];
				fft_steps[time_index*n_sites+i][1]=fftout_full[i][1];
			}

		last_run=mode;
	}
}
