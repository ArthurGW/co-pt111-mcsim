fftw_complex * fft_mult ( int index1, int index2 )
{
	int i,j,k,offset1,offset2,mult_row;
	fftw_complex *step1, *step2, *result;
	double step1_re,step1_im,step2_re,step2_im;
	
	step1=(fftw_complex *) calloc(n_sites,2*sizeof(double));
	step2=(fftw_complex *) calloc(n_sites,2*sizeof(double));
	result=(fftw_complex *) calloc(n_sites+1,2*sizeof(double));
	
	offset1=index1*n_sites;
	offset2=index2*n_sites;
	
	for(i=0;i<n_sites;i+=dim)	/* Populate temporary matrices to multiply */
		for(j=0;j<dim;j++) {
			step1[i+j][0]=fft_steps[offset1+i+j][0];
			step1[i+j][1]=fft_steps[offset1+i+j][1];
			
			step2[i+j][0]=fft_steps[offset2+i+j][0];
			step2[i+j][1]=fft_steps[offset2+i+j][1];
		}
	
	for(i=0;i<n_sites;i+=dim)	/* Perform multiplication */
		for(j=0;j<dim;j++) {
			result[i+j][0]=0.;
			result[i+j][1]=0.;
			
			for(k=0;k<dim;k++) {
				step1_re=step1[i+k][0];
				step1_im=step1[i+k][1];
				step2_re=step2[(mult_row=k*dim)+j][0];
				step2_im=step2[mult_row+j][1];
				
				result[i+j][0]+=step1_re*step2_re-step1_im*step2_im;
				result[i+j][1]+=step1_im*step2_re+step1_re*step2_im;
			}
		}
	
	result[n_sites][0]=(double)index1;
	result[n_sites][1]=(double)index2;
	
	return result;
}