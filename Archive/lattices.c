#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <time.h>
#include <math.h>

int main ( int argc, char *argv[] )
{
	fftw_plan fwd;
	int i,j,k,now,move,atom=0,adatoms[363][2],x,y,xm,ym,orig,nns,mnns,nn[6],moves=0,rejections=0,repeats,n,xn[6],yn[6],rep;
	fftw_complex *fftout,*steps,*ISF, temp_ISF;
	double *occupations, step1_re, step2_re, step1_im, step2_im;
	FILE *out;
	long double lookup[3];
	
	occupations=(double *) fftw_malloc(33*33*sizeof(double));
	fftout=(fftw_complex *) fftw_malloc(33*33*2*sizeof(double));
	steps=(fftw_complex *) fftw_malloc(11*33*17*2*sizeof(double));
	ISF=(fftw_complex *) fftw_malloc(11*17*2*sizeof(double));
	
	for(i=0;i<11*17;i++) {
		ISF[i][0]=0.;
		ISF[i][1]=0.;
	}
	repeats=atoi(argv[1]);
	rep=atoi(argv[2]);
	printf("Repeats: %d\t Repulsion: %d\n",repeats,rep);

	
	now=(int)time(NULL);
	
	srandom(now);

	for(i=0;i<33;i++) {
		for(j=0;j<33;j++) {
			occupations[33*i+j]=floor(((double)random()/RAND_MAX)+.5);
		}
	}

	for(i=0;i<3;i++) {
		lookup[i]=(long double)exp(-0.0341308965*(i+1)*rep);
	}
	
	fwd=fftw_plan_dft_r2c_2d(33,33,occupations,fftout,FFTW_MEASURE); /* Plan forward transform */
	
	for(i=0;i<33;i++)
		for(j=0;j<33;j++)
			occupations[33*i+j]=0.;
			
	for(i=0;i<33;i++)
		for(j=(3-i%3)%3;j<33;j+=3) {
			occupations[33*i+j]=1.;
			adatoms[atom][0]=i;
			adatoms[atom++][1]=j;
		}
	printf("Atoms: %d\n",atom);
			
	fftw_execute(fwd);

	for(k=0;k<repeats;k++) {
		if(!(k&131071))printf("Repeat: %d\n",k);
	
		for(i=0;i<33*17;i++) {
			steps[i][0]=fftout[i][0];
			steps[i][1]=fftout[i][1];
		}

		for(j=1;j<11;j++){

			for(i=0;i<atom;i++) {

				/* Select random adatom */
				orig=(int)floor(atom*random()/RAND_MAX);
				y=adatoms[orig][0];
				x=adatoms[orig][1];
	
				xn[0]=(x-1+33)%33;
				xn[1]=x;
				xn[2]=(x+1)%33;
				xn[3]=(x+1)%33;
				xn[4]=x;
				xn[5]=(x-1+33)%33;
				yn[0]=(y-1+33)%33;
				yn[1]=(y-1+33)%33;
				yn[2]=y;
				yn[3]=(y+1)%33;
				yn[4]=(y+1)%33;
				yn[5]=y;
	
				move=(int)floor(6*random()/RAND_MAX);	/* Select random move */
		
				if( occupations[ yn[move]*33+xn[move] ] ) { rejections++; continue; }	/* If no move selected, go to next move */
				else { xm=xn[move]; ym=yn[move]; }

				nns=0; mnns=0;	/* Reset numbers of neighbours */

				/* Calculate current number of occupied neighbours */
				for(n=0;n<6;n++) nns+=occupations[ yn[n]*33+xn[n] ];

				xn[0]=(xm-1+33)%33;
				xn[1]=xm;
				xn[2]=(xm+1)%33;
				xn[3]=(xm+1)%33;
				xn[4]=xm;
				xn[5]=(xm-1+33)%33;
				yn[0]=(ym-1+33)%33;
				yn[1]=(ym-1+33)%33;
				yn[2]=ym;
				yn[3]=(ym+1)%33;
				yn[4]=(ym+1)%33;
				yn[5]=ym;
	

				/* Calculate potential new number of occupied neighbours */
				for(n=0;n<6;n++) mnns+=occupations[ yn[n]*33+xn[n] ];
				mnns--;	/* Reduce by one due to vacated space */
	

				/* Simple test interaction model */
				if( mnns<=nns || ( ((double)random()/RAND_MAX) <= lookup[mnns-nns-1] ) ) {
					occupations[33*y+x]=0;	/* Move adatom */
					occupations[33*ym+xm]=1;
					adatoms[orig][0]=ym;
					adatoms[orig][1]=xm;
					moves++;
				}
				else rejections++;
			}
			
			fftw_execute(fwd);

			for(i=0;i<33*17;i++) {
				steps[j*33*17+i][0]=fftout[i][0];
				steps[j*33*17+i][1]=fftout[i][1];
			}
		}

	for(j=0;j<11;j++) {

		for(i=0;i<17;i++) {	/* Perform multiplication */
			step1_re=(double)steps[i][0];	/* Take elements from lower timestep */
			step1_im=(double)steps[i][1];

			step2_re=(double)steps[j*17+i][0];	/* Take elements from higher timestep */
			step2_im=(double)(-1*steps[j*17+i][1]);	/* Take complex conjugate */
			
			ISF[j*17+i][0]+=((step1_re*step2_re-step1_im*step2_im)/(atom*repeats));
			ISF[j*17+i][1]+=((step1_im*step2_re+step1_re*step2_im)/(atom*repeats));
		}
	}
	}
	atom=0;
	for(i=0;i<33;i++) {
		for(j=0;j<33;j++) {
			if(occupations[33*i+j])atom++;
		}
	}
	printf("Atoms: %d\tMoves: %d\tRejections: %d\t Prob.: %.9f\n",atom,moves,rejections,(double)moves/(moves+rejections));
	
	for(i=0;i<33;i++) {
		for(j=0;j<33;j++) {
			printf("%.0f ",occupations[33*i+j]);
		}
		printf("\n");
	}

	out=fopen("/u/u/ajgw20/MC/sqrt3xsqrt3.csv2","w");
	for(j=0;j<11;j++) {
		fprintf(out,"%d\t",j);
		for(i=0;i<33;i+=2) fprintf(out,"%f\t%f\t\t",ISF[i*17+i/2][0],ISF[i*17+i/2][1]);
		fprintf(out,"\n");
	}
	fclose(out);
	
	fftw_free(fftout); fftw_free(steps); fftw_free(occupations); fftw_free(ISF);
	
	return 0;
}

