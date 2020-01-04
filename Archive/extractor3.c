#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define dim 32

int main ( int argc, char *argv[] )
{
	int i,j,k, reps[17]={0,1,2,5,10,20,50,100,200,300,400,500,600,700,800,900,1000};
	double thetas[5]={0.065,0.1,0.2,0.3,0.4},meansquare[2],errorval,val,zeroval[dim/2][2];
	char inname[20], outname[20], zeroname[20], *buffer, *p;
	FILE *in, *out, *zero;
	
	sprintf(outname,"seriesthetas.csv2");
	out=fopen(outname,"w");

	for(i=0;i<17;i++) {
		fprintf(out,"%d",reps[i]);
		sprintf(zeroname,"32-0.065-%d.csv2",reps[i]);
		printf("zeroname: %s\n",zeroname);
		zero=fopen(zeroname,"r");
		for(k=0;k<(dim/2);k++) {
				buffer=(char *)malloc(100*sizeof(char));
				fgets(buffer, 99, zero);
				p=buffer;
				while (*p!='\t')
					p++;
				p++;
				zeroval[k][0]=strtod(p,NULL);
				while (*p!='\t')
					p++;
				p++;
				zeroval[k][1]=strtod(p,NULL);
				zeroval[k][1]=zeroval[k][1]*zeroval[k][1];
				free(buffer);
		}
		fclose(zero);
		
		for(j=0;j<5;j++){
			meansquare[0]=0.;
			meansquare[1]=0.;
			sprintf(inname,"32-%.3f-%d.csv2",thetas[j],reps[i]);
			printf("inname: %s\n",inname);
			in=fopen(inname,"r");
			for(k=0;k<(dim/2);k++) {
				buffer=(char *)malloc(100*sizeof(char));
				fgets(buffer, 99, in);
				p=buffer;
				while (*p!='\t')
					p++;
				p++;
				val=strtod(p,NULL)-zeroval[k][0];
				while (*p!='\t')
					p++;
				p++;
				errorval=strtod(p,NULL);
				errorval=sqrt(errorval*errorval+zeroval[k][1]);
				meansquare[0]+=val*val;
				meansquare[1]+=val*val*errorval*errorval;
				free(buffer);
			}
			fclose(in);
			if( meansquare[0]==0. ) fprintf(out,"\t0.00000\t0.00000");
			else fprintf(out,"\t%.5f\t%.5f",sqrt(2*meansquare[0]/dim),sqrt(2*meansquare[1]/(dim*meansquare[0])));
		}
		fprintf(out,"\n");	
	}
	fclose(out);
	
	return 0;
}
