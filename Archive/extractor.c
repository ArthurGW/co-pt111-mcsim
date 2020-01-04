#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define ISF_VECT_SCL 1.960985419	/* Factor used to calculate length of <11-2> inverse azimuth */
#define TIMEERROR 0.44995

int main ( int argc, char *argv[] )
{
	int i,j,linecount=0,dim;
	double vect_len,result[4],error_fract;
	char inname[20], outname[20], *buffer, *p;
	FILE *in, *out;
	
	sprintf(inname,"%s.prm",argv[1]);
	sprintf(outname,"%s.csv",argv[1]);
	dim=atoi(argv[2]);
	
	vect_len=ISF_VECT_SCL*(1.-(2./dim));
	
	in=fopen(inname,"r");
	out=fopen(outname,"w");
	
	for(i=0;i<dim/2;i++) {
		fprintf(out, "%d\t%.5f", i, i*vect_len/((dim/2.)-1.) );
		for(j=0;j<4;j++){
			if(j==2) {
				error_fract = (result[0]!=0.) ? (result[1]/result[0])+TIMEERROR : 0.;
				fprintf(out, "\t%.7f", result[0]*error_fract);
			}
			buffer=(char *)malloc(50*sizeof(char));
			fgets(buffer, 49, in);
			linecount++;
			p=buffer;
			while (*p!='=')
				p++;
			p++;
			result[j]=strtod(p,NULL);
			fprintf(out, "\t%.7f", result[j]);
			free(buffer);
		}
		fprintf(out, "\n");
	}
	
	fclose(in);
	fclose(out);
	
	return 0;
}
