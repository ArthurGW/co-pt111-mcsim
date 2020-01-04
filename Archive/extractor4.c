#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define dim 32

int main ( int argc, char *argv[] )
{
	int i,j,k, rep;
	double thetas[5]={0.065,0.1,0.2,0.3,0.4},val[2],zeroval[2];
	char inname[20], outname[20], zeroname[20], *buffer, *p;
	FILE *in, *out, *zero;
	
	rep=atoi(argv[1]);

	sprintf(outname,"fractionsrep%d.csv2",rep);
	out=fopen(outname,"w");

	for(i=1;i<(dim/2);i++) {
		sprintf(zeroname,"32-0.065-%d.csv2",rep);
		if( i==1 ) printf("zeroname: %s\n",zeroname);
		buffer=(char *)malloc(100*sizeof(char));
		k=0;
		zero=fopen(zeroname,"r");
		while(k++<i+1) fgets(buffer, 99, zero);
		fclose(zero);
		p=buffer;
		zeroval[0]=strtod(p,NULL);
		fprintf(out,"%.5f\t1.00000",zeroval[0]);
		while (*p!='\t')
			p++;
		p++;
		zeroval[0]=strtod(p,NULL);
		while (*p!='\t')
			p++;
		p++;
		zeroval[1]=strtod(p,NULL)/zeroval[0];
		free(buffer);
		fprintf(out,"\t%.5f",sqrt(2)*zeroval[1]);
		
		for(j=1;j<5;j++){
			sprintf(inname,"32-%.3f-%d.csv2",thetas[j],rep);
			if( i==1 ) printf("inname: %s\n",inname);
			buffer=(char *)malloc(100*sizeof(char));
			k=0;
			in=fopen(inname,"r");
			while(k++<i+1) fgets(buffer, 99, in);
			fclose(in);
			p=buffer;
			while (*p!='\t')
				p++;
			p++;
			val[0]=strtod(p,NULL);
			while (*p!='\t')
				p++;
			p++;
			val[1]=strtod(p,NULL)/val[0];
			free(buffer);
			
			fprintf(out,"\t%.5f\t%.5f",val[0]/zeroval[0],(val[0]/zeroval[0])*sqrt(val[1]*val[1]+zeroval[1]*zeroval[1]));
		}
		fprintf(out,"\n");	
	}
	fclose(out);
	
	return 0;
}
