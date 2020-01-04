double separation1 ( int coord1, int coord2 )
{
	unsigned int i, j, /*jmin=0, jmax,*/ num_mirrors, sep, new_sep; 
	signed int twodx, twodx_m, dy, dy_m;

	signed int mirrors[8];
	signed int current_mirrors[6]={0,0,0,0,0,0};

	static double lookup[6*dim-5][dim];

	if( !coord1 & !coord2 ) {
		printf("Creating separation lookup table...");

		mirrors[0]=dim;
		mirrors[1]=2*dim;
		mirrors[2]=3*dim-1;
		mirrors[3]=2*dim-2;
		mirrors[4]=2*dim-1;
		mirrors[5]=-2;
		mirrors[6]=dim-1;
		mirrors[7]=-2*dim-2;

		for(i=0;i<6*dim-5;i++) {
			/*if(i<dim) jmax=(int)floor(0.5*i);*/

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

				for(i=0;i<2*num_mirrors;i+=2) {
					twodx_m=twodx+current_mirrors[i];
					dy_m=dy+current_mirrors[i+1]/2;
	
					new_sep=(twodx_m*twodx_m)+(3*dy_m*dy_m);
					sep = ( sep <= new_sep ) ? sep : new_sep;
				}

				lookup[i][j]=sqrt(sep);
			}
		}

		printf("done!\n");
	}

	dy=(int)floor(coord2/dim) - (int)floor(coord1/dim);
	twodx=2 * ( (coord2&(dim-1)) - (coord1&(dim-1)) ) + dy;

	if(dy<0) { dy*=-1; dx*=-1; }

	return SQRT1OV8*lookup[(0.5*twodx+1.5*(dim-1))*2][dy];

}
