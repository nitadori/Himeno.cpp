#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

static const float omega=0.8;

struct Matrix {
	float* m;
	int mnums;
	int mrows;
	int mcols;
	int mdeps;

	float &operator()(int n, int r, int c, int d){
		return m[(n) * mrows * mcols * mdeps + (r) * mcols * mdeps + (c) * mdeps + (d)];
	}
	const float &operator()(int n, int r, int c, int d) const {
		return m[(n) * mrows * mcols * mdeps + (r) * mcols * mdeps + (c) * mdeps + (d)];
	}

	int newMat(int mnums,int mrows, int mcols, int mdeps)
	{
		this->mnums= mnums;
		this->mrows= mrows;
		this->mcols= mcols;
		this->mdeps= mdeps;
		this->m= NULL;
		this->m= (float*)malloc(mnums * mrows * mcols * mdeps * sizeof(float));

		return(this->m != NULL) ? 1:0;
	}

	void mat_set_init(){
		for(int i=0; i<mrows; i++)
			for(int j=0; j<mcols; j++)
				for(int k=0; k<mdeps; k++)
					(*this)(0,i,j,k) = (float)(i*i)
						/(float)((this->mrows - 1)*(this->mrows - 1));
	}

	void mat_set(int l, float val){
		for(int i=0; i<mrows; i++)
			for(int j=0; j<mcols; j++)
				for(int k=0; k<mdeps; k++)
					(*this)(l,i,j,k)=  val;
	}

	void clearMat()
	{
		if(this->m) free(this->m);
		this->m= NULL;
		this->mnums= 0;
		this->mcols= 0;
		this->mrows= 0;
		this->mdeps= 0;
	}
};

static double second()
{

	struct timeval tm;
	double t ;

	static int base_sec = 0,base_usec = 0;

	gettimeofday(&tm, NULL);

	if(base_sec == 0 && base_usec == 0)
	{
		base_sec = tm.tv_sec;
		base_usec = tm.tv_usec;
		t = 0.0;
	} else {
		t = (double) (tm.tv_sec-base_sec) + 
			((double) (tm.tv_usec-base_usec))/1.0e6 ;
	}

	return t ;
}

static double fflop(int mx,int my, int mz)
{
  return((double)(mz-2)*(double)(my-2)*(double)(mx-2)*34.0);
}

static double mflops(int nn,double cpu,double flop)
{
  return(flop/cpu*1.e-6*(double)nn);
}

void set_param(int is[],char *size)
{
	if(!strcmp(size,"XS") || !strcmp(size,"xs")){
		is[0]= 32;
		is[1]= 32;
		is[2]= 64;
		return;
	}
	if(!strcmp(size,"S") || !strcmp(size,"s")){
		is[0]= 64;
		is[1]= 64;
		is[2]= 128;
		return;
	}
	if(!strcmp(size,"M") || !strcmp(size,"m")){
		is[0]= 128;
		is[1]= 128;
		is[2]= 256;
		return;
	}
	if(!strcmp(size,"L") || !strcmp(size,"l")){
		is[0]= 256;
		is[1]= 256;
		is[2]= 512;
		return;
	}
	if(!strcmp(size,"XL") || !strcmp(size,"xl")){
		is[0]= 512;
		is[1]= 512;
		is[2]= 1024;
		return;
	} else {
		printf("Invalid input character !!\n");
		exit(6);
	}
}

float jacobi(int nn, Matrix a,Matrix b,Matrix c,
       Matrix p,Matrix bnd,Matrix wrk1,Matrix wrk2)
{
	int    imax,jmax,kmax;
	float  gosa,s0,ss;

	imax= p.mrows-1;
	jmax= p.mcols-1;
	kmax= p.mdeps-1;

	gosa = 0.0;
	for(int n=0 ; n<nn ; n++){
		gosa = 0.0;

		for(int i=1 ; i<imax; i++)
			for(int j=1 ; j<jmax ; j++)
				for(int k=1 ; k<kmax ; k++){
					s0= a(0,i,j,k)*p(0,i+1,j,  k)
						+ a(1,i,j,k)*p(0,i,  j+1,k)
						+ a(2,i,j,k)*p(0,i,  j,  k+1)
						+ b(0,i,j,k)
						*( p(0,i+1,j+1,k) - p(0,i+1,j-1,k)
								- p(0,i-1,j+1,k) + p(0,i-1,j-1,k) )
						+ b(1,i,j,k)
						*( p(0,i,j+1,k+1) - p(0,i,j-1,k+1)
								- p(0,i,j+1,k-1) + p(0,i,j-1,k-1) )
						+ b(2,i,j,k)
						*( p(0,i+1,j,k+1) - p(0,i-1,j,k+1)
								- p(0,i+1,j,k-1) + p(0,i-1,j,k-1) )
						+ c(0,i,j,k) * p(0,i-1,j,  k)
						+ c(1,i,j,k) * p(0,i,  j-1,k)
						+ c(2,i,j,k) * p(0,i,  j,  k-1)
						+ wrk1(0,i,j,k);

					ss= (s0*a(3,i,j,k) - p(0,i,j,k))*bnd(0,i,j,k);

					gosa+= ss*ss;
					wrk2(0,i,j,k)= p(0,i,j,k) + omega*ss;
				}

		for(int i=1 ; i<imax ; i++)
			for(int j=1 ; j<jmax ; j++)
				for(int k=1 ; k<kmax ; k++)
					p(0,i,j,k)= wrk2(0,i,j,k);

	} /* end n loop */

	return (gosa);
}

int main(int argc, char *argv[])
{
	int    nn;
	int    imax,jmax,kmax,mimax,mjmax,mkmax,msize[3];
	float  gosa,target;
	double cpu0,cpu1,cpu,flop;
	char   size[10];

	Matrix a,b,c,p,bnd,wrk1,wrk2; // これはmainに置く

	if(argc == 2){
		strcpy(size, argv[1]);
	} else {
		printf("For example: \n");
		printf(" Grid-size= XS (32x32x64)\n");
		printf("\t    S  (64x64x128)\n");
		printf("\t    M  (128x128x256)\n");
		printf("\t    L  (256x256x512)\n");
		printf("\t    XL (512x512x1024)\n\n");
		printf("Grid-size = ");
		scanf("%s",size);
		printf("\n");
	}

	set_param(msize,size);

	mimax= msize[0];
	mjmax= msize[1];
	mkmax= msize[2];
	imax= mimax-1;
	jmax= mjmax-1;
	kmax= mkmax-1;

	target = 60.0;

	printf("mimax = %d mjmax = %d mkmax = %d\n",mimax,mjmax,mkmax);
	printf("imax = %d jmax = %d kmax =%d\n",imax,jmax,kmax);

	/*
	 *    Initializing matrixes
	 */
	p.   newMat(1,mimax,mjmax,mkmax);
	bnd. newMat(1,mimax,mjmax,mkmax);
	wrk1.newMat(1,mimax,mjmax,mkmax);
	wrk2.newMat(1,mimax,mjmax,mkmax);
	a.   newMat(4,mimax,mjmax,mkmax);
	b.   newMat(3,mimax,mjmax,mkmax);
	c.   newMat(3,mimax,mjmax,mkmax);

	p.mat_set_init();
	bnd. mat_set(0,1.0);
	wrk1.mat_set(0,0.0);
	wrk2.mat_set(0,0.0);
	a.   mat_set(0,1.0);
	a.   mat_set(1,1.0);
	a.   mat_set(2,1.0);
	a.   mat_set(3,1.0/6.0);
	b.   mat_set(0,0.0);
	b.   mat_set(1,0.0);
	b.   mat_set(2,0.0);
	c.   mat_set(0,1.0);
	c.   mat_set(1,1.0);
	c.   mat_set(2,1.0);

	/*
	 *    Start measuring
	 */
	nn= 3;
	printf(" Start rehearsal measurement process.\n");
	printf(" Measure the performance in %d times.\n\n",nn);

	cpu0= second();
	gosa= jacobi(nn,a,b,c,p,bnd,wrk1,wrk2);
	cpu1= second();
	cpu= cpu1 - cpu0;
	flop= fflop(imax,jmax,kmax);

	printf(" MFLOPS: %f time(s): %f %e\n\n",
			mflops(nn,cpu,flop),cpu,gosa);

	nn= (int)(target/(cpu/3.0));

	printf(" Now, start the actual measurement process.\n");
	printf(" The loop will be excuted in %d times\n",nn);
	printf(" This will take about one minute.\n");
	printf(" Wait for a while\n\n");

	cpu0 = second();
	gosa = jacobi(nn,a,b,c,p,bnd,wrk1,wrk2);
	cpu1 = second();
	cpu = cpu1 - cpu0;

	printf(" Loop executed for %d times\n",nn);
	printf(" Gosa : %e \n",gosa);
	printf(" MFLOPS measured : %f\tcpu : %f\n",mflops(nn,cpu,flop),cpu);
	printf(" Score based on Pentium III 600MHz using Fortran 77: %f\n",
			mflops(nn,cpu,flop)/82.84);

	/*
	 *   Matrix free
	 */ 
	p.   clearMat();
	bnd. clearMat();
	wrk1.clearMat();
	wrk2.clearMat();
	a.   clearMat();
	b.   clearMat();
	c.   clearMat();

	return (0);
}

