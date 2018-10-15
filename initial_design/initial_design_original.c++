#include <iostream>
#include <fstream>
#include<math.h>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define	 int8	      unsigned int
#define  Min(a,b)     (a<b)?a:b
#define	 MAX	      (int8) (pow(2,31)-1)
#define  autoseed(x)  (int8)(rand()*x+x)

int rc(int n ,double seed) ;
double runif(double seed);
void djinformation(int **A, int n, int k, int* );
int rc2(int n,int del,double seed);
double mcorr(int **A, int n,int k);
double corrd(int **A, int n,int a, int b);
int **LHD(int n,int k,double seed);
int **SLHD(int n,int k,double seed);
double score(int **A, int n, int k, int p) ;
double UB(int n,int k,int p);
double LB(int n,int k,int p);
int column(int **A, int n,int k);
int element(int **A, int n, int k,int p);
int main()
{
	srand(time(0));
   std:: ofstream out;  //¼g¨ì·sÀÉ®×¤º
    out.open("parameters"); // change the name of output file
	double seed=rand();
    const int k=5;                // define k: number of factors  <= change here
	const int n=10;                // define n: number of runs    <= change here
	int nsearch=5;             // number of different initials
	int p=100;
	double tfac=0.95;
	const double eps=0.0000001;
	int nterp=150;//10*19*4;   // parameter for Simulated Annealing
	double t0;
	double bestyet=100000;
	double crit;
	double cor;
	double newcrit;
	double newcritbest;
	double critbest;
    double corbest;
	double nbestyet=0;
	int itotal;
	double temp;
	int **xbest;
    int **xtry;
	int **x;
	int *Mminf;
    Mminf	=new int[2];
 	int i;


// Initialize xtry, xbest, and x;
	xtry	=new int*[k];
    xbest	=new int*[k];
	x		=new int*[k];
	for(i=0;i<k;i++)
	{
		xbest[i] =new int[n];
		xtry[i]  =new int[n];
		x[i]     =new int[n];
	}

/////calculate starting temperature.
	double dbar=(n+1)*k*pow(3,(-1));
	double c=0.5;
	double p0=0.99;
	double pf=0.01;
	double t1=(pow(dbar,(-p))*n*(n-1))*pow(4*c*(p-1),(-1));
	double t2=pow((1-c),(1-p));
	double t3=pow((1+c),(1-p));
	double phi0=pow((t1*(t2-t3)),(pow(p,(-1))));
	double delt0=pow((pow(phi0,p)-pow((1-c)*dbar,(-p))+pow((1-c)*dbar-1,(-p))-pow((1+c)*dbar,(-p))+pow((1+c)*dbar+1,(-p))),pow(p,(-1)));
	t0=-delt0*pow((log(p0)),(-1));
	//ut<<"start"<<"\n";
	printf("start");

	for(int loop=0;loop<5;loop++)
	{
		nterp=nterp+100;
		out<<"nterp="<<nterp<<"\n";
		out<<"t0="<<t0<<"\n";
///////////////////////////////////
		x=LHD(n,k,seed);
////repeated search loop //////////
		for(int isearch=1;isearch<nsearch;isearch++)
		{
//		cout<<"isearch="<<isearch<<"\n";
			seed=rand();
///// generate LHD //////////////
//       x=LHD(n,k,seed);
//////initialized the best design ////////////////
			for(int n2=0;n2<k;n2++)
			{
				for(int n1=0;n1<n;n1++)
				{
				 *(*(xbest+n2)+n1)=*(*(x+n2)+n1);
				 *(*(xtry+n2)+n1)=*(*(x+n2)+n1);
				}
			}
			crit=score(x,n,k,p);
			cor=mcorr(x,n,k);
			critbest=crit;
			corbest=cor;
			newcritbest=pow(corbest,2)+(critbest-LB(n,k,p))*pow((UB(n,k,p)-LB(n,k,p)),(-1));
		//	cout<<"tempsc="<<newcritbest<<"\n";
///////////initialize tempertures and counts///////////////////
			temp=t0;
			itotal=0;
			int ichange=1;
////////// variable temperature loop ////////////////////////
     		while(ichange==1)
			{
				ichange=0;
	//// constant temperature loop /////////////////////
				int ipert=1;
				while (ipert<nterp)
				{
					itotal=itotal+1;
       //////// switch to be tried is elements
		 /////// change two component in a column //////////
					 int ind;
					 int tran1;
					 int tran2;
					 seed=rand();
					 //ind=rc((k-1),seed);
					 ind=column(x,n,k); ind=ind-1;
					 seed=rand();
					 //tran1=rc(n,seed);
					 tran1=element(x,n,k,p);
					 seed=rand();
					 tran2=rc2(n,tran1,seed);
       /////////perturb x to xtry////////////////////
 					 *(*(xtry+ind+1)+tran2)=*(*(x+ind+1)+tran1);
				     *(*(xtry+ind+1)+tran1)=*(*(x+ind+1)+tran2);
       //////////////////////////////////////////////
					 double  crittry=score(xtry,n,k,p);
					 double  cortry=mcorr(xtry,n,k);
					 double  newcrittry=pow(cortry,2)+(crittry-LB(n,k,p))*pow((UB(n,k,p)-LB(n,k,p)),(-1));
////// is try better than best? //////////////////////////////
					if(newcrittry<newcritbest)
					{
	//////////yes: replace x, xbest by xtry ; set iterp=1;ichange=1////////////////////////
						ichange=1;
						for(int nn2=0;nn2<k;nn2++)
						{
						   for(int nn1=0;nn1<n;nn1++)
						   {
							  *(*(xbest+nn2)+nn1)=*(*(xtry+nn2)+nn1);
						   }
						}
						 *(*(x+ind+1)+tran1)=*(*(xtry+ind+1)+tran1);
						 *(*(x+ind+1)+tran2)=*(*(xtry+ind+1)+tran2);
						critbest=crittry;
						corbest=cortry;
						newcritbest=newcrittry;
						ipert=1;
						crit=crittry;
						cor=cortry;
						newcrit=newcrittry;

					}
					else
					{
//////////No:, increase ipert by 1. is xtry as good or better than x?
						ipert=ipert+1;
						//if(crittry<crit)
						if(newcrittry<newcrit)
						{
			  ////// xtry is better than x; replace x by xtry ///////////
 						 *(*(x+ind+1)+tran1)=*(*(xtry+ind+1)+tran1);
						 *(*(x+ind+1)+tran2)=*(*(xtry+ind+1)+tran2);
						 ichange=1;
						 crit=crittry;
						 cor=cortry;
						 newcrit=newcrittry;
		  //////////////////////////////////////////////////////////
					}
					else
					{
		 ///////// xtry is worst than x////////////////////////////
					    double delta1=crittry-crit;
						double delta2=cortry-cor;
						double range=UB(n,k,p)-LB(n,k,p);
					//    double prob=exp(-delta1*pow(temp,(-1)));
			        	double prob=exp(-delta1*pow(temp,(-1))-delta2*range*pow(temp,(-1)));
 					    seed=seed+isearch+ipert;
					    double q=runif(seed);
					     if(prob>=q)
						 {///// replce x by xtry by prob///////////
 				           	 *(*(x+ind+1)+tran1)=*(*(xtry+ind+1)+tran1);
				             *(*(x+ind+1)+tran2)=*(*(xtry+ind+1)+tran2);
							 ichange=1;
							 crit=crittry;
							 cor=cortry;
							 newcrit=newcrittry;
						 }///////////////////////////////////
                         else
						 {///// reset x try to x for the next pertubation
 				           	 *(*(xtry+ind+1)+tran1)=*(*(x+ind+1)+tran1);
				             *(*(xtry+ind+1)+tran2)=*(*(x+ind+1)+tran2);
						 }//////////////////////////////////////////
		 //////////////////////////////////////////////////////////
					}
				}
			}
	//// end of constant temperature loop ////////////
			temp=temp*tfac;
		}
///////// End of variable temperature loop///////////////////
////// result of this search//////////////////////////////
        if(newcritbest<(bestyet+eps))
		{
			nbestyet=nbestyet+1;
		}
    	if(newcritbest<(bestyet-eps))
		{
			bestyet=newcritbest;
			nbestyet=1;
		}

		for(int ii=0;ii<n;ii++)
		{
			for(int jj=0;jj<k;jj++)
			{
				out<<*(*(xbest+jj)+ii)<<"  ";
			}
			out<<"\n";
		}
	double combine=pow(mcorr(xbest,n,k),2)+(score(xbest,n,k,p)-LB(n,k,p))*pow((UB(n,k,p)-LB(n,k,p)),(-1));

	djinformation(xbest,n, k, Mminf);
	out<<"d="<<*(Mminf+0)<<",j="<<*(Mminf+1)<<"\n";
	out<<"score="<<newcritbest<<"\n";
	out<<"isearch="<<isearch<<"\n";
	out<<"Number of design examed="<<itotal<<"\n";
	out<<"phi_p="<<score(xbest,n,k,p)<<"\n";
	out<<"correlation="<<mcorr(xbest,n,k)<<"\n";

//    out<<score(xbest,n,k,p)<<"  "<<mcorr(xbest,n,k)<<"\n";
	/////end of search loop////////////////////////////
   }
}

 // delete xtry, xbest, and x;
	for(i=0; i<k; i++)
	{
		delete [] xbest[i];
		delete [] xtry[i];
		delete [] x[i];
	}

	delete []xbest;
	delete []xtry;
	delete []x;
    delete []Mminf;
}

int **LHD(int n, int k, double seed)
{
	int te;
	int **LHD;
	int *r;
	r=new int [n];
    LHD=new int*[k];
	for(int iii=0;iii<k;iii++)
	{
		LHD[iii]=new int[n];
	}
     for(int j2=0;j2<n;j2++)  // first dimention is 1,2,3,.....
	{
		*(*(LHD+0)+j2)=(j2+1);
//        cout<<*(*(LHD+0)+j2)<<"  ";
	}
//	cout<<"\n";
	for(int j1=0;j1<(k-1);j1++)
	{
	   for(int cc=0;cc<n;cc++)
	   {
		*(r+cc)=cc+1;
	   }
	   for(int c=0;c<n;c++)
	   {
			seed=seed+j1*c;
			seed=seed+10;
			te=rc(n-c,seed);
			*(*(LHD+j1+1)+c)=*(r+te);
//		cout<<*(*(LHD+j1+1)+c)<<"  "; // print sample
		    for(int c1=0;c1<(n-c);c1++)
			{
			    if(c1>=te)
				{
				*(r+c1)=*(r+c1+1);
				}
			}
	   }//cout<<"\n";
	}
   return(LHD);
   for(int iiii=0;iiii<k;iiii++)
   {
	delete [] LHD[iiii];
   }
   delete [] LHD;
   delete [] r;
}


int rc2(int n,int del,double seed)
{
   int rctwo;

   rctwo= rc( n-1, seed);
   if (rctwo >= del)  rctwo++;

   return(rctwo);
}

int rc(int n,double seed ) // choose randomly from 0 to (n-1)
{
   int r;
   double u; // assume I can find number from U(0,1)
   u = (double) (rand()*autoseed(seed)%MAX)/MAX;
   r=(int)(n*u);
   return(r);
}

double runif(double seed)
{
   double runif;
   runif= (double) (rand()*autoseed(seed)%MAX)/MAX;
   return(runif);
}

double score(int **A, int n, int k,int p) // try to input a LHD metrix A[][]
{
//	int A[3][6];

	//int n=6;  // number of points
	//int k=3;  // number of dimention
	int i;
    int dim=n*(n-1)/2;
	double score=0;
	int count=0;
	int c=0;

	int *d;
	int *j;
	d=new int[dim];
	j=new int[dim];

	for(i=0;i<dim;i++)
	{
		d[i] = 0;
		j[i] = 1;
	}
 	count= c = 0;
   	for(int k1=0;k1<(n-1);k1++)
	{
		for(int k2=(k1+1);k2<n;k2++)
		{
 			d[count++] = 0;
		//	cout<< "count="<<count<<"\n ";
			for(int k3=0;k3<k;k3++)
			{
				d[count-1] += abs(A[k3][k1]-A[k3][k2]);
			}
			for(int k4=0;k4<(count-1);k4++)
			{
				if(d[count-1]== d[k4])
				{
					j[k4]++;
                    count--;
				}
			}
		}

	}
 	for(i=0; i<count; i++)
	{
		score += (j[i])* pow(d[i],(-p));
	}

	score=pow(score,(pow(p,(-1))));

	delete [] d;
	delete [] j;

	return(score);
}

void djinformation(int **A, int n, int k, int* Mminf) // try to input a LHD metrix A[][]
{
    const int dim= (int)(n*(n-1)*0.5);
	double score=0;

	int *d;
	int *j;
	d=new int[dim];
	j=new int[dim];

	for(int i=0;i<dim;i++)
	{
		d[i] = 0;
		j[i] = 1;
	}

 	int count=0;
	int c=0;

   	for(int k1=0;k1<(n-1);k1++)
	{
		for(int k2=(k1+1);k2<n;k2++)
		{
 			d[count]=0;
			count++;

			for(int k3=0;k3<k;k3++)
			{
				d[count-1] += abs(*(*(A+k3)+k1)-*(*(A+k3)+k2));
			}
			for(int k4=0;k4<(count-1);k4++)
			{
				if(d[count-1]==d[k4])
				{
					j[k4]++;
                    count--;
				}
			}
		}
	}

	int Mmd=100000;
	int Mmj=0;

	for(int k5=0;k5<count;k5++)
	{
		if(Mmd>d[k5])
		{
			Mmd= d[k5];
			Mmj= j[k5];
		}
	}

	Mminf[0] = Mmd;
	Mminf[1] = Mmj;

	delete [] d;
	delete [] j;

}


double corrd(int **A, int n,int a, int b)
{
	double cor=0;
	double cov=0;
	double var=0;
	const double p2tm1 = pow(2,(-1));

	for(int i1=0;i1<n;i1++)
	{
 		cov += (A[a][i1]-(n+1)*p2tm1) * (A[b][i1]-(n+1)*p2tm1);
		var += (A[0][i1]-(n+1)*p2tm1) * (A[0][i1]-(n+1)*p2tm1);
	}
	cor = cov * pow(var,(-1));
	return(cor);
}


double mcorr(int **A, int n,int k)
{
	double mcor=0;

 	for(int j1=1;j1<k;j1++)
	{
		for(int j2=0;j2<j1;j2++)
		{
			mcor += pow(corrd(A,n,j1,j2),2);
		}
	}

	mcor = 2 * mcor / (double)(k*(k-1));
	mcor = pow(mcor,(0.5));
	return(mcor);
}

double UB(int n, int k,int p)
{
	double ub=0;
	for(int i=1;i<(n-1);i++)
	{
		ub += (n-i)* pow((i*k),(-p));
	}
	ub = pow(ub,(pow(p,(-1))));
	return(ub);
}


double LB(int n,int k,int p)
{
	double lb=0;
	lb=3*pow((n+1)*k,(-1))*pow(n*(n-1)*pow(2,(-1)),(pow(p,(-1))));
	return(lb);
}


int column(int **A, int n,int k)
{
	int col=0;
	int i;
	double maxcorr=0;


	double *correlation;
	correlation=new double[k];

 	for(i=0;i<k;i++)
	correlation[i] = 0;

	for(int p1=0;p1<k;p1++)
	{
		for(int p2=0;p2<k;p2++)
		{
			if(p2!=p1)
			{
				correlation[p1] += corrd(A,n,p1,p2);
			}
		}

		if (correlation[p1] < 0 ) correlation[p1] *= -1;

		if(correlation[p1] >= maxcorr)
		{
			maxcorr = correlation[p1];
            col=p1;
		}
	}
	delete [] correlation;

	return(col);

}

int element(int **A, int n, int k,int p)
{
	int ele=0;
	double distan=0;
	double maxscore=0;
	double *dist;
	dist=new double[n];
	int pp1, pp2, pp3;
	int i;

	for(i=0; i<n; i++)	dist[i] = 0;

	for(pp1=0;pp1<n-1;pp1++)
	{
		for(pp2=pp1+1;pp2<n;pp2++)
		{
			distan = 0;
			for(pp3=0;pp3<k;pp3++)
			{
				distan += abs(A[pp3][pp1]-A[pp3][pp2]);
			}
			dist[pp1] += pow(distan,(-p));
			dist[pp2] += pow(distan,(-p));
		}
	}

	for(i=0; i<n; i++)
	{
		if(dist[i]>maxscore)
        {
			maxscore = dist[i];
			ele=i;
		}
	}

	delete [] dist;
	return ele;
}
