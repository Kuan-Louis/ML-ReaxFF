#include <iostream>
#include <fstream>
#include<math.h>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
# include <time.h>
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
    const int k=10;                // define k: number of factors  <= change here
	const int n=1000;                // define n: number of runs    <= change here
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

