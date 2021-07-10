#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include<conio.h>

//     **** SUBROUTINE PROGRAMME ****
float g(float xx[], int nn)
{  int i;	float res;
   res=xx[1];
   for (i=2;i<=nn;i++)
      if (xx[i]>res)	res=xx[i];
   return(res);
}


int main()
{
//srand(time(NULL));
int A1[201][201],B1[201][201],B2[201][201],deg1[201],deg2[201],i,j,j1,ii,n,m,mm,j2,j3,i1,k1,p,i2,i3,tn,l1,k2,k3,k4,k5,nn,i4,it;

float x0[201],x[201],y0[201],y[201],z0[201],z[201],k[5][1201],ax0[201],ax[201],ay0[201],ay[201],az0[201],az[201],dm[201];
float RAND,sum3,prand,avger,sum1,w0,h,l,w,a,c,b,eta,error,rf,rp,np,swp,pinter,temp,t0,sigma,beta,rho,dr,eps,sum2;
float x20[201],x2[201],y20[201],y2[201],z20[201],z2[201];

FILE *fopen(),*fp1,*fp2;
fp1=fopen("L_RIS_DM_ts.out","w");
fp2=fopen("L_RIS_DM_MSF_0.out","w");
it=0;

n=200;
nn=(3*n)*(3*n+1);
k1=2; // k SHOULD BE LESS OR EQUAL TO n/4

float* v = calloc(nn+1, sizeof(float));
float* znorm = calloc(nn+1, sizeof(float));
float* gsc = calloc(nn+1, sizeof(float));
float* cum = calloc(nn+1, sizeof(float));

mm=200000;// Time Iteration
tn=150000;// Transient

sigma=10.0;	beta=8.0/3.0;	rho=28.0;
prand=0.1;// Probability Of The Random Network

dr=3.0;
eps=3.5;
eta=3.0;

			//Rewiring Of The Layer-1,3 Starts
			for(i=1;i<=n;i++) for(j=1;j<=n;j++) A1[i][j]=0;
			for(i=1;i<=n;i++)
			{
				for(j=i+1;j<=n;j++)
			  {
				  if((float)rand()/RAND_MAX<=prand)
			    {
				    A1[i][j]=1;
				    A1[j][i]=1;
			    }
			  }
			}

			for(i=1;i<=n;i++)
			{
				deg1[i]=0;
			  for(j=1;j<=n;j++)
					deg1[i]=deg1[i]+A1[i][j];
			}

			for(i=1;i<=n;i++)
			{
				m=1;
			  for(j=1;j<=n;j++)
			  {
			  	if(A1[i][j]==1)
			    {
			    	B1[i][m]=j;
				    m=m+1;
			    }
			  }
			}
			//Rewiring Of The Layer-1,3 Ends

			//Rewiring Of The Layer-2 Starts
			for(i=1;i<=n;i++) for(j=1;j<=n;j++) A1[i][j]=0;
			for(i=1;i<=n;i++)
			{
				for(j=i+1;j<=n;j++)
			  {
				  if((float)rand()/RAND_MAX<=prand)
			    {
				    A1[i][j]=1;
				    A1[j][i]=1;
			    }
			  }
			}

			for(i=1;i<=n;i++)
			{
				deg2[i]=0;
				for(j=1;j<=n;j++)
					deg2[i]=deg2[i]+A1[i][j];
			}
			for(i=1;i<=n;i++)
			{
				m=1;
				for(j=1;j<=n;j++)
			  {
				  if(A1[i][j]==1)
			    {
				    B2[i][m]=j;
			      m=m+1;
			    }
			  }
			}
//	Rewiring Of The Layer-2 Ends

		for(i2=12*it;i2<=12*it+6;i2=i2+6)//(i2=0;i2<=n;i2=i2+6)
		{
      for(i=1;i<=n;i++)
			{
				if(i>i2)
        	dm[i]=eta;
        else
        	dm[i]=0.0;
      }

			for(i=1;i<=n;i++)
			{
				x0[i]=0.01*(-1.0+2.0*(float)rand()/RAND_MAX);		y0[i]=0.01*(-1.0+2.0*(float)rand()/RAND_MAX);		z0[i]=0.01*((float)rand()/RAND_MAX);
				x20[i]=0.01*(-1.0+2.0*(float)rand()/RAND_MAX);	y20[i]=0.01*(-1.0+2.0*(float)rand()/RAND_MAX);	z20[i]=0.01*((float)rand()/RAND_MAX);
			}

h=0.01;
w0=0.0;
error=0.0;

for(i=3*n+1;i<=nn;i++)
   v[i]=0.0;
for(i=1;i<=3*n;i++)
{
   v[(3*n+1)*i]=1.0;
   cum[i]=0.0;
}

for(ii=1;ii<=mm;ii++)
{			//	printf("%d\n",ii);

				for(j=1;j<=n;j++)
				{
					sum1=0.0;
		      for(m=1;m<=deg1[j];m++)
						sum1=sum1+x0[B1[j][m]];
		      l=sigma*(y0[j]-x0[j])+eps*(sum1-deg1[j]*x0[j])+dm[j]*(x20[j]-x0[j]);
		      k[1][6*j-5]=l*h;
		      l=x0[j]*(rho-z0[j])-y0[j]+dm[j]*(y20[j]-y0[j]);
		      k[1][6*j-4]=l*h;
		      l=x0[j]*y0[j]-beta*z0[j]+dm[j]*(z20[j]-z0[j]);
		      k[1][6*j-3]=l*h;

					sum1=0.0;
		      for(m=1;m<=deg2[j];m++)
          	sum1=sum1+x20[B2[j][m]];
		      l=sigma*(y20[j]-x20[j])+eps*(sum1-deg2[j]*x20[j])+2.0*dm[j]*(x0[j]-x20[j]);
		      k[1][6*j-2]=l*h;
		      l=x20[j]*(rho+dr-z20[j])-y20[j]+2.0*dm[j]*(y0[j]-y20[j]);
		      k[1][6*j-1]=l*h;
		      l=x20[j]*y20[j]-beta*z20[j]+2.0*dm[j]*(z0[j]-z20[j]);
		      k[1][6*j]=l*h;
				}

				for(j=1;j<=n;j++)
			  {
      		x[j]=x0[j]+k[1][6*j-5]/2.0;
      		y[j]=y0[j]+k[1][6*j-4]/2.0;
      		z[j]=z0[j]+k[1][6*j-3]/2.0;
      		x2[j]=x20[j]+k[1][6*j-2]/2.0;
      		y2[j]=y20[j]+k[1][6*j-1]/2.0;
      		z2[j]=z20[j]+k[1][6*j]/2.0;
   			}//	w=w0 +h/2.0


   			for(j=1;j<=n;j++)
			  {
					sum1=0.0;
		      for(m=1;m<=deg1[j];m++)
						sum1=sum1+x[B1[j][m]];
		      l=sigma*(y[j]-x[j])+eps*(sum1-deg1[j]*x[j])+dm[j]*(x2[j]-x[j]);
		      k[2][6*j-5]=l*h;
		      l=x[j]*(rho-z[j])-y[j]+dm[j]*(y2[j]-y[j]);
		      k[2][6*j-4]=l*h;
		      l=x[j]*y[j]-beta*z[j]+dm[j]*(z2[j]-z[j]);
		      k[2][6*j-3]=l*h;

					sum1=0.0;
		      for(m=1;m<=deg2[j];m++)
          	sum1=sum1+x2[B2[j][m]];
		      l=sigma*(y2[j]-x2[j])+eps*(sum1-deg2[j]*x2[j])+2.0*dm[j]*(x[j]-x2[j]);
		      k[2][6*j-2]=l*h;
		      l=x2[j]*(rho+dr-z2[j])-y2[j]+2.0*dm[j]*(y[j]-y2[j]);
		      k[2][6*j-1]=l*h;
		      l=x2[j]*y2[j]-beta*z2[j]+2.0*dm[j]*(z[j]-z2[j]);
		      k[2][6*j]=l*h;
			  }

   			for(j=1;j<=n;j++)
			  {
      		x[j]=x0[j]+k[2][6*j-5]/2.0;
      		y[j]=y0[j]+k[2][6*j-4]/2.0;
      		z[j]=z0[j]+k[2][6*j-3]/2.0;
      		x2[j]=x20[j]+k[2][6*j-2]/2.0;
      		y2[j]=y20[j]+k[2][6*j-1]/2.0;
      		z2[j]=z20[j]+k[2][6*j]/2.0;
   			}//	w=w0+h/2.0


		   	for(j=1;j<=n;j++)
   			{
					sum1=0.0;
		      for(m=1;m<=deg1[j];m++)
						sum1=sum1+x[B1[j][m]];
		      l=sigma*(y[j]-x[j])+eps*(sum1-deg1[j]*x[j])+dm[j]*(x2[j]-x[j]);
		      k[3][6*j-5]=l*h;
		      l=x[j]*(rho-z[j])-y[j]+dm[j]*(y2[j]-y[j]);
		      k[3][6*j-4]=l*h;
		      l=x[j]*y[j]-beta*z[j]+dm[j]*(z2[j]-z[j]);
		      k[3][6*j-3]=l*h;

					sum1=0.0;
		      for(m=1;m<=deg2[j];m++)
          	sum1=sum1+x2[B2[j][m]];
		      l=sigma*(y2[j]-x2[j])+eps*(sum1-deg2[j]*x2[j])+2.0*dm[j]*(x[j]-x2[j]);
		      k[3][6*j-2]=l*h;
		      l=x2[j]*(rho+dr-z2[j])-y2[j]+2.0*dm[j]*(y[j]-y2[j]);
		      k[3][6*j-1]=l*h;
		      l=x2[j]*y2[j]-beta*z2[j]+2.0*dm[j]*(z[j]-z2[j]);
		      k[3][6*j]=l*h;
			  }

			  for(j=1;j<=n;j++)
			  {
				  x[j]=x0[j]+k[3][6*j-5];
				  y[j]=y0[j]+k[3][6*j-4];
				  z[j]=z0[j]+k[3][6*j-3];
				  x2[j]=x20[j]+k[3][6*j-2];
				  y2[j]=y20[j]+k[3][6*j-1];
				  z2[j]=z20[j]+k[3][6*j];
			  }//       w=w0+h


			 	for(j=1;j<=n;j++)
			  {
					sum1=0.0;
		      for(m=1;m<=deg1[j];m++)
						sum1=sum1+x[B1[j][m]];
		      l=sigma*(y[j]-x[j])+eps*(sum1-deg1[j]*x[j])+dm[j]*(x2[j]-x[j]);
		      k[4][6*j-5]=l*h;
		      l=x[j]*(rho-z[j])-y[j]+dm[j]*(y2[j]-y[j]);
		      k[4][6*j-4]=l*h;
		      l=x[j]*y[j]-beta*z[j]+dm[j]*(z2[j]-z[j]);
		      k[4][6*j-3]=l*h;

					sum1=0.0;
		      for(m=1;m<=deg2[j];m++)
          	sum1=sum1+x2[B2[j][m]];
		      l=sigma*(y2[j]-x2[j])+eps*(sum1-deg2[j]*x2[j])+2.0*dm[j]*(x[j]-x2[j]);
		      k[4][6*j-2]=l*h;
		      l=x2[j]*(rho+dr-z2[j])-y2[j]+2.0*dm[j]*(y[j]-y2[j]);
		      k[4][6*j-1]=l*h;
		      l=x2[j]*y2[j]-beta*z2[j]+2.0*dm[j]*(z[j]-z2[j]);
		      k[4][6*j]=l*h;
				}

				for(j=1;j<=n;j++)
				{
				  x0[j]=x0[j]+(k[1][6*j-5]+2.0*(k[2][6*j-5]+k[3][6*j-5])+k[4][6*j-5])/6.0;
				  y0[j]=y0[j]+(k[1][6*j-4]+2.0*(k[2][6*j-4]+k[3][6*j-4])+k[4][6*j-4])/6.0;
				  z0[j]=z0[j]+(k[1][6*j-3]+2.0*(k[2][6*j-3]+k[3][6*j-3])+k[4][6*j-3])/6.0;
				  x20[j]=x20[j]+(k[1][6*j-2]+2.0*(k[2][6*j-2]+k[3][6*j-2])+k[4][6*j-2])/6.0;
				  y20[j]=y20[j]+(k[1][6*j-1]+2.0*(k[2][6*j-1]+k[3][6*j-1])+k[4][6*j-1])/6.0;
				  z20[j]=z20[j]+(k[1][6*j]+2.0*(k[2][6*j]+k[3][6*j])+k[4][6*j])/6.0;
				}

/////////////////////////////////////////////////////////////////////////////////////////////////
   if(ii>tn)
   {
      for(j1=0;j1<=3*n-1;j1++)
      {
         for(j=1;j<=n;j++)
         {
            ax0[j]=v[(3*j-2)*3*n+1+j1];
            ay0[j]=v[(3*j-1)*3*n+1+j1];
            az0[j]=v[(3*j)*3*n+1+j1];
         }

         for(j=1;j<=n;j++)
         {
            sum1=0.0;
            for(m=1;m<=deg1[j];m++)
               sum1=sum1+ax0[B1[j][m]];
            l=sigma*(ay0[j]-ax0[j])+eps*(sum1-deg1[j]*ax0[j])-dm[j]*ax0[j];
            k[1][3*j-2]=l*h;
            l=ax0[j]*(rho-z0[j])-x0[j]*az0[j]-ay0[j]-dm[j]*ay0[j];
            k[1][3*j-1]=l*h;
            l=ax0[j]*y0[j]+x0[j]*ay0[j]-beta*az0[j]-dm[j]*az0[j];
            k[1][3*j]=l*h;
         }

         for(j=1;j<=n;j++)
         {
            ax[j]=ax0[j]+k[1][3*j-2]/2.0;
            ay[j]=ay0[j]+k[1][3*j-1]/2.0;
            az[j]=az0[j]+k[1][3*j]/2.0;
         }//	w=w0 +h/2.0
       
         for(j=1;j<=n;j++)
         {
            sum1=0.0;
            for(m=1;m<=deg1[j];m++)
               sum1=sum1+ax[B1[j][m]];
            l=sigma*(ay[j]-ax[j])+eps*(sum1-deg1[j]*ax[j])-dm[j]*ax[j];
            k[2][3*j-2]=l*h;
            l=ax[j]*(rho-z0[j])-x0[j]*az[j]-ay[j]-dm[j]*ay[j];
            k[2][3*j-1]=l*h;
            l=ax[j]*y0[j]+x0[j]*ay[j]-beta*az[j]-dm[j]*az[j];
            k[2][3*j]=l*h;
         }

         for(j=1;j<=n;j++)
         {
            ax[j]=ax0[j]+k[2][3*j-2]/2.0;
            ay[j]=ay0[j]+k[2][3*j-1]/2.0;
            az[j]=az0[j]+k[2][3*j]/2.0;
         }//	w=w0+h/2.0

         for(j=1;j<=n;j++)
         {
            sum1=0.0;
            for(m=1;m<=deg1[j];m++)
               sum1=sum1+ax[B1[j][m]];
            l=sigma*(ay[j]-ax[j])+eps*(sum1-deg1[j]*ax[j])-dm[j]*ax[j];
            k[3][3*j-2]=l*h;
            l=ax[j]*(rho-z0[j])-x0[j]*az[j]-ay[j]-dm[j]*ay[j];
            k[3][3*j-1]=l*h;
            l=ax[j]*y0[j]+x0[j]*ay[j]-beta*az[j]-dm[j]*az[j];
            k[3][3*j]=l*h;
         }

         for(j=1;j<=n;j++)
         {
            ax[j]=ax0[j]+k[3][3*j-2];
            ay[j]=ay0[j]+k[3][3*j-1];
            az[j]=az0[j]+k[3][3*j];
         }//       w=w0+h

         for(j=1;j<=n;j++)
         {
            sum1=0.0;
            for(m=1;m<=deg1[j];m++)
               sum1=sum1+ax[B1[j][m]];
            l=sigma*(ay[j]-ax[j])+eps*(sum1-deg1[j]*ax[j])-dm[j]*ax[j];
            k[4][3*j-2]=l*h;
            l=ax[j]*(rho-z0[j])-x0[j]*az[j]-ay[j]-dm[j]*ay[j];
            k[4][3*j-1]=l*h;
            l=ax[j]*y0[j]+x0[j]*ay[j]-beta*az[j]-dm[j]*az[j];
            k[4][3*j]=l*h;
         }

         for(j=1;j<=n;j++)
         {
            ax0[j]=ax0[j]+(k[1][3*j-2]+2.0*(k[2][3*j-2]+k[3][3*j-2])+k[4][3*j-2])/6.0;
            ay0[j]=ay0[j]+(k[1][3*j-1]+2.0*(k[2][3*j-1]+k[3][3*j-1])+k[4][3*j-1])/6.0;
            az0[j]=az0[j]+(k[1][3*j]+2.0*(k[2][3*j]+k[3][3*j])+k[4][3*j])/6.0;
         }

         for(j=1;j<=n;j++)
         {
            v[(3*j-2)*3*n+1+j1]=ax0[j];
            v[(3*j-1)*3*n+1+j1]=ay0[j];
            v[(3*j)*3*n+1+j1]=az0[j];
         }
      }

      znorm[1]=0.0;
      for(j=1;j<=3*n;j++)
	 znorm[1]=znorm[1]+pow(v[3*n*j+1],2.0);	//	write(*,*)v(n*j+1)**2.0
      znorm[1]=sqrt(znorm[1]);
      for(j=1;j<=3*n;j++)
         v[3*n*j+1]=v[3*n*j+1]/znorm[1];
      for(j=2;j<=3*n;j++)
      {
         for(k2=1;k2<=j-1;k2++)
         {
            gsc[k2]=0.0;
            for(l1=1;l1<=3*n;l1++)
               gsc[k2]=gsc[k2]+v[3*n*l1+j]*v[3*n*l1+k2];
         }/**/
         for(k2=1;k2<=3*n;k2++)
         {
            for(l1=1;l1<=j-1;l1++)
            v[3*n*k2+j]=v[3*n*k2+j]-gsc[l1]*v[3*n*k2+l1];
         }
         znorm[j]=0.0;
         for(k3=1;k3<=3*n;k3++)
            znorm[j]=znorm[j]+pow(v[3*n*k3+j],2.0);
         znorm[j]=sqrt(znorm[j]);
	 for(k4=1;k4<=3*n;k4++)
	    v[3*n*k4+j]=v[3*n*k4+j]/znorm[j];
      }
      for(k5=1;k5<=3*n;k5++)
         cum[k5]=cum[k5]+log(znorm[k5])/log(2.0);	//	write(*,*)cum[k5];

   }
/////////////////////////////////////////////////////////////////////////////////////////////////
   w0=w0+h;
     
/*   if((g(x0,n)>500.0)||(g(y0,n)>500.0))
   {
      error=1000000.0;
      break;
   }
   if(ii>tn)
   {
      fprintf(fp1,"%f\t",w0); for(j=1;j<=n;j++) fprintf(fp1,"%f\t%f\t%f\t%f\t%f\t%f\t",x01[j],y01[j],z01[j],x02[j],y02[j],z02[j]);
      fprintf(fp1,"\n");
   }*/

}//Time Iteraction
t0=w0-tn*h;

for(i=1;i<=n;i++)
{
   for(j=i+1;j<=n;j++)
   {
      if(cum[j]>cum[i])
         {
            temp=cum[i];
            cum[i]=cum[j];
            cum[j]=temp;
         }
   }
}

printf("%d\t%f\t%f\t%f\t%f\n",i2,eps,eta,dr,cum[1]/t0);
fprintf(fp2,"%d\t%f\t%f\t%f\t%f\n",i2,eps,eta,dr,cum[1]/t0);

}

return 0;
fclose(fp1);
fclose(fp2);
}
