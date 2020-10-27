#include<stdio.h>
#include<math.h> 

void make_g(double *g, int imax){
  g[0]=0.0;
  g[1]=0.0;
  for(int i=2;i<=imax;i++){
    g[i]=g[i-1]+1.0/(i-1.0);
  }
  for(int i=2;i<=imax;i++){
    g[i]/=log(2.0);
  }
}

void make_LC(int *LC, int L){
  for(int r=0;r<=L;r++){
    int tmp=1;
    for(int i=1;i<=L;i++){
      tmp*=i;
    }
    for(int i=1;i<=r;i++){
      tmp/=i;
    }
    for(int i=1;i<=L-r;i++){
      tmp/=i;
    }
    LC[r]=tmp;
  }
}

void make_parameter(double *parameter, double q, int L){
  for(int r=0;r<=L;r++){
    parameter[r]=pow(q,r)*pow(1.0-q,L-r);
  }
}


void make_marginal_pmf(double *marginal_pmf, int L, double q,double *parameter, int *LC, int imax){
  marginal_pmf[0]=0.0;
  double tmp;
  for(int i=1;i<=imax;i++){
    tmp=0.0;
    for(int r=0;r<=L;r++){
      tmp+=LC[r]*pow(parameter[r],2.0)*pow(1.0-parameter[r],i-1.0);
    }
    marginal_pmf[i]=tmp;
  }
}

  double H(double q){
    return -q*log2(q)-(1.0-q)*log2(1.0-q);
  }
  
