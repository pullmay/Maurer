#include<stdio.h>
#include<math.h>  
#include<time.h>  
#include<omp.h>
#include"kansuu.h"

/* コンパイル方法 */
/* g++ -fopenmp final.c -o final */


/* パラメータを定義 */
#define L 4        //4,5,6,7,8
#define q 0.4      //(0,1)
#define Kmax 16000   //1000*2^L
#define imax 1000000 //1 million（無限和をどこまで計算するかを表す）
#define jmax imax   //i million as the same as imax 

int main(){

  //-----------------------------------------------------------------------------// 

  clock_t start,end;

  int i,j,k,K,r,r1,r2;
  double tmp1,tmp2,tmp3;

  /* 時間計測 */
  printf("timer start\n");
  start = clock();

  //関数gはあらかじめ計算しておいて配列に入れておく．g(0)からg(imax)まで
  static double g[1000001];//要素数: imax+1=1000001 
  make_g(g,imax);

  //周辺分布の配列，および，結合分布を計算する関数のための配列を用意しておく
  //二項係数
  int LC[L+1];
  make_LC(LC,L);
  //pa = q^r * (1-q)^(L-r)
  double pa[L+1];
  make_parameter(pa,q,L);
  //marginal_pmf = \sum_{i=1}^{imax} \binom{L}{r} Q(r)^2 * (1-Q(r))^{i-1}
  static double marginal_pmf[1000001];//imax+1
  make_marginal_pmf(marginal_pmf,L,q,pa,LC,imax);

  //-----------------------------------------------------------------------------//

  //Var[g]
  printf("calculating var_g\n");
  double var_g;
  tmp1=0.0;
  tmp2=0.0;
  for(i=imax;i>=1;i--){
    tmp1+=g[i]*g[i]*marginal_pmf[i];
    tmp2+=g[i]*marginal_pmf[i];
  }
  var_g=tmp1-tmp2*tmp2;

  end = clock();
  printf("%.2f秒かかりました\n",(double)(end-start)/CLOCKS_PER_SEC);
  printf("var_g: %f\n", var_g);

  //-----------------------------------------------------------------------------//

  //case1
  double J1[Kmax+1];
  J1[0]=0.0;
  J1[1]=0.0;
  J1[2]=0.0;
  for(k=3;k<=Kmax;k++){
    J1[k]=0.0;
    for(r=0;r<=L;r++){
      J1[k]+=LC[r]*pa[r]*pa[r]*pow(1.0-pa[r],k-3)*g[k-2];
    }
    J1[k]+=J1[k-1];
  }
  double I1_const=0.0;
  for(r=0;r<=L;r++){
    I1_const+=LC[r]*pa[r]*(-1.0/log(2.0))*log(pa[r]);
  }

  //case2
  double IJ2[Kmax+1];
  IJ2[0]=0.0;
  IJ2[1]=0.0;
  for(k=2;k<=Kmax;k++){
    IJ2[k]=0.0;
    for(r=0;r<=L;r++){
      IJ2[k]+=-g[k-1]/log(2)*LC[r]*pa[r]*pa[r]*pow(1.0-pa[r],k-2)*log(pa[r]);
    }
  }

  //case3
  static double I3[9][1000001];
  for(r=0;r<=L;r++){
    I3[r][0]=0.0;
    I3[r][imax]=0.0;
    for(i=imax-1;i>=1;i--){//大きい方から
      I3[r][i]=I3[r][i+1]+g[i]*pow(1.0-pa[r],i);
    }
  }
  
  
  //case5
  static double J5[9][1200001];//kmax+imax(=1000*2^L + 1000000) //0で初期化される
  for(r=0;r<=L;r++){
    J5[r][0]=0.0;
    J5[r][jmax]=0.0;
    for (j=1200001;j>jmax;j--){
      J5[r][j] = 0.0;
    }
    for(j=jmax-1;j>=1;j--){//大きい方から
      J5[r][j]=J5[r][j+1]+g[j]*pow(1.0-pa[r],j);
    }
  }

  // for (j=jmax-100; j<=jmax+100; j++){
  //   printf("J5[0][%d]:%0.16e\n",j, J5[0][j]);
  // }

  printf("--------------- start calculating Covariance ---------------\n");

  // Cov[g(A_1),g(A_k)]=cov[k]の計算
  double cov[Kmax+1];//並列化する
  //#pragma omp parallel
  //{
    double tmp1omp,tmp2omp,tmp3omp;
    int iomp,jomp,r1omp,r2omp,romp;
    #pragma omp parallel for private(tmp1omp,tmp2omp,tmp3omp,iomp,jomp,r1omp,r2omp,romp)
    //並列化する
    for(k=Kmax;k>=2;k--){

      tmp1omp=0.0;

      //case1
      printf("case 1\n");
      tmp1omp+=J1[k]*I1_const;

      //case2
      printf("case 2\n");
      tmp1omp+=IJ2[k];
   
      //case3
      printf("case 3\n");
      tmp3omp=0.0;
      //#pragma omp parallel for private(tmp2omp,r1omp,r2omp) reduction(+:tmp3omp)
      for(jomp=jmax;jomp>=k;jomp--){
        tmp2omp=0.0;
        for(r1omp=0;r1omp<=L;r1omp++){
          for(r2omp=0;r2omp<=L;r2omp++){
              tmp2omp+=LC[r1omp]*LC[r2omp]*pa[r1omp]*pa[r1omp]*pa[r2omp]*pa[r2omp]*pow(1.0-pa[r1omp],-2)*pow(1.0-pa[r2omp]/(1.0-pa[r1omp]),jomp-k)*pow(1.0-pa[r2omp],k-2)*g[jomp]*I3[r1omp][jomp-k+2];
          }
          tmp2omp-=LC[r1omp]*pa[r1omp]*pa[r1omp]*pa[r1omp]*pa[r1omp]*pow(1.0-pa[r1omp],-2)*pow(1.0-pa[r1omp]/(1.0-pa[r1omp]),jomp-k)*pow(1.0-pa[r1omp],k-2)*g[jomp]*I3[r1omp][jomp-k+2];
        }
        tmp3omp = tmp3omp + tmp2omp;
      }
      tmp1omp+=tmp3omp;
      

      //case4
      printf("case 4\n");
      tmp1omp+=0.0;

      //case5
      printf("case 5\n");
      tmp3omp=0.0;
      //#pragma omp parallel for private(tmp2omp,r1omp,r2omp) reduction(+:tmp3omp)
      for(iomp=imax;iomp>=1;iomp--){
        tmp2omp=0.0;
        for(r1omp=0;r1omp<=L;r1omp++){
          for(r2omp=0;r2omp<=L;r2omp++){
            tmp2omp+=LC[r1omp]*LC[r2omp]*pa[r1omp]*pa[r1omp]*pa[r2omp]*pa[r2omp]*pow(1.0-pa[r1omp],-3)*g[iomp]*pow(1.0-pa[r2omp]/(1.0-pa[r1omp]),iomp-1)*J5[r1omp][k+iomp];
          }
          tmp2omp-=LC[r1omp]*pa[r1omp]*pa[r1omp]*pa[r1omp]*pa[r1omp]*pow(1.0-pa[r1omp],-3)*pow(1.0-pa[r1omp]/(1.0-pa[r1omp]),iomp-1)*g[iomp]*J5[r1omp][k+iomp];
        }
        tmp3omp = tmp3omp + tmp2omp;
      }
      tmp1omp+=tmp3omp;
      
      cov[k]=tmp1omp-L*L*H(q)*H(q);
      printf("%d %0.16e\n",k,cov[k]);
    }
  //}

  printf("--------------- finish calculating Covariance ---------------\n");

  //Var[f]の計算
  //var_f[0]=var_f[1]=0.0
  double var_f[Kmax+1];
  for(K=2;K<=Kmax;K++){
    tmp1=0.0;
    for(k=2;k<=K;k++){
      tmp1+=(K-(k-1))*cov[k];
    }
    var_f[K]=var_g/K+2.0*tmp1/(K*K);
    // printf("%0.16e\n",var_f[K]);
  }

  //-----------------------------------------------------------------------------//

  // 配列をファイルに書き込み
  // FILE *fp;
  // fp=fopen("./cov_data/cov_q040L4_final","wb");
  // fwrite(cov,8,Kmax+1,fp);
  // fclose(fp);

  // fp=fopen("./var_data/var_q040L4_final","wb");
  // fwrite(var_f,8,Kmax+1,fp);
  // fclose(fp);

}