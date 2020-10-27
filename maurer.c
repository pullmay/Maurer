#include <stdio.h>
#include <math.h> 
#include "kansuu.h"
#include "MT.h"


int main()
{
  int seed = 111;
  init_genrand(seed);//メルセンヌツイスタの初期化 //プログラムごとの変える
  //parameter
  int L = 8;
  int Q = 10*pow(2,L);//10*2^L
  int K = 1000*pow(2,L);//1000*2^L
  int imax = 1000000;
  double q = 0.33;
  double hq = -q*log2(q)-(1-q)*log2(1-q);

  //関数gはあらかじめ計算しておいて配列に。g(0)からg(imax)まで。
  static double g[1000001];
  make_g(g, imax);

  //
  int length_of_data = L*(Q+K);
  int block_size = L;
  int num_block = Q + K;
  int num_init_block = Q;
  int num_test_block = K;

  int La = 1<<L;//2^L
  // int x[La];//table
  double f;//test function
  unsigned char tmp;

  int M=4000000;  //不偏分散を計算のための検定統計量fの数
  int N=1;        //不偏分散の平均を計算するための不偏分散の数
  double mu;
  double var;
  double array_mu[N+1];
  double array_var[N+1];

  //make mask
  int mask=0;
  for(int i=0;i<L;i++){//下位L-bitを取り出すため
    mask=(mask<<1)+1;//00...011111...111（下位L-bitが1）
  }

  // double array_f[N+1][M+1]={}; //0で初期化 //8GB
  double array_f1[N+1];
  double array_f2[N+1];
  // double var_f[N+1];
  for (int n=0; n<N; n++)
  {
    // double array_f[M+1]={}; //配列を0で初期化
    mu = 0.0; //初期化
    var = 0.0; //初期化
    int i,j;
    int x[La];
    for (int m=0; m<M; m++)
    {
      // printf("%d/%d out of %d/%d\n",m,M,n,N);
      f = 0.0; //初期化
      
      for(i=0;i<La;i++){x[i]=0;} //配列を0で初期化
      //init
      for (i=0; i<num_init_block; i++)
      {
        tmp = 0;
        // tmp = genrand_int32()&mask; //in case of q=0.50
        for(j=0;j<L;j++)
        {
          if(genrand_res53()<q){tmp=(tmp<<1)+1;}else{tmp=(tmp<<1)+0;}
        }
        x[tmp] = i;
      }

      //test
      for (i=0; i<num_test_block; i++)
      {
        tmp = 0;
        // tmp = genrand_int32()&mask; //in case of q=0.50
        for(j=0;j<L;j++)
        {
          if(genrand_res53()<q){tmp=(tmp<<1)+1;}else{tmp=(tmp<<1)+0;}
        }
        f += g[i + num_init_block - x[tmp]];
        x[tmp] = i + num_init_block;
      }
      array_f1[n] += (f/num_test_block)/M; 
      array_f2[n] += pow(f/num_test_block,2.0)/(M-1.0);
    }

    array_var[n]=array_f2[n]-M/(M-1.0)*pow(array_f1[n],2.0); //variance of f
  }

  printf("-----------------\n");
  
  double var_average = 0.0;
  for (int n=0; n<N; n++)
  {
    printf("%.16e\n", array_var[n]); //不偏分散(N個)
    var_average += array_var[n];
  }
  var_average = var_average/N; //不偏分散の平均値
  printf("-----------------\n");
  printf("The average of the Unbiased Variance of this Simulation : %.16e\n", var_average);
  printf("-----------------\n");
  printf("L=%d\n",L);
  printf("Q=%d\n",Q);
  printf("K=%d\n",K);
  printf("q=%f\n",q);
  printf("seed=%f\n",seed);
  printf("M=%d\n",M);
  printf("N=%d\n",N);


  // 配列をファイルに書き込み
  // 不偏分散(N個)をファイルに書き込む
  // FILE *fp;
  // fp=fopen("var_mersenne_seed_1","wb");
  // fwrite(array_var, 8, N, fp);
  // fclose(fp);

}
