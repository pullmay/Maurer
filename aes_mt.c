#include<stdio.h>
#include<math.h>
#include"aes.h"
#include "kansuu.h"
#include "MT.h"

/* AESで生成した系列をMTで変換し，変換後の系列に対して検定を行う */
/* 最後にファイルを保存するので，ディレクトリを用意しておくこと */
/* 検定本数はarraysizeで指定 */

//MTによるフリッピング
void flipping(int n, double p, char *s){
  for (int i=0;i<n;i++){
    if (s[i]==1 && genrand_res53()<=p){
      s[i]=0;
    }
  }
}

//系列における1の個数をカウントする関数
int count(int n, char *s){
  int cnt = 0;
  for (int i=0;i<n;i++){
    if (s[i]==1) {
      cnt+=1;
    }
  }
  return cnt;
}

double maurerTest(int L, int Q, int K, double q, char *source){
  //返り値：f(x^n)
  int imax = 1000000;
  double hq = -q*log2(q)-(1-q)*log2(1-q);
  double f=0;//test statistic value
  unsigned char tmp;
  int i,j;
  int La = 1<<L;//2^L
  int x[La];


  static double g[1000001];
  make_g(g,imax);

  int mask=0;
  for(int i=0;i<L;i++){//下位L-bitを取り出すため
    mask=(mask<<1)+1;//00...011111...111（下位L-bitが1）
  }

  for(i=0;i<La;i++){x[i]=0;} //配列を0で初期化
  //-----init-----
  for (i=0; i<Q; i++)
  {
    tmp = 0;
    // tmp = genrand_int32()&mask; //in case of q=0.50
    // for(j=0;j<L;j++)
    // {
    //   if(genrand_res53()<q){tmp=(tmp<<1)+1;}else{tmp=(tmp<<1)+0;}
    // }
    for (j=0;j<L;j++){
      tmp=(tmp<<1)+source[i*L+j];
    }
    x[tmp] = i;
  }
  //-----test-----
  for (i=0; i<K; i++)
  {
    tmp = 0;
    // tmp = genrand_int32()&mask; //in case of q=0.50
    for(j=0;j<L;j++)
    {
      // if(genrand_res53()<q){tmp=(tmp<<1)+1;}else{tmp=(tmp<<1)+0;}
      tmp=(tmp<<1)+source[i*L+j];
    }
    f += g[i + Q - x[tmp]];
    x[tmp] = i + Q;
  }

  return f/K;

}

int main(){
  init_genrand(11111);//
  //-----parameter-----
  int L = 8;
  int Q = 10*pow(2,L);//10*2^L
  int K = 1000*pow(2,L);//1000*2^L
  double q = 0.33;
  double hq = -q*log2(q)-(1-q)*log2(1-q);
  int n = L*(Q+K);//length of sequence
  double f,p1,p2;
  int arraysize=300000;//--------------------------------------------------------------------検定本数
  double array_p_proposed[arraysize];
  double array_p_yamamoto[arraysize];
  //-----sigma-----
  double sigma_yamamoto = 0.00349225;
  double sigma_proposed = 0.0034886003389633727;
  double sigma = (0.3862500+(0.3640569)/K)*sqrt(3.3704039/K);//q=0.50


  //-----AES-----
  unsigned int plaintext[16];  //平文(128bit)
  unsigned int PT[16];  //平文(128bit)
  unsigned int ciphertext[16]; //暗号文(128bit)
  unsigned char key[16];       //鍵(128bit)
  char rand[n];//2068480
  unsigned char y[128];

  /*AESの仕様書では「鍵長は128bit,192bit,256bitから選択」となっているがこのプログラムではサポートしていない*/

  unsigned char key_exp[176];  //拡張鍵を入れとく変数


  /*以下サンプル*/


  //鍵を与える
  key[0]=0x00;
  key[1]=0x01;
  key[2]=0x02;
  key[3]=0x03;
  key[4]=0x04;
  key[5]=0x05;
  key[6]=0x06;
  key[7]=0x07;
  key[8]=0x08;
  key[9]=0x09;
  key[10]=0x0a;
  key[11]=0x0b;
  key[12]=0x0c;
  key[13]=0x0d;
  key[14]=0x0e;
  key[15]=0x0f;
  
  //鍵はあらかじめ拡張しておく
  keyexpancion(key,key_exp);
  /*鍵拡張は最初に1回しとけばよい。平文ごとにやる必要はない。*/


  //平文を与える
  PT[0]=0x00000000;
  PT[1]=0x00000000;
  PT[2]=0x00000000;
  PT[3]=0x00000000;
  

  //暗号化
  aes(key_exp,plaintext,ciphertext);

  for (int l=0; l<arraysize; l++){
  //1000回
    //---generation---
    for (int k=0; k<16160; k++){
      //平文の更新
      if (PT[3]==pow(2,32)-1){
        PT[2] += 1;
        PT[3] = 0;
      } else {
        PT[3] += 1;
      }
      //平文を暗号化
      aes(key_exp, PT, ciphertext);
      //ビットに変換
      for(int i=0;i<4;i++){
        for(int j=31;j>=0;j--){
          y[32*i+(31-j)]=(ciphertext[i]>>j)&1;//and演算
          // printf("%d",(ciphertext[i]>>j)&1);
        }
      }
      // printf("\n");
      //rand
      for (int m=0; m<128; m++){
        rand[k*128+m] = y[m];
        // printf("%d",rand[k+128+m]);
      }
    }

    //変換
    flipping(n, 1-2*q, rand);//---------------------flip
    f = maurerTest(L,Q,K,q,rand);
    p1 = erfc(fabs((f-L*hq)/(sqrt(2)*sigma_proposed)));
    p2 = erfc(fabs((f-L*hq)/(sqrt(2)*sigma_yamamoto)));
    array_p_proposed[l] = p1;
    array_p_yamamoto[l] = p2;
    printf("now:%d\n",l);
  }

  // 配列をファイルに書き込み
  // FILE *fp;
  // FILE *gp;
  // fp=fopen("./file/proposed_100000","wb");
  // gp=fopen("./file/yamamoto_100000","wb");
  // fwrite(array_p_proposed, 8, arraysize, fp);
  // fwrite(array_p_yamamoto, 8, arraysize, gp);
  // fclose(fp);
  // fclose(gp);

  return 0;

}
