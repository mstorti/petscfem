#include <cstdio>
extern "C" void dsyev2_(int *N,double *a, double *w, double *work,int *lwork,int *info);

//  void dsyev_(int *jobz,int *uplo,int *n,double *a,int *lda,
//  	    double *w,double *work, int *lwork, int *info);

int main() {
  int N=3;
  double a[N][N],w[N],work[5*N];
  for (int j=0; j<N; j++) {
    for (int k=0; k<N; k++) a[j][k] = 0.01;
    a[j][j] = 1.;
  }
  int info,lwork = 5*N;
  dsyev2_(&N,&a[0][0],w,work,&lwork,&info);
  for (int j=0; j<N; j++) printf("%f\n",w[j]);
}
