#include <cstdio>

extern "C" void dsyev_(const char*jobz,const char *uplo,int *n,double *a,int *lda,
		       double *w,double *work, int *lwork, int *info);

int main() {
  int N=3;
  double a[N][N],w[N],work[5*N];
  for (int j=0; j<N; j++) {
    for (int k=0; k<N; k++) a[j][k] = 0.01;
    a[j][j] = 1.;
  }
  int info,lwork = 5*N;
  dsyev_("V","U",&N,&a[0][0],&N,w,work,&lwork,&info);
  for (int j=0; j<N; j++) printf("%f\n",w[j]);
  for (int j=0; j<N; j++) {
    for (int k=0; k<N; k++) 
      printf("%f ",a[k][j]);
    printf("\n");
  }
}

