#include <iostream>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iomanip>
using namespace std;

class householder
{
private:
 int i, j, k, m, n, sgn;
 double eps, s, sum;
 double *w;
 double **a, **b, **c,**p;
//p保留最终正交矩阵
public:
 householder()
 {
  eps = 1e-10;
 }
 void solution();
 double compute_s();
 int sign();
 ~householder()
 {
  delete[] w;
  for (i = 0; i < n; i++)
  {
   delete[] a[i];
  }
  delete[] a;
  for (i = 0; i < n; i++)
  {
   delete[] b[i];
  }
  delete[] b;
  for (i = 0; i < n; i++)
  {
   delete[] c[i];
  }
  delete[] c;
 }
};
void multiply(double **a,double **b,double **c,int n,int m){
	int i,j,k;
	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
		 c[i][j] = 0.0;
			 for (k = 0; k < m; k++){
			  c[i][j] += a[i][k] * b[k][j];
			 }
		}   
   }
}
int main()
{
	cout.precision(15);
	cout.setf(ios::fixed);
	householder tridiagonal;
	tridiagonal.solution();
	return 0;
}

void householder::solution(){
 cout << "input dimension of matrix:";
 cin >> n;
 w = new double[n];
 a = new double*[n];
 for (i = 0; i < n; i++)
 {
  a[i] = new double[n];
 }
 b =  new double*[n];
 for (i = 0; i < n; i++)
 {
  b[i] = new double[n];
 }
 c =  new double*[n];
 for (i = 0; i < n; i++)
 {
  c[i] = new double[n];
 }
 p = new double*[n];
 for(i=0;i<n;i++){
	 p[i]=new double[n];
 }
for (i = 0; i < n; i++){
	for (j = 0; j < n; j++){
		if (i == j){
			p[i][j] = 1.0;
		}
		else{
			p[i][j] = 0.0;
		}
	}   
}
FILE *fp = fopen("./ctest/data500.txt","r");
if(!fp){
	cout<<" open matrix file failed!"<<endl;
}
int tmp;
for(int i=0;i<n;i++){
	for(int j=0;j<n;j++){
		if(fscanf(fp,"%d",&tmp) == EOF) break; 
		a[i][j] = tmp;
	}
	if(feof(fp)) break;
}
cout<<"data loaded ."<<endl;
fclose(fp);

cout<<"start calculating ... "<<endl;
for (k = 0; k < (n-2); k++){
	cout<<"k: "<<k<<endl;
	s = compute_s();
	if (s > eps){
		for (i = 0; i < n; i++){
			w[i] = 0.0;
		}
		w[k+1] = sqrt(0.5 * (1 + fabs(a[k][k+1]) / s));
		for (i = (k+2); i < n; i++){
			w[i] = (a[k][i] / (2 * s * w[k+1])) * sign();
		}
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++)
			{
				c[i][j] = w[i] * w[j];
			}
		}
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++){
				if (i == j){
					b[i][j] = 1.0;
				}
				else{
					b[i][j] = 0.0;
				}
			}   
		}
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++){
			  b[i][j] -= 2 * c[i][j];
			}   
		}//正交矩阵生成
		for (i = 0; i < n; i++){
			for (j = 0; j < n; j++){
				double temp = 0;
				 for (m = 0; m < n; m++){
			    	  temp += b[i][m] * p[m][j];
				 }
				 p[i][j] = temp;
			}   
		} 
		multiply(a, b, c, n, n);//c = a*b
		multiply(b, c, a, n, n);//a = b*c
	}
 }
	FILE *fr = fopen("./p500.txt","w");//写入正交矩阵
	for(int i=0;i<n;++i){
		for(int j=0;j<n;++j){
			fprintf(fr,"%.15lf ",p[i][j]);
		}
		fprintf(fr,"%s","\n");
	}
	fclose(fr);
	cout<<"正交矩阵已写入p.txt文件!\n"<<endl; 
	
	cout << "the result is :" << endl;
	FILE *ft = fopen("./ab500.txt","w");//写入三对角矩阵元素
	for(int i=0;i<n;++i){
		for(int j=0;j<n;++j){
			if(j==i) fprintf(ft,"%.15lf ",a[i][j]);
		}
	}
	fprintf(ft,"%.15lf ",0);
	for(int i=0;i<n;++i){
		for(int j=0;j<n;++j){
			if(j==i+1) fprintf(ft,"%.15lf ",a[i][j]);
		}
	}
	fclose(ft);
	cout<<"三对角元素已写入ab.txt文件!\n"<<endl; 
	char stop = getchar();
}

double householder::compute_s()
{
 sum = 0.0;
 for (i = (k+1); i < n; i++)
 {
  sum += a[k][i] * a[k][i];
 }
 sum = sqrt(sum);
 return(sum);
}

int householder::sign()
{
 if (fabs(a[k][k+1]) < eps)
 {
  sgn = 1;
 }
 else if (a[k][k+1] < 0)
 {
  sgn = -1;
 }
 else
 {
  sgn = 1;
 }
 return(sgn);
}

