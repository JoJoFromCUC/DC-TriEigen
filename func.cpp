#include "func.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>

using namespace std;

double EPS = 1e-14;
double ZERO = 1e-8; 
double GAP = 1e-4; //deflate condition

//转置
void transpose(vector<vector<double> > &A,vector<vector<double> > &T){
    if(A.empty() || A[0].empty()) return;
    int m = A.size();
    int n = A[0].size();
    T.clear();
    T.resize(n,vector<double> (m,0));
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            T[j][i] = A[i][j];
        }
    }
}

//上三角转置
void transpose(vector<vector<double> > &A){
    int m = A.size();
    for(int i=0;i<m;i++){
        for(int j=i+1;j<m;j++){
            swap(A[i][j],A[j][i]);
        }
    }
}

//(level-3)reducer优化
void multiply(vector<vector<double> > &A,vector<vector<double> > &B,vector<vector<double> > &C){
    C.clear();
    if(A.empty() || A[0].empty() || B.empty() || B[0].empty()) return ;
    C.resize(A.size(),vector<double> (B[0].size(),0));
    //cilk::reducer< cilk::op_add<double> > tmp (0);
    for(int i=0;i<A.size();i++){
        for(int j=0;j<B[0].size();j++){
            //double tmp = 0;
            C[i][j] = 0;
            for(int k=0;k<A[0].size();k++){
                C[i][j] += A[i][k]*B[k][j];
            }
            // C[i][j] = tmp.get_value();
            // tmp.set_value(0);
        }
    }
}

//(level-2)
void multiply(const vector<vector<double> > &X,const vector<double> &v,vector<double> &res){
    res.clear();
    if(X.empty() || v.empty()) return;
    int m = X[0].size();
    res.resize(m,0);
    for(int i=0;i<m;i++){
        for(int j=0;X[i].size();j++){
            res[i] += X[i][j]*v[j];
        }
    }
}

//(level-1)
double dotProduct(const vector<double> &a, const vector<double> &b){
    double res=0;
    //cilk::reducer< cilk::op_add<double> >res (0);
    for(int i=0;i<a.size();i++){
        res += a[i]*b[i];
    }
    return res;
}

//res= A*A'*v
void rightMultiply(SparseMatrix &A, const vector<double> &v, vector<double> &res){
    int m=A.rows;
    int n=A.cols;
    res.clear();
    res.resize(m,0);
    vector<double> w(n,0);
    A.moveFirst();
    while(A.hasNext()){
        Cell c = A.next();
        w[c.col] += c.value*v[c.row];
    }
    A.moveFirst();
    while(A.hasNext()){
        Cell c=A.next();
        res[c.row]+=c.value*w[c.col];
    }
}

//res= A'*A*v
// void leftMultiply(SparseMatrix &A, const vector<double> &v, vector<double> &res){
//     int m=A.rows;
//     int n=A.cols;
//     res.clear();
//     res.resize(n,0);
//     vector<double> w(m,0);
//     A.moveFirst();
//     while(A.hasNext()){
//         Cell c=A.next();
//         w[c.row]+=c.value*v[c.col];
//     }
//     A.moveFirst();
//     while(A.hasNext()){
//         Cell c=A.next();
//         res[c.col]+=c.value*w[c.row];
//     }
// }

//C= B'*A
void rightMultiply(const vector<vector<double> > &B,SparseMatrix &A, vector<vector<double> > &C){
    int m=B[0].size();
    int k=B.size();
    int n=A.cols;
    for(int i=0;i<C.size();i++){
        fill(C[i].begin(),C[i].end(),0);
    }
    A.moveFirst();
    while(A.hasNext()){
        Cell c=A.next();
        for(int i=0;i<m;i++){
            C[c.col][i]+=c.value*B[c.row][i];
        }
    }
}

//C = A'*B
// void leftMultiply(const vector<vector<double> > &B,SparseMatrix &A, vector<vector<double> > &C){
//     int r=B[0].size();
//     int n=B.size();
//     int m=A.rows;
//     C.clear();
//     C.resize(m,vector<double>(r,0));
//     A.moveFirst();
//     while(A.hasNext()){
//         Cell c = A.next();
//         for(int i=0;i<r;i++){
//             C[c.row][i]+=c.value*B[c.col][i];
//         }
//     }
// }

//2范数
double norm(const vector<double> &v){
    //cilk::reducer< cilk::op_add<double> >r (0);
	double r = 0;
    for(int i=0;i<v.size();i++)
        r += v[i]*v[i];
    return sqrt(r);
}

//归一化
double normalize(vector<double> &v){
    //cilk::reducer< cilk::op_add<double> >r (0);
	double r = 0;
    for(int i=0;i<v.size();i++){
        r += v[i]*v[i];
    }
    double R = sqrt(r);
    if(R > EPS){
        for(int i=0;i<v.size();i++){
            v[i] /= R;
        }
    }
    return R;
}

//数乘
void multiply(vector<double> &v, double d){
    for(int i=0;i<v.size();i++){
        v[i] *= d;
    }
}

//random standard vector gen
void randUnitVector(int n, vector<double> &v){
    srand(time(NULL));
    v.clear();v.resize(n);
    while(true){
        double r=0;
        //cilk::reducer< cilk::op_add<double> >r (0);
        for(int i=1;i<n;i++){
            // v[i]= i % 10-5;
            v[i] = rand()*1.0/RAND_MAX - 0.5;
            r += v[i]*v[i];
        }
        r = sqrt(r);
        if(r > EPS){
            for(int i=0;i<n;i++){
                v[i]/=r;
            }
            break;
        }
    }
}

//print matrix
void print(vector<vector<double> > &X){
    cout.precision(6);
    cout.setf(ios::fixed);
    for(int i=0 ;i < X.size();i++){
        for(int j=0;j < X[i].size();j++){
            cout<<X[i][j]<<' ';
        }
        cout<<endl;
    }
    cout<<endl;
}
//deal with 0 problem
double prezero(double num){
	return num = num==0?ZERO:num;
}

//mpi收发数据
/*
bool read_data(double local_a[],int local_n,int n,string vec_name,int myrank,MPI_Comm comm){
    double* a = NULL;
    if(myrank == 0){
        a = new double[n];
        //准备数据
        MPI_Scatter(a,local_n,MPI_DOUBLE,local_a,local_n,MPI_DOUBLE,0,comm);
        delete [] a;
    }else{
        MPI_Scatter(a,local_n,MPI_DOUBLE,local_a,local_n,MPI_DOUBLE,0,comm);
    }
}
bool print_data(double local_b[],int local_n,int n,string title,int myrank,MPI_Comm comm){
    double *b = NULL;
    if(myrank == 0){
        b = new double[n];
        MPI_Gather(local_b,local_n,MPI_DOUBLE,b,local_n,MPI_DOUBLE,0,comm);
        //操作数据
        delete [] b;
    }else{
        MPI_Gather(local_b,local_n,MPI_DOUBLE,b,local_n,MPI_DOUBLE,0,comm);
    }
}
int main(){
    int myrank=0,comm_sz=0;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
    //maxn 数据总大小
    int local_n = maxn/comm_sz;
    double local_a[maxn]，local_b[maxn],local_c[maxn];
    MPI_Finalize();
    return 0;
}*/
