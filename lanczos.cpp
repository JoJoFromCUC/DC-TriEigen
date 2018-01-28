#include "func.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <fstream>

void lanczos(SparseMatrix &A, vector<vector<double> > &P, vector<double> &alpha, vector<double> &beta, unsigned int rank){
    //P'*A*A'*P = T = diag(alpha) + diag(beta,1) + diag(beta, -1)
    //P=[p1,p2, ... , pk]
    rank=min(A.cols,min(A.rows,rank));
    vector<double> p;
    unsigned int m=A.rows;
    unsigned int n=A.cols;
    vector<double> prevP(m,0);
    randUnitVector(m,p);	//生成随机归一化P向量
    P.clear();
    P.resize(m,vector<double>(rank,0));
    vector<double> v;
    alpha.clear();
	alpha.resize(rank);
    beta.clear();
	beta.resize(rank);
    beta[0]=0; //补0使协调
    for(int i=0;i<rank;i++){
        for(int j=0;j<p.size();j++){
            P[j][i]=p[j];	//P的每一列都为p
        }
        rightMultiply(A, p, v);	//v=A*A‘*P
        alpha[i]=dotProduct(p,v);
        if(i+1<rank){
            for(int j=0;j<m;j++){
                v[j]=v[j]-beta[i]*prevP[j]-alpha[i]*p[j];
            }
            beta[i+1]=norm(v);
			beta[i+1] = prezero(beta[i+1]);
            prevP=p;
            for(int j=0;j<m;j++){
                p[j]=v[j]/beta[i+1];
            }
        }
    }
    FILE *fk = fopen("./ctest/P3000.txt","w");//写入正交矩阵P
	for(int i=0;i<m;++i){
		for(int j=0;j<rank;++j){
			fprintf(fk,"%.16lf ",P[i][j]);
		}
		fprintf(fk,"%s","\n");
	}
	fclose(fk);
    FILE *fab = fopen("./ctest/ab3000.txt","w");//写入alpha
	for(int i=0;i<m;++i){
		fprintf(fab,"%.16lf ",alpha[i]);	
	}
    for(int j=0;j<m;++j){
        fprintf(fab,"%.16lf ",beta[j]);//写入beta
    }
	fclose(fab);
    cout<<"分解结果已写入文件\n"<<endl;
}
int main(){
    int DIM = 3000;
    SparseMatrix A;
    A.rows = DIM;
    A.cols = DIM;
    cout<<"loading data..."<<endl;
    FILE *fp;
    fp = fopen("./ctest/data3000.txt","r");
    if(!fp){
        cout<<"open file failed!"<<endl;
        return -1;
    }
    for(int i=0;i<A.rows;i++){
        int tmp;
            for(int j=0;j<A.cols;j++){
                if(fscanf(fp,"%d",&tmp) == EOF) break; 
                    A.cells.push_back(Cell(i,j,tmp));
            }
        if(feof(fp)) break;
    }
    vector<vector<double> > P(DIM,vector<double>(DIM,0));
    vector<double> alpha;
    vector<double> beta;
    lanczos(A,P,alpha,beta,DIM);//三对角化
    fclose(fp);
    return 0;
}
