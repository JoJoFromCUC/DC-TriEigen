/*
author: JoeZhu
time: Oct 2nd,2017
institution: UCAS
license：GPL 3.0
*/
#include <iostream>
#include <time.h>
#include <fstream>
#include <iomanip>
#include "func.h"
#define DIM 50
using namespace std;

int main(int argc,char *argv[]){
    cout.precision(15);
    cout.setf(ios::fixed);
    SparseMatrix A;
    A.rows = DIM;
    A.cols = DIM;
    int rank = DIM;
    int len = 2*DIM-1;

	// int com_size = 0,myrank =0;//MPI初始化
    // MPI_Status status;
    // MPI::Init();
    // MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    // MPI_Comm_size(MPI_COMM_WORLD,&com_size);
	//读入矩阵数据
	cout<<"loading data..."<<endl;
	FILE *fp;
	fp = fopen("./ctest/data50.txt","r");
	if(!fp){
		cout<<" open matrix file failed!"<<endl;
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
	fclose(fp);
	/* cout<<"A = "<<endl;//输出A
    	A.moveFirst();
    	for(int i=0;i<A.rows;i++){
        	for(int j=0;j<A.cols;j++){
           		cout<<A.next().value<<" ";
        	}
        	cout<<endl;
    	}
		cout<<endl;  */
	vector<double> alpha(DIM,0);
	vector<double> beta(DIM,0);
	vector<vector<double> > Q;
    vector<double> D; 
	ifstream in("./ctest/ab50.txt");//读取主对角、副对角元素
	if(!in){
		cout<<" open tri elements file failed!"<<endl;
		return -1;
	}
	D.clear();
	double temp;
	while(in >> temp){
		D.push_back(temp);
	}
	in.close();
	cout<<"D size:"<<D.size()<<endl;;
	for(int i=0;i<DIM;i++){
		alpha[i] = D[i];
		beta[i] = D[i+DIM];
	}
	//alpha[DIM-1] = D[DIM-1];
	D.clear();
	cout<<beta[0]<<endl<<beta[DIM-1]<<endl;
	cout<<"alpha.size: "<<alpha.size()<<endl
		<<"beta.size: "<<beta.size()<<endl;	
	//DCSub(alpha, beta, Q, D, 0, m-1);
	clock_t start,finish;//计时
	start = clock();	
	cout<<"rank = "<<rank<<endl;
	cout<<endl;
	cout<<"start calculating ..."<<endl;

	resolve(A,rank,alpha,beta);//开始分解 
	
	finish = clock();
	cout<<"used time:"<<(double)(finish-start)/CLOCKS_PER_SEC<<" s"<<endl;
	char stop = getchar();
	//MPI::Finalize();
    return 0;
}
