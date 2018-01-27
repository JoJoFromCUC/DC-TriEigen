﻿/*
author: JoeZhu
time: Oct 2nd,2017
institution: UCAS
license：GPL 3.0
*/
#include <iostream>
#include <time.h>
#include <fstream>
#include "func.h"
#define DIM 1000

int main(){
    cout.precision(8);
    cout.setf(ios::fixed);
    vector<vector<double> > A(DIM,vector<double>(DIM,0));
    // A.rows = DIM;
    // A.cols = DIM;
    int rank = DIM;
	//int len = 2*DIM-1;
	
    //读入矩阵数据
	FILE *fp;
	fp = fopen("./ctest/data1000.txt","r");
	if(!fp){
		cout<<"open file failed!"<<endl;
		return -1;
	}
	static int A_rows = A.size();
	static int A_cols = A[0].size();
    for(int i=0;i<A_rows;i++){
		int tmp;
        for(int j=0;j<A_cols;j++){
			if(fscanf(fp,"%d",&tmp) == EOF) break; 
            A[i][j] = tmp;
        }
		if(feof(fp)) break;
    }
	fclose(fp);  
	
	/* cout<<"A = "<<endl;//输出A
    A.moveFirst();
    for(int i=0;i<A.rows;i++){
        for(int j=0;j<A.cols;j++){
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
	cout<<endl;  */
	
	vector<vector<double> > Q;
    vector<double> D; 
	vector<double> alpha(DIM);
	vector<double> beta(DIM-1);
	
	alpha.clear();
	beta.clear();
	for(int i=0;i<DIM;i++){
		for(int j=0;j<DIM;j++){
			if(j == i) alpha.push_back(A[i][j]);
			else if(j == i+1) beta.push_back(A[i][j]); 
		}
	} 
	
	/* ifstream in("./ctest/data1.txt");//读取主对角、副对角元素
	if(!in){
		cout<<"文件打开失败！"<<endl;
		return -1;
	}
	D.clear();
	int temp;
	while(in >> temp){
		D.push_back(temp);
	}
	in.close();
	alpha.clear();
	beta.clear();
	for(int i=0;i<DIM-1;i++){
		alpha.push_back(D[i]);
		beta.push_back(D[i+DIM]);
	}
	alpha.push_back(D[DIM-1]); 
	D.clear();
	D.resize(DIM);
	*/
	// cout<<"alpha:"<<endl;
	// for(int i=0;i<DIM;i++){
	// 	cout<<alpha[i]<<" ";
	// }
	// cout<<endl<<"beta:"<<endl;
	// for(int i=0;i<DIM-1;i++){
	// 	cout<<beta[i]<<" ";
	// }
	cout<<endl<<"data loadig ..."<<endl<<"alpha.size: "<<alpha.size()<<endl
		<<"beta.size: "<<beta.size()<<endl;
		
	clock_t start,finish;//计时
	start = clock();
	
	cout<<"rank = "<<rank<<endl;
    cout<<endl;
	cout<<"start calculating ..."<<endl;

	//DCTridiagonal(alpha,beta,Q,D);
	
	vector<vector<double> > U,V;
    vector<double> s;
	resolve(A,rank,U,s,V);//开始svd分解 
	
	cout<<endl;
	cout<<"eigen values :"<<endl;//输出特征值
    for(int i=0;i<s.size();i++){
        for(int j=0;j<s.size();j++){
            if(i == j) cout<<s[i]<<' ';
        }
        cout<<endl;
    }
    cout<<endl;
	
	// fp = fopen("./result.txt","w");//写入特征向量矩阵
	// int rows = U.size(),cols = U[0].size();
	// fprintf(fp,"%s","特征向量矩阵为:\n");
	// for(int i=0;i<rows;++i){
	// 	for(int j=0;j<cols;++j){
	// 		fprintf(fp,"%f ",U[i][j]);
	// 	}
	// 	fprintf(fp,"%s","\n");
	// }
	// fclose(fp);
    // cout<<"特征向量矩阵已写入result.txt文件!\n"<<endl; 
	
    //print(U);

	/* fp = fopen("./ctest/result.txt","a");//写入矩阵V
	rows = V.size(),cols = V[0].size();
	fprintf(fp,"%s","V矩阵为:\n");
	for(int i=0;i<rows;++i){
		for(int j=0;j<cols;++j){
			fprintf(fp,"%f ",V[i][j]);
		}
		fprintf(fp,"%c",'\n');
	}
	fclose(fp);
    cout<<"V 矩阵已写入result.txt文件!"<<endl; */
    //print(V);
	
	finish = clock();
	cout<<"the time wasted is:"<<(double)(finish-start)/CLOCKS_PER_SEC<<" s"<<endl;
	char stop = getchar();
    return 0;
}