#include "func.h"
#ifndef DC_H
#define DC_H
using namespace std;

vector<double> secularEquationSolver(vector<double> &z, vector<double> &D, double sigma,int start,int end);

void DCSub(vector<double> &alpha, vector<double> &beta, vector<vector<double> > &Q, vector<double> &D, int start, int end);

void DCTridiagonal(vector<double> alpha, vector<double> &beta, vector<vector<double> > &Q, vector<double> &D);

void resolve(SparseMatrix &A, int r,vector<double> &alpha,vector<double> &beta);

#endif