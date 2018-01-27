#include "func.h"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
using namespace std;

const double EPS = 1e-14;//精度
const double ZERO = 1e-8;//0

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

//(level-3)
void multiply(vector<vector<double> > &A,vector<vector<double> > &B,vector<vector<double> > &C){
    C.clear();
    if(A.empty() || A[0].empty() || B.empty() || B[0].empty()) return ;
    int m = A.size();
    int n = A[0].size();
    int l = B[0].size();
     C.resize(m,vector<double> (l,0));
    for(int i=0;i<m;i++){
        for(int j=0;j<l;j++){
            C[i][j]=0;
            for(int k=0;k<n;k++){
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
}

//(level-2)
void multiply(const vector<vector<double> > &X,const vector<double> &v,vector<double> &res){
    res.clear();
    if(X.empty() || v.empty()) return;
    int m = X.size();
    int n = X[0].size();
    res.resize(m,0);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            res[i] += X[i][j]*v[j];
        }
    }
}

//(level-1)
double dotProduct(const vector<double> & a, const vector<double> & b){
    double res=0;
    int n = a.size();
    for(int i=0;i<n;i++)
        res+=a[i]*b[i];
    return res;
}

//res= A*A'*v
void rightMultiply(vector<vector<double> > &A, const vector<double> & v, vector<double> & res){
    int m=A.size();
    int n=A[0].size();
    res.clear();
    res.resize(m,0);
    // vector<double> w(n,0);
    // A.moveFirst();
    // while(A.hasNext()){
    //     Cell c = A.next();
    //     w[c.col] += c.value*v[c.row];
    // }
    // A.moveFirst();
    // while(A.hasNext()){
    //     Cell c=A.next();
    //     res[c.row]+=c.value*w[c.col];
    // }
    vector<vector<double> > T;//转置
    vector<vector<double> > C;
    transpose(A,T);
    multiply(A,T,C);
    multiply(C,v,res);
}

//res= A'*A*v
// void leftMultiply(vector<vector<double> > &A, const vector<double> & v, vector<double> & res){
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
void rightMultiply(vector<vector<double> > &B,vector<vector<double> > &A, vector<vector<double> > &C){
    int m=B[0].size();
    int k=B.size();
    //int n=A.cols;
    //for(int i=0;i<C.size();i++) fill(C[i].begin(),C[i].end(),0);
    // A.moveFirst();
    // while(A.hasNext()){
    //     Cell c=A.next();
    //     for(int i=0;i<m;i++){
    //         C[c.col][i]+=c.value*B[c.row][i];
    //     }
    // }
    vector<vector<double> > T;//转置
    transpose(B,T);
    multiply(T,A,C);
}

//C = A'*B
// void leftMultiply(const vector<vector<double> > & B,SparseMatrix &A, vector<vector<double> > & C){
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
    double r = 0;
    for(int i=0;i<v.size();i++)
        r += v[i]*v[i];
    return sqrt(r);
}

//归一化
double normalize(vector<double> &v){
    double r = 0;
    for(int i=0;i<v.size();i++)
        r += v[i]*v[i];
    r = sqrt(r);
    if(r > EPS){
        for(int i=0;i<v.size();i++)
            v[i] /= r;
    }
    return r;
}

//数乘
void multiply(vector<double> &v, double d){
    for(int i=0;i<v.size();i++)
        v[i] *= d;
}

//随机单位向量
void randUnitVector(int n, vector<double> &v){
    v.clear();v.resize(n);
    while(true){
        double r=0;
        for(int i=0;i<n;i++){
            v[i]= i % 8;
            r+=v[i]*v[i];
        }
        r=sqrt(r);
        if(r>EPS){
            for(int i=0;i<n;i++)
                v[i]/=r;
            break;
        }
    }
}

//打印输出
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
//避免除0操作
double prezero(double num){
	return num = num==0?ZERO:num;
}

//特征方程求解
vector<double> secularEquationSolver(vector<double> &z, vector<double> &D, double sigma,int start,int end){
 
    int n=z.size();
    vector<double> res(n);
    //sort : d_0 < d_1 < ... < d_{n-1}
    vector<int> index;
    vector<double> d(n);
    merge_sort(D,index);//归并从小到大排序
    if(sigma<0)
        reverse(index.begin(),index.end());
    vector<double> b(n);
    for(int i=0;i<n;i++){
        b[i]=z[index[i]];
        d[i]=D[index[i]];
    }
    vector<double> lambda(n);
    for(int i=0;i<n;i++){
        vector<double> delta(d.size());
        for(int j=0;j<delta.size();j++){
		    delta[j]=(d[j]-d[i])/sigma;
			delta[j]= prezero(delta[j]);//去0
		}
        double gamma=0;
        if(i+1<n){
            //gamma>1/delta[i+1]
            double A=b[i]*b[i];
            double B=-A/delta[i+1]-1;
			A = prezero(A);
            for(int j=0;j<delta.size();j++)
                if(j!=i)
                    B-=b[j]*b[j]/delta[j];
            double C=1.0;
            for(int j=0;j<delta.size();j++)
                if(j!=i)
                    C+=b[j]*b[j]/delta[j];
            C/=delta[i+1];
            C-=b[i+1]*b[i+1]/delta[i+1];
            gamma=(-B+sqrt(B*B-4*A*C))/(2*A);
			gamma = prezero(gamma);
        }
			
        //牛顿法迭代求解
        double diff=1;
        while(diff*diff>EPS){
            double g=0;
            for(int j=0;j<n;j++){
				double dg = delta[j]*gamma-1;
				dg = prezero(dg);
                g-=b[j]*b[j]/(dg*dg);
            }
			g = prezero(g); 
            double f=1;
            for(int j=0;j<n;j++){
				double idg = delta[j]-1/gamma;
				idg = prezero(idg);
                f+=b[j]*b[j]/(idg);
            }
            //f+g(newGamma-gamma)=0
            double newGamma=-f/g+gamma;
            diff=fabs(newGamma-gamma);
            gamma=newGamma;
        }
        lambda[i]=double(1)/gamma*sigma+d[i];
    }

    for(int i=0;i<n;i++)
        res[index[i]]=lambda[i];
    return res;
}

//分治过程
void DCSub(vector<double> &alpha, vector<double> &beta, vector<vector<double> > &Q, vector<double> &D, int start, int end){
    if(start==end){
        Q[start][start]=1;
        D[start]=alpha[start];
        return;
    }else{
        int mid=(start+end)/2;  //划分
        alpha[mid]-=beta[mid+1];  //统一协调秩1修补矩阵
        alpha[mid+1]-=beta[mid+1];
		DCSub(alpha,beta,Q,D,start,mid);  //递归
        DCSub(alpha,beta,Q,D,mid+1,end);
		//_Cilk_sync;
        int n=end-start+1;//矩阵规模
        vector<double> z(n,0);
        for(int i=start;i<=mid;i++)  //构造向量z=(q1',q2')
            z[i-start]=Q[mid][i];	//子矩阵最后一行
        for(int i=mid+1;i<=end;i++)
            z[i-start]=Q[mid+1][i];	//子矩阵第一行

        //计算矩阵 D+beta[mid+1]*z*z'的特征值
        vector<double> d(n,0);
        for(int i=0;i<n;i++)
            d[i]=D[i+start];//获得子矩阵特征值
		
		cout<<start<<" : "<<end<<endl;
		cout<<"z[],d[] completed ."<<endl;		
		
        //计算特征方程 1 + sum_j \frac{z^2_j}{d_j-\lambda} =0 中lambda的值
        vector<double> lambda = secularEquationSolver(z, d, beta[mid+1],start,end);//整合子矩阵特征值和修正矩阵特征值
		cout<<"lambda completed ."<<endl;	
        //对块内每个特征值计算局部特征向量 p = (D-\lambda I)^{-1} *z
        vector<vector<double> > P(n,vector<double>(n));
        for(int i=0;i<n;i++){
            vector<double> p(n);
            for(int j=0;j<n;j++){
				double tem = D[j+start]-lambda[i];
				p[j]= z[j]*double(1)/prezero(tem);
				/* if(_isnan(p[j])){//溢出或非确定性操作检测
					//cout<<j<<" ****IND***** "<<lambda[i];
				} */
				/*  if(start==250 && end == 499&& i>240){
					cout<<D[j+start]<<" "<<lambda[i]<<" "<<z[j]<<" "<<p[j]<<endl;
				}  */
			}
            normalize(p);
            for(int j=0;j<n;j++){
				P[j][i]=p[j];
			}
        }
        cout<<"P[][] completed."<<endl;
        vector<vector<double> > oldQ(n,vector<double>(n));
        for(int i=0;i<n;i++)
            for(int j=0;j<n;j++){
                oldQ[i][j]=Q[i+start][j+start];
            }
        
        for(int i=0;i<n;i++){	//更新当前Q矩阵
            for(int j=0;j<n;j++){
                Q[i+start][j+start]=0;
                for(int k=0;k<n;k++){
                    Q[i+start][j+start]+=oldQ[i][k]*P[k][j];
                }
            }
        }
		cout<<"Q updating completed."<<endl;
        //更新特征值
        for(int i=0;i<n;i++){
			D[i+start]=lambda[i];
		}	
		cout<<endl;
    }
}

//分治三对角入口
void DCTridiagonal(vector<double> alpha, vector<double> &beta, vector<vector<double> > &Q, vector<double> &D){
    int m=alpha.size();
    Q.clear();
    Q.resize(m,vector<double>(m,0));
    D.clear();
    D.resize(m,0);
    DCSub(alpha, beta, Q, D, 0, m-1);
}

//主程序入口
void resolve(vector<vector<double> > &A, int r, vector<vector<double> > &U, vector<double> &s, vector<vector<double> > &V){
    int m=A.size();
    int n=A[0].size();
    //lanczos: A*A'=P*T*P'
    if(m<=n){
        int l=m;
        vector<vector<double> > P(m,vector<double>(l,0));
        vector<double> alpha(l,0); 
        vector<double> beta(l,0);
        cout<<"lanczos..."<<endl;
        lanczos(A,P,alpha,beta,l);
        cout<<"lanczos completed. "<<endl;
        vector<vector<double> > W;
		vector<double> D(l,0);//特征值向量
		vector<vector<double> > Q;//特征向量矩阵
		DCTridiagonal(alpha,beta,Q,D);//调用分治法分解
		/*cout<<"Q       :"<<endl;
		 for(int i=0;i<m;i++){
			for(int j=0;j<m;j++){
				cout<<Q[i][j]<<" ";
			}
			cout<<endl;
		} */
		vector<int> index;	//排序后索引
		cout<<"sorting..."<<endl;
		merge_sort(D,index);
		cout<<"sorting completed."<<endl;
		reverse(index.begin(),index.end());	//逆序
		W.resize(l,vector<double>(l));
		for(int i=0;i<l;i++)
			for(int j=0;j<l;j++)
				W[i][j]=Q[i][index[j]];	//改变原Q的顺序赋给W
       
        U.clear();
		U.resize(m,vector<double>(l));
		cout<<"calculating U..."<<endl;
        multiply(P,W,U);//U=P*W
		cout<<"calculating U completed ."<<endl;
        for(int i=0;i<U.size();i++)
            U[i].resize(r);
        V.clear();V.resize(n,vector<double>(r));
		cout<<"calculating V..."<<endl;
        rightMultiply(U,A,V);//V=U'*A
		cout<<"calculating V completed ."<<endl;
        s.clear();s.resize(r,0);
        for(int i=0;i<r;i++){	//归一化
            for(int j=0;j<n;j++){
                s[i]+=V[j][i]*V[j][i];
            }
            s[i]=sqrt(s[i]);
            if(s[i]>EPS){
                for(int j=0;j<n;j++){
                    V[j][i] /= s[i];
                }
            }
        }
    }
}

void lanczos(vector<vector<double> > &A, vector<vector<double> > &P, vector<double> &alpha, vector<double> &beta, unsigned int rank){
    //P'*A*A'*P = T = diag(alpha) + diag(beta,1) + diag(beta, -1)
    //P=[p1,p2, ... , pk]
    vector<double> p;
    unsigned int m=A.size();
    unsigned int n=A[0].size();
    rank=min(n,min(m,rank));
    vector<double> prevP(m,0);
    randUnitVector(m,p);	//生成随机归一化P向量
    P.clear();
    P.resize(m,vector<double>(rank,0));
    vector<double> v;
    alpha.clear();
	alpha.resize(rank);
    beta.clear();
	beta.resize(rank);
    beta[0]=0;
    cout<<"initial completed..."<<endl;
    for(int i=0;i<rank;i++){
        for(int j=0;j<p.size();j++){
            P[j][i]=p[j];	//P的每一列都为p
        }
        rightMultiply(A, p, v);	//v=A'*P
        alpha[i]=dotProduct(p,v);
        if(i+1<rank){
            for(int j=0;j<m;j++)
                v[j]=v[j]-beta[i]*prevP[j]-alpha[i]*p[j];
            beta[i+1]=norm(v);
			beta[i+1] = prezero(beta[i+1]);
            prevP=p;
            for(int j=0;j<m;j++)
                p[j]=v[j]/beta[i+1];
        }
    }
}
