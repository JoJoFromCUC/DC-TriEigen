#include "dc.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>


//特征方程求解
vector<double> secularEquationSolver(vector<double> &z, vector<double> &D, double sigma,int start,int end){//sigma=beta[start+end/2]
    int n=z.size(),num;
    vector<double> res(n);
    vector<int> index;
    vector<double> d(n);
    vector<double> b(n);
    vector<double> lambda(n);
    
    merge_sort(D,index);   //排序
    if(sigma<0){
        reverse(index.begin(),index.end());
    }
    vector<double> new_b;
    vector<double> new_d;
    for(int i=0;i<n;i++){
        b[i]=z[index[i]];
        d[i]=D[index[i]];
        if(b[i]==0){  //直接得出特征值
           continue;
        } 
        new_b.push_back(b[i]); //获取不重复的d和不为0的b
        new_d.push_back(d[i]);
    }
    // for(int i=0;i<n;i++){
    //     cout<<new_b[i]<<" "<<new_d[i]<<endl;
    // }
    num = new_d.size();
    cout<<"need to computer number:"<<num<<endl;
    vector<double> tem_lambda(num); //除去相同或相似元素后的局部特征向量

    for(int i=0;i<num;i++){
        vector<double> delta(num);
        for(int j=0;j<num;j++){
		    delta[j]=(new_d[j]-new_d[i])/sigma;//d[j],d[i]太过接近时需要处理
            //if(delta[j] <ZERO  && j!=i) cout<<"i j"<<i<<j<<" "<<new_d[i]<<" "<<new_d[j]<<endl; //delta[j] = prezero(delta[j]);
            delta[j] = prezero(delta[j]);//加入扰动
		}
        double gamma=0;
        if(i+1<num){
            //gamma>1/delta[i+1]
            double A = new_b[i]*new_b[i];//A过小导致gamma无穷大
            double B = -A/delta[i+1]-1;
            //if(A<EPS) A = ZERO;//cout<<"++++++++++++++++++++"<<i<<"++++++++++++++++++++"<<endl;
			A = prezero(A);
            for(int j=0;j<delta.size();j++){
                if(j!=i){//避免delta[j]=0
                    B -= new_b[j]*new_b[j]/delta[j];
                }
            }
            double C=1;
            for(int j=0;j<delta.size();j++){
                if(j!=i){
                    C += new_b[j]*new_b[j]/delta[j];
                }
            }
            C /= delta[i+1];
            C -= new_b[i+1]*new_b[i+1]/delta[i+1];
            gamma = (-B+sqrt(B*B-4*A*C))/(2*A);
            cout<<"gamma"<<i<<" :"<<gamma<<endl;
			gamma = prezero(gamma);
        }
        //牛顿法迭代求解
        double diff=1;
        //int count=0;
        //出现不收敛情况
        while(diff*diff>EPS){
            double g=0;
            for(int j=0;j<num;j++){
				double dg = delta[j]*gamma-1;
				//dg = prezero(dg);
                g -= new_b[j]*new_b[j]/(dg*dg);
            }
			g = prezero(g); 
            double f=1;
            for(int j=0;j<num;j++){
				double idg = delta[j]-1/gamma;
				//idg = prezero(idg);
                f += new_b[j]*new_b[j]/(idg);
            }
            //f+g(newGamma-gamma)=0
            double newGamma = -f/g + gamma;
            diff=fabs(newGamma-gamma);
            gamma=newGamma;
            //count++;
            //if(count > 40000) break; 
        }
        //cout<<"gamma"<<i<<" :"<<gamma<<endl;

        tem_lambda[i] = 1/gamma*sigma + new_d[i];
    }
    int pos = 0;
    for(int i=0;i<n;i++){
        if(b[i]==0){  //直接得出特征值
           lambda[i] = d[i];
           continue;
        } 
        lambda[i] = tem_lambda[pos];
        pos++;
    }
    for(int i=0;i<n;i++){
        res[index[i]] = lambda[i];
    }
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
        alpha[mid]-=beta[mid+1];  //统一协调秩1修正矩阵
        alpha[mid+1]-=beta[mid+1];
	    DCSub(alpha,beta,Q,D,start,mid);  //递归
        DCSub(alpha,beta,Q,D,mid+1,end);
	    //cilk_sync;
        int n=end-start+1;//矩阵规模
        vector<double> z(n,0);
        for(int i=start;i<=mid;i++)  //构造向量z=(q1',q2')
            z[i-start]=Q[mid][i];	//子矩阵最后一行
        for(int i=mid+1;i<=end;i++)
            z[i-start]=Q[mid+1][i];	//子矩阵第一行

        /*
          z中对0元素处理的处理  
        */    

        //计算矩阵 D+beta[mid+1]*z*z'的特征值
        vector<double> d(n,0);
        for(int i=0;i<n;i++){
            d[i]=D[i+start];//获得子矩阵特征值
        }

		/*
        d[i]出现相同或相近元素的处理
        执行givens变换 for z[i],z[j]
        */
        bool flag = false;
        vector<vector<double> > givens(n,vector<double>(n,0));//生成givens矩阵
        for(int i=0;i<n;i++){
            givens[i][i] = 1;
            for(int j=i+1;j<n;j++){
                if( fabs(d[i]-d[j]) < GAP){
                    flag = true;
                    givens[i][i] = givens[j][j] = givens[i][j] = 0.7071; //givens(i,j,PI/4);
                    givens[j][i] = -0.7071;
                    z[i] = sqrt(z[i]*z[i]+z[j]*z[j]);
                    z[j] = 0;
                } 
            }
        }
        cout<<start<<" : "<<end<<endl;
        cout<<"z[],d[] completed ."<<endl;		
		
        //计算特征方程 1 + sum_j \frac{z^2_j}{d_j-\lambda} =0 中lambda的值
        vector<double> lambda = secularEquationSolver(z, d, beta[mid+1],start,end);//求解合并后的特征值
		cout<<"lambda completed ."<<endl;	
        //对块内每个特征值计算局部特征向量 p = (D-\lambda I)^{-1} *z
        vector<vector<double> > P(n,vector<double>(n,0));//局部特征向量矩阵
        vector<double> p(n);
        for(int i=0;i<n;i++){
            if(z[i] == 0){//z[i]==0,则直接得出特征向量等于ei
                // for(int j=0;j<n;j++){
                //     if(j==i) P[j][i] = 1;
                //     else P[j][i] = 0;
                // }
                P[i][i] = 1;
                continue;
            }

            /*
            update z[]

            求解特征向量存在问题,需要使用新的公式先更新z的值
            */
            for(int j=0;j<n;j++){ 
                double tem = D[j+start]-lambda[i];
                p[j]= z[j]/prezero(tem);
                //cout<<D[j+start]<<" "<<lambda[i]<<" "<<z[j]<<" "<<p[j]<<endl;    
            }
            normalize(p);        //归一化
            for(int j=0;j<n;j++){
                P[j][i]=p[j];
            }
        }
        cout<<"P[][] completed."<<endl;
        vector<vector<double> > oldQ(n,vector<double>(n));
        for(int i=0;i<n;i++){
            for(int j=0;j<n;j++){
                oldQ[i][j]=Q[i+start][j+start];
            }
        }  
        if(flag){
            vector<vector<double> > tempQ(n,vector<double>(n));
            for(int i=0;i<n;i++){	//multiply givens;
                for(int j=0;j<n;j++){
                    //temp[i+start][j+start]=0;
                    for(int k=0;k<n;k++){
                        tempQ[i][j]+=oldQ[i][k]*givens[k][j];
                    }
                }
            }
            for(int i=0;i<n;i++){
                for(int j=0;j<n;j++){
                    oldQ[i][j]=tempQ[i][j];
                }
            }  
        }     
        for(int i=0;i<n;i++){	//update Q
            for(int j=0;j<n;j++){
                Q[i+start][j+start]=0;
                for(int k=0;k<n;k++){
                    Q[i+start][j+start]+=oldQ[i][k]*P[k][j];
                }
            }
        }
		cout<<"Q updating completed."<<endl;
        //Update D
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
void resolve(SparseMatrix &A, int r,vector<double> &alpha,vector<double> &beta){
    int m=A.rows;
    int n=A.cols;
    //lanczos: A*A'=P*T*P'
    if(m==n){
        int l=m;
        vector<vector<double> > P(m,vector<double>(r,0));
        vector<vector<double> > U,V;
	    vector<double> s;
        
        //lanczos(A,P,alpha,beta,l);//三对角化
        vector<vector<double> > W;
        vector<double> D(l,0);
        vector<vector<double> > Q(l,vector<double>(l,0));

        DCTridiagonal(alpha,beta,Q,D);//调用分治法分解
        
        vector<int> index;	//排序后索引
        cout<<"sorting..."<<endl;
        merge_sort(D,index);
        cout<<"sorting completed."<<endl;
        reverse(index.begin(),index.end());	//逆序
        W.resize(l,vector<double>(l));
        for(int i=0;i<l;i++){
            for(int j=0;j<l;j++){
                W[i][j]=Q[i][index[j]];	//改变原Q的顺序赋给W
            }
        }
        //print(W);
		// for(int i=0;i<m;i++){
		// 	cout<<W[0][i]<<endl;
		// }
        //load file P
        FILE *fc ;
        cout<<"read standard matrix :"<<endl;
        fc = fopen("./ctest/P500.txt","r");
		if(!fc){
			cout<<"open P file failed!"<<endl;
		}
		for(int i=0;i<m;i++){
			double tmp;
            for(int j=0;j<r;j++){
                if(fscanf(fc,"%lf",&tmp) == EOF) break; 
                P[i][j] = tmp;               
            }
			if(feof(fc)) break;
		}
		fclose(fc);
        // print(P);
        U.clear();
        U.resize(m,vector<double>(l));
        cout<<"calculating U..."<<endl;
        multiply(P,W,U); //U=P*W 特征向量矩阵
        cout<<"calculating U completed ."<<endl;
        for(int i=0;i<U.size();i++){
            U[i].resize(r);
        }
        V.clear();V.resize(n,vector<double>(r));
        cout<<"calculating V..."<<endl;
        rightMultiply(U,A,V); //V=U'*A
        cout<<"get eigenvalues ready ..."<<endl;
        s.clear();s.resize(r,0);
        for(int i=0;i<r;i++){	//归一化
            for(int j=0;j<n;j++){
                s[i]+=V[j][i]*V[j][i];
            }
            s[i]=sqrt(s[i]);
            // if(s[i]>EPS){
            //     for(int j=0;j<n;j++){
            //         V[j][i] /= s[i];
            //     }
            // }
        }
        cout<<"eigen values :"<<endl;//输出特征值
        for(int i=0;i<s.size();i++){
            cout<<s[i]<<endl;
        }
        // FILE *fp = fopen("./result.txt","w");//写入特征向量矩阵U
        // int rows = U.size(),cols = U[0].size();
        // fprintf(fp,"%s","特征向量矩阵为:\n");
        // for(int i=0;i<rows;++i){
        //     for(int j=0;j<cols;++j){
        //         fprintf(fp,"%.15lf ",U[i][j]);
        //     }
        //     fprintf(fp,"%s","\n");
        // }
        // fclose(fp);
        // cout<<"特征向量矩阵已写入result.txt文件!\n"<<endl; 
        //print(U);
        cout<<endl;
        
    }
    else{
        cout<<"the input matrix is illegal !"<<endl;
    }
}