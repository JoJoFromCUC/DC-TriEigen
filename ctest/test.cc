#include <iostream>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <unistd.h>

using namespace std;
int main(){
	vector<int> b(50);
	cout<< numeric_limits<double>::min()<<endl;
	/* FILE *fp;
	int a[100];
	fp = fopen("./data.txt","r");
	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j) {
			if (fscanf(fp, "%d", &a[i][j]) == EOF) break;
		}
		if (feof(fp)) break;
	}
	for (int i = 0; i < 10; ++i) {
		for (int j = 0; j < 10; ++j) {
			cout << a[i][j] << " ";
		}
		cout << endl;
	} */
	/* ifstream in("./data1.txt");
	if(!in) { 
        cerr<<"Can't open the file."<<endl; 
        return -1; 
    } 
	int temp;
	b.clear();
	while(in >> temp){
		b.push_back(temp);
	}
	for(int i=0;i<b.size();i++){
		cout<<b[i]<<" ";
	}
	cout<<endl;
	in.close();  */
	fstream ofile;
	ofile.open("./data5000.txt",ios::out|ios::trunc);
	for(int i=0;i<5000;i++){
		for(int j=0;j<5000;j++){
			if(j == i) ofile<<4<<" ";
			else if(j==i-1||j==i+1) ofile<<1<<" ";
			else ofile<<0<<" ";
		}
		ofile<<"\n";
	}
	ofile.close();
	clock_t start, finish;
	start = clock();
	sleep(3);//<unistd.h>
    //cout<<time(NULL);
    //cout<<srand(time(NULL));
	finish = clock();
	cout << (double)(finish - start) / CLOCKS_PER_SEC;
    system("pause");
    return 0;
}