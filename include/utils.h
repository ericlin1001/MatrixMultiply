#ifdef UTILS_H
#else
#define UTILS_H
#include<assert.h>
#include<math.h>
#include<sys/time.h>
#include<unistd.h>
#include <fstream>
#include<iostream>
using namespace std;
#define MAX_BUFFER 100
#define MATH_PI M_PI
#define MATH_EXP M_E
#define Trace(m) {cout<<#m"="<<(m)<<endl;}
#define ASSERT(cond) if(!(cond)){cerr<<"Error: condition("#cond") fails!!"<<endl;};
#define Test(m) cout<<#m"={"; m; cout<<"}"<<endl;
void Tagg(const char *str){
	static int count=0;
	count++;
	printf("******Tag%d:%s\n",count,str);
}
double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;
	if ( phase == 0 ) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2.0 * U1 - 1.0;
			V2 = 2.0 * U2 - 1.0;
			S = V1 * V1 + V2 * V2;
		} while(S >= 1 || S == 0);
		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);
	phase = 1 - phase;
	return X;
}
inline double NormD(double u,double t){
	return gaussrand()*t+u;
}
double drand(){
	//[0,1);
	double r=rand();
	r/=((double)RAND_MAX+1);
	return r;
}
double drand(double min,double max){
	return drand()*(max-min)+min;
}

template<class T>
void printVec(const vector<T>&arr){
	cout<<"(";
	for(int i=0;i<arr.size();i++){
		if(i!=0)cout<<',';
		cout<<arr[i];
	}
	cout<<")";
}
class Tic{
	//accuration in milliseconds
	private:
		static long lastTime;
		Tic(){}
		inline static long getTimeMs(){
			timeval timeStart;
			gettimeofday(&timeStart,NULL);
			long res=((long)timeStart.tv_sec)*1000+(long)timeStart.tv_usec/1000;
			return res;
		}
	public:
		static long mtic(){
			//in milliseconds.
			long currentTime=getTimeMs();
			long dur=currentTime-lastTime;
			lastTime=currentTime;
			return dur;
		}
		static double tic(const char *tag="begin"){
			double res;
			res=dtic();
			if(strcmp(tag,"begin")==0){
				cout<<"Tic::"<<tag<<endl;
			}else{
				cout<<"Tic::"<<tag<<" used:"<<res<<"(seconds)."<<endl;
			}
			return res;
		}
		inline static double dtic(){
			//in seconds.
			return (double)mtic()/1000.0;
		}
		static void test(){
			Tic::mtic();
			usleep(1234000);//sleep for 1234 milliseconds.(1.234seconds)
			cout<<Tic::dtic()<<"seconds"<<endl;
		}

};
long Tic::lastTime=0;

void printArr(int *arr,int size){
	cout<<"(";
	for(int i=0;i<size;i++){
		if(i!=0)cout<<',';
		cout<<arr[i];
	}
	cout<<")";
}
class RandomPermutation{
	vector<int>p;
	int i;
	int n;
	public:
		RandomPermutation(int tn):n(tn),i(0){
			p.resize(n);
			for(int i=0;i<n;i++){
				p[i]=i;
			}
		}
		void generate(){
			i=0;
		}
		int next(){
			//at most invoked n times before next generate.
			assert(i<n&&i>=0);
			int r=rand()%(n-i);
			i++;
			swap(p[n-i],p[r]);
			return p[n-i];
		}
};
void calStatistics(const vector<double>&arr,double &min,double &max,double &mean,double &std){
	min=arr[0];
	max=arr[0];
	mean=0.0;
	for(int i=0;i<arr.size();i++){
		double x=arr[i];
		if(x<min)min=x;
		if(x>max)max=x;
		mean+=x;
	}
	mean/=(double)arr.size();
	std=0.0;
	for(int i=0;i<arr.size();i++){
		double x=arr[i];
		std+=pow(x-mean,2);
	}
	std/=(double)arr.size();
	std=sqrt(std);
}
class Save{
	private:
			ofstream out;
			int ix;
	public:
			Save(){
			}
		Save(const char *title){
			init(title);
		}
		Save(const char *title,const char *x,const char *y){
			init(title);
			setXY(x,y);
		}
		void setXY(const char *x,const char *y){
			out<<x<<"\t"<<y<<endl;
			ix=0;
		}
		void init(const char *title){
			char buff[150];
			sprintf(buff,"%s",title);
			out.open(buff);
			if(!out.is_open()){
				cerr<<"Error:Can't open file:"<<buff<<endl;
			}
		}
		void add(int x,double y){
			out<<x<<"\t"<<y<<endl;
		}
		void add(double y){
			//start from 0
			out<<ix<<"\t"<<y<<endl;
			ix+=1;
		}
		~Save(){
			out.close();
		}
};
#endif
