#include "include/template.h"
#include "include/BasicDE.h"
#include "include/SignalHandleHelper.h"
#include "include/IDHelper.h"
#include "pecFunction.h"
#include<iostream>
#include<fstream>
#undef DEBUG
//#define DEBUG
using namespace std;
DefFunction(TestF,-100,100,100)
	return xs[0];
	EndDef

	class ParallelDE:public EA
{
	private:
		//about function:f
		Function *f;
		vector<vector<double> >range;
		int numDim;
		//algorithm related parameters.
		int PopSize;
		double F,CR;
		//
		vector<vector<double> >x;//x,trail x.
		vector<vector<double> >tmpX;
		vector<double>fx;
		vector<double>tmpFx;
		vector<double>tx;
		//
		int bestI;
		MPIHelper*mpi;
		SearchParam *param;
		Save *save;
	private:
		inline void updateX(){
			updateX_schema2();
		}
		void updateX_schema3(){
			//main process
			vector<vector<double> >txs;
			vector<double> ftxs;
			txs.resize(PopSize);

			RandomPermutation perm(PopSize);
			for(int i=0;i<PopSize;i++){
				perm.generate();
				int a=perm.next(); int b=perm.next(); int c=perm.next();
				if(a==i){a=perm.next();}
				if(b==i){b=perm.next();}
				if(c==i){c=perm.next();}
				int randDim=rand()%numDim;
				for(int j=0;j<numDim;j++){
					if(j==randDim||drand()<=CR){
						tx[j]=x[a][j]+F*(x[b][j]-x[c][j]);
						if(tx[j]<range[j][0] || tx[j]>range[j][1]){
							tx[j]=drand(range[j][0],range[j][1]);
						}
					}else{
						tx[j]=x[i][j];
					}
				}
				txs[i]=tx;
			}
			evaluatePopulation(txs,ftxs);
			for(int i=0;i<PopSize;i++){
				vector<double>&tx=txs[i];
				double &ftx=ftxs[i];
				if(ftx<fx[i]){
					x[i]=tx;
					fx[i]=ftx;
					if(ftx<fx[bestI]){
						bestI=i;
					}
				}
			}
		}
		void updateX_schema2(){
			//main process
			vector<vector<double> >txs;
			vector<double> ftxs;
			txs.resize(PopSize);

			RandomPermutation perm(PopSize);
			for(int i=0;i<PopSize;i++){
				perm.generate();
				int a=perm.next(); int b=perm.next(); int c=perm.next();
				if(a==i){a=perm.next();}
				if(b==i){b=perm.next();}
				if(c==i){c=perm.next();}
				int randDim=rand()%numDim;
				for(int j=0;j<numDim;j++){
					if(j==randDim||drand()<=CR){
						tx[j]=x[a][j]+F*(x[b][j]-x[c][j]);
						if(tx[j]<range[j][0]){
							//	tx[j]=drand(range[j][0],range[j][1]);
							tx[j]=range[j][0];
						}
						if(tx[j]>range[j][1]){
							tx[j]=range[j][1];
						}
					}else{
						tx[j]=x[i][j];
					}
				}
				txs[i]=tx;
			}
			evaluatePopulation(txs,ftxs);
			for(int i=0;i<PopSize;i++){
				vector<double>&tx=txs[i];
				double &ftx=ftxs[i];
				if(ftx<fx[i]){
					x[i]=tx;
					fx[i]=ftx;
					if(ftx<fx[bestI]){
						bestI=i;
					}
				}
			}
		}
		void schema1_updateX(){
			//main process
			vector<vector<double> >txs;
			vector<double> ftxs;
			txs.resize(PopSize);

			RandomPermutation perm(PopSize);
			for(int i=0;i<PopSize;i++){
				perm.generate();
				int a=bestI; int b=perm.next(); int c=perm.next();
				int randDim=rand()%numDim;
				for(int j=0;j<numDim;j++){
					if(j==randDim||drand()<=CR){
						tx[j]=x[a][j]+F*(x[b][j]-x[c][j]);
						if(tx[j]<range[j][0] || tx[j]>range[j][1]){
							tx[j]=drand(range[j][0],range[j][1]);
						}
					}else{
						tx[j]=x[i][j];
					}
				}
				txs[i]=tx;
			}
			evaluatePopulation(txs,ftxs);
			for(int i=0;i<PopSize;i++){
				vector<double>&tx=txs[i];
				double &ftx=ftxs[i];
				if(ftx<fx[i]){
					x[i]=tx;
					fx[i]=ftx;
					if(ftx<fx[bestI]){
						bestI=i;
					}
				}
			}
		}
	private:
		double getBestFx()const{
			return fx[bestI];
		}
		void update(int maxGeneration){
#define SaveData if(save!=NULL){save->add(getBestFx());}
			SaveData;
			//cout<<"0:F(g="<<0<<")="<<getBestFx()<<endl;
			cout<<"0:"<<getBestFx()<<endl;
			for(int g=1;g<=maxGeneration;g++){
				//	if(maxGeneration<30||g%(maxGeneration/30)==0){
				//		cout<<g/(maxGeneration/30)<<":F(g="<<g<<")="<<getBestFx()<<endl;
				//	}
				cout<<g<<":"<<getBestFx()<<endl;
				updateX();
				SaveData;
			}
		}
	public:
		ParallelDE(){}
		ParallelDE(MPIHelper *h):mpi(h),save(0){
		}
		ParallelDE(MPIHelper *h,SearchParam*p):mpi(h),param(p){
			initParam(param);
		}
		~ParallelDE(){
			//	endEvaluate();
		}
		void setParam(MPIHelper *h,SearchParam*p){
			mpi=h;
			param=p;
			initParam(param);
		}
		void addSave(Save *s){
			save=s;
		}
		void initParam(SearchParam *param){
			this->param=param;
			//
			PopSize=param->getInt("PopSize");
			F=param->getDouble("F");
			CR=param->getDouble("CR");
			param->getBiVector("Range",range);
#ifdef DEBUG
			cout<<"Range:";
			for(int i=0;i<range.size();i++){
				printf("[%g,%g],",range[i][0],range[i][1]);
			}
			cout<<endl;
#endif
			setName(param->getString("Name"));
		}

		void calulateBestI(){
			bestI=0;
			for(int i=0;i<PopSize;i++){
				if(fx[i]<fx[bestI]){ 
					bestI=i;
				}
			}
		}
#define MESS_END 0
#define MESS_EVAL_ARRAY 2

		void generateSplitTask(int numTask,int numProcesses,vector<int>&task){
			//uniformly distribute the tasks among all processes.
			task.resize(numProcesses+1);
			int numXPerProcesses=numTask/numProcesses;
			int numRemain=numTask-numXPerProcesses*numProcesses;
			task[0]=0;
			int i;
			for(i=0;i<numRemain;i++){
				task[i+1]=task[i]+numXPerProcesses+1;
			}
			for(;i<numProcesses;i++){
				task[i+1]=task[i]+numXPerProcesses;
			}
		}
		void evaluatePopulation(vector<vector<double> >&xs,vector<double>&fx){
			fx.resize(xs.size());
			int numSlaves=mpi->getNumProcesses()-1;
			ASSERT(numSlaves>0);
			//
			vector<int>task;
			generateSplitTask(xs.size(),numSlaves,task);

			for(int i=1;i<=numSlaves;i++){
				int dest=i;
				int len=task[i]-task[i-1];
				if(len>=1){
					mpi->send(MESS_EVAL_ARRAY,dest);
					mpi->send(len,dest);//len
					//cout<<"client("<<dest<<") "<<"end send len:"<<len<<endl;
					for(int j=task[i-1];j<task[i];j++){
						mpi->send(&xs[j][0],xs[j].size(),dest);
					}
					//cout<<"client("<<dest<<") "<<"end send ."<<endl;
				}
			}

			for(int i=1;i<=numSlaves;i++){
				int dest=i;
				//cout<<"client("<<dest<<") "<<"recv from ."<<endl;
				for(int j=task[i-1];j<task[i];j++){
					mpi->recv(&fx[j],1,dest);
				}
				//cout<<"client("<<dest<<") "<<"end recv ."<<endl;
			}
		}
		void endEvaluate(){
			for(int i=1;i<mpi->getNumProcesses();i++){
				mpi->send(MESS_END,i);
			}
		}
		virtual double getMin(Function *f,int MaxFEs,vector<double>&out_x,double &out_fx){
			if(save!=0){
				save->setXY("Generation",f->getName());
			}
			if(mpi->getNumProcesses()<=1){
				BasicDE de;
				de.initParam(param);
				cout<<"Warning:NumProcesses<=1, Use BasicDE()"<<endl;
				return de.getMin(f,MaxFEs,out_x,out_fx);
			}
			//allocating space.
			this->f=f;
			numDim=f->getNumDim();
			tx.resize(numDim);
			x.resize(PopSize);
			fx.resize(PopSize);
			for(int i=0;i<PopSize;i++){
				x[i].resize(numDim);
			}
			//
			if(mpi->isMaster()){
				ASSERT(range.size()>=f->getNumDim());
				ASSERT(range[0].size()>=2);
				//population initializing....
				for(int i=0;i<PopSize;i++){
					for(int d=0;d<numDim;d++){
						x[i][d]=drand(range[d][0],range[d][1]);
					} 
				}
				evaluatePopulation(x,fx);
				calulateBestI();
				cout<<"*******************Start update per generation"<<endl;
				//update, main process.
				update(MaxFEs/PopSize-1);
				cout<<"*******************end Start update per generation"<<endl;
				endEvaluate();//stop evaluating....
				calulateBestI();
				out_x=x[bestI];
				out_fx=fx[bestI];
				return out_fx;
			}else{//slavery processes,only evaluate the f(x).
				bool isEnd=false;
#ifdef DEBUG
				cout<<mpi->getName()<<" starts..."<<endl;
#endif
				while(!isEnd){
					int type;
					int len;
#ifdef DEBUG
					cout<<mpi->getName()<<" is waiting...."<<endl;
#endif
					mpi->recv(type,0);
					switch(type){
						case MESS_END:
#ifdef DEBUG
							cout<<mpi->getName()<<" is exiting..."<<endl;
#endif
							isEnd=true;
							break;
						case MESS_EVAL_ARRAY:
#ifdef DEBUG
							cout<<mpi->getName()<<":MESS_EVAL_ARRAY"<<endl;
#endif
							mpi->recv(len,0);
#ifdef DEBUG
							cout<<mpi->getName()<<" recv len:"<<len<<endl;
#endif
							for(int i=0;i<len;i++){
								mpi->recv(&x[i][0],numDim,0);
							}
#ifdef DEBUG
							cout<<mpi->getName()<<" end recv"<<len<<endl;
#endif
							for(int i=0;i<len;i++){
								fx[i]=f->evaluate(&x[i][0]);
							}
#ifdef DEBUG
							cout<<mpi->getName()<<" end evaluation"<<len<<endl;
#endif
							for(int i=0;i<len;i++){
								mpi->send(&fx[i],1,0);
							}
#ifdef DEBUG
							cout<<mpi->getName()<<" end send"<<len<<endl;
#endif
							break;
						default:
							cerr<<"Error:Unknown mess_type"<<endl;
							//assert(false);
							break;
					}
				}
				return -1;
			}
		}
};
///////////////////////
void saveConfigData(int id,const char *f,const char *algorithm,const char *param,int run,int MaxRun,int numOfProcesses,int MaxFEs,int PopSize,int NumDim,double F,double CR,const char *state,double usedTime,double absError,
		vector<double>&x,double fx){
	ofstream runConfig;
	char buff[1000];
	sprintf(buff,"Run-configuration-%d.txt",id);
	if(strcmp(state,"start")==0){
		cout<<"Save file:"<<buff<<endl;
	}
	runConfig.open(buff,ofstream::app|ofstream::out);
	runConfig<<"ID\tFunction\tAlgorithm\tParamemterFile\tRun\tMaxRun\tNumOfProcesses\tMaxFEs\tPopSize\tNumDim\tF-parameter\tCR-paramter\tState\tUsedTime\tAbsError\tX\tFx"<<endl;
	sprintf(buff,"%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%g\t%s\t%g\t%g",id,f,algorithm,param,run,MaxRun,numOfProcesses,MaxFEs,
			PopSize,NumDim,F,CR,state,usedTime,absError);
	runConfig<<buff<<"\t";
	if(x.size()==0){
		runConfig<<"-1";
	}else{
		for(int i=0;i<x.size();i++){
			sprintf(buff,"%g",x[i]);
			if(i!=0){
				runConfig<<",";
			}
			runConfig<<buff;
		}
	}
	sprintf(buff,"\t%g",fx);
	runConfig<<buff;

	runConfig<<endl;
	runConfig.flush();
	runConfig.close();
}

ParallelDE *de=0;
void IntHandler(int s){
	if(de!=NULL){
		delete de;
		de=NULL;
	}
	cout<<"My:Caught signal SIGINT"<<endl;
	exit(1);
}
int MainProgram(MPIHelper*mpi,int run,int MaxRun,int configID,ParallelDE*de,Function*f,SearchParam*param,bool isFindMin=true){
	/*************Shared Data*************/
	vector<double>x;
	double fx=-1;
	srand(time(NULL));
	Save save;
	/**************end Shared Data***********/

	/********Master data*************/
	char state[50]="start";
	double usedTime;
	double absError;
	/**********end Master data*************/

	//
	//
	if(mpi->isMaster()){
		cout<<"Runing "<<de->getName()<<" "<<endl;
		cout<<"FunName(MyBestF,Optima)"<<endl;
		//
		usedTime=-1;
		absError=-1;
		//
		char saveFileName[100];
		sprintf(saveFileName,"Data-%d.txt",configID);
		save.init(saveFileName);
		de->addSave(&save);
		saveConfigData(configID,f->getName(),de->getName(),param->getName(),run,MaxRun,mpi->getNumProcesses(),param->getInt("MaxFEs"),
				param->getInt("PopSize"),param->getInt("NumDim"),param->getDouble("F"),param->getDouble("CR"),state,usedTime,absError,x,fx);
		Tic::tic("begin");
	}
	cout<<mpi->getName()<<" starts computing..."<<endl;
	if(isFindMin){
		de->getMin(f,param->getInt("MaxFEs"),x,fx);
	}else{
		de->getMax(f,param->getInt("MaxFEs"),x,fx);
	}
	if(mpi->isMaster()){
		usedTime=Tic::tic("end");
		absError=fabs(fx-f->getFBest());
		strcpy(state,"end");

		saveConfigData(configID,f->getName(),de->getName(),param->getName(),run,MaxRun,mpi->getNumProcesses(),param->getInt("MaxFEs"),
				param->getInt("PopSize"),param->getInt("NumDim"),param->getDouble("F"),param->getDouble("CR"),state,usedTime,absError,x,fx);

		printf("%s(%g,%g)\n",f->getName(),fx,f->getFBest());
		cout<<"ends successfully!"<<endl;
	}
	return 0;
}
MPIHelper *g_mpi=0;
class Matrix{
	int numRows,numCols;
	vector<vector<double > >values;
	public:
	Matrix(int numRow,int numCol){
		this->numRows=numRow;
		this->numCols=numCol;
		values.resize(numRows);
		for(int r=0;r<numRows;r++){
			values[r].resize(numCols);
		}
	}
	void set(int r,int c,double v){
		values[r][c]=v;
	}
	double get(int r,int c)const{
		return values[r][c];
	}
	int getNumRows()const{return numRows;}
	int getNumCols()const{return numCols;}
	vector<double> getRow(int r)const{
		return values[r];
	}
	vector<double> getCol(int c)const{
		vector<double>col;
		col.resize(numRows);
		for(int r=0;r<numRows;r++){
			col[r]=get(r,c);
		}
		return col;
	}
	Matrix multiply(const Matrix&right){
		assert(this->numCols==right.numRows);
		Matrix v(numRows,right.numCols);
		for(int r=0;r<v.getNumRows();r++){
			for(int c=0;c<v.getNumCols();c++){
				double value=0;
				for(int k=0;k<this->numCols;k++){
					value+=this->get(r,k)*right.get(k,c);
				}
				v.set(r,c,value);
			}
		}
		return v;
	}
	Matrix pmultiply(const Matrix&right){
		//for A(m*n)*B(n*p) parallel for n.(speedup by n).
		assert(this->numCols==right.numRows);
		Matrix v(numRows,right.numCols);
		for(int r=0;r<v.getNumRows();r++){
			for(int c=0;c<v.getNumCols();c++){
				int index=r*v.getNumCols()+c;
				int assignID=index%(g_mpi->getNumProcesses());
				if(assignID==g_mpi->getID()){
					double value=0;
					for(int k=0;k<this->numCols;k++){
						value+=this->get(r,k)*right.get(k,c);
					}
					v.set(r,c,value);
				}
			}
		}
		for(int r=0;r<v.getNumRows();r++){
			for(int c=0;c<v.getNumCols();c++){
				int index=r*v.getNumCols()+c;
				int assignID=index%(g_mpi->getNumProcesses());
				if(g_mpi->isMaster()){
					if(assignID!=0){
						double value;
						g_mpi->recv(value,assignID);
						v.set(r,c,value);
					}
				}else{
					if(assignID==g_mpi->getID()){
						g_mpi->send(v.get(r,c),0);
					}
				}
			}
		}
		return v;
	}
	void print(){
		cout<<"["<<endl;
		for(int r=0;r<getNumRows();r++){
			for(int c=0;c<getNumCols();c++){
				cout<<get(r,c)<<"\t";
			}
			cout<<endl;
		}
		cout<<"]"<<endl;
	}
};

Matrix initMatrix(int row,int cols){
	Matrix a(row,cols);
	for(int r=0;r<a.getNumRows();r++){
		for(int c=0;c<a.getNumCols();c++){
			a.set(r,c,r*a.getNumCols()+c);
		}
	}
	return a;
}
void serialTest(bool debug=true){
	if(g_mpi->isMaster()){
		//const int NumRows=2,NumCols=3;
		Matrix a=initMatrix(1,3);
		Matrix b=initMatrix(3,1);
		if(debug){
		cout<<"A=";
		a.print();
		cout<<"B=";
		b.print();
		cout<<"A*B=";
		}
		Matrix c=b.multiply(a);
		if(debug){
		c.print();
		cout<<endl;
		}
	}
}
void parallelTest(bool debug=true){
	Matrix a=initMatrix(1,3);
	Matrix b=initMatrix(3,1);
		if(debug){
	if(g_mpi->isMaster()){
		cout<<"A=";
		a.print();
		cout<<"B=";
		b.print();
		cout<<"A*B=";
		}
	}
	Matrix c=b.pmultiply(a);
		if(debug){
	if(g_mpi->isMaster()){
		c.print();
		cout<<endl;
		}
	}
}
int main(int argc,char *argv[]){
	const int Times=10000;
	MPIHelper mpi(argc,argv);

	g_mpi=&mpi;
	if(g_mpi->isMaster())
	Tic::tic("begin");
	for(int i=0;i<Times;i++)
	serialTest(false);
	if(g_mpi->isMaster())
	Tic::tic("serial");
	//
	for(int i=0;i<Times;i++)
	parallelTest(false);
	if(g_mpi->isMaster())
	Tic::tic("parallel");
}
