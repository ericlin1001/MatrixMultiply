//#define DEBUG
#include "include/template.h"
#include "include/BasicDE.h"
#include "include/SignalHandleHelper.h"
#include "include/IDHelper.h"
#include "pecFunction.h"
#include<iostream>
#include<fstream>
#undef DEBUG
using namespace std;
DefFunction(MyF1,-100,100,0)
	return xs[0];
	/*
	   for(int i=0;i<size;i++){
	   double x=xs[i];
	   res+=x*x;
	   }
	   return res;
	 */
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
	private:
		void updateX(){
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
					if(j==randDim||drand()<CR){
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
#define SaveData s.add(getBestFx());
			Save s(f->getName(),"Generation","F");
			SaveData;
			for(int g=1;g<=maxGeneration;g++){
				if(g%10==0){
					Trace(g);
					cout.flush();
				}
				updateX();
				SaveData;
				//Tic::tic("update one X");
			}
		}
	public:
		ParallelDE(MPIHelper *h):mpi(h){
		}
		~ParallelDE(){
			cout<<"ParallelDE end1."<<endl;
			endEvaluate();
			cout<<"ParallelDE end2."<<endl;
		}
		void initParam(SearchParam *param){
			PopSize=param->getInt("PopSize");
			F=param->getDouble("F");
			CR=param->getDouble("CR");
			param->getBiVector("Range",range);
			setName(param->getString("Name"));
		}

		void calulateBestI(){
			for(int i=0;i<PopSize;i++){
				if(fx[i]<fx[bestI]){ 
					bestI=i;
				}
			}
		}
#define MESS_END 0
#define MESS_EVAL_ARRAY 2
		void generateSplitTask(int numTask,int numProcesses,vector<int>&task){
			task.resize(numProcesses+1);
			int numXPerProcess=(numTask-1)/numProcesses+1;
			task[0]=0;
			for(int i=1;i<numProcesses+1;i++){
				task[i]=task[i-1]+numXPerProcess;
			}
			task[numProcesses]=numTask;
		}
		void evaluatePopulation(vector<vector<double> >&xs,vector<double>&fx){
			fx.resize(xs.size());
			int numSlaves=mpi->getNumProcesses()-1;
			ASSERT(numSlaves>0);
			//
			vector<int>task;
			//TODO: what happnens if numSlaves >PopSize;
			generateSplitTask(xs.size(),numSlaves,task);

			for(int i=1;i<=numSlaves;i++){
				int dest=i;
				mpi->send(MESS_EVAL_ARRAY,dest);
				int len=task[i]-task[i-1];
				mpi->send(len,dest);//len
				//cout<<"client("<<dest<<") "<<"end send len:"<<len<<endl;
				for(int j=task[i-1];j<task[i];j++){
					mpi->send(&xs[j][0],xs[j].size(),dest);
				}
				//cout<<"client("<<dest<<") "<<"end send ."<<endl;
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
				//population initializing....
				bestI=0;
				for(int i=0;i<PopSize;i++){
					for(int d=0;d<numDim;d++){
						x[i][d]=drand(range[d][0],range[d][1]);
					} 
				}

				evaluatePopulation(x,fx);
				calulateBestI();
				//update, main process.
				update(MaxFEs/PopSize-1);
				endEvaluate();//stop evaluating....
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
							//cout<<mpi->getName()<<" recv len:"<<len<<endl;
							for(int i=0;i<len;i++){
								mpi->recv(&x[i][0],numDim,0);
							}
							//	   cout<<mpi->getName()<<" end recv"<<len<<endl;
							for(int i=0;i<len;i++){
								fx[i]=f->evaluate(&x[i][0]);
							}
							//cout<<mpi->getName()<<" end evaluation"<<len<<endl;
							for(int i=0;i<len;i++){
								mpi->send(&fx[i],1,0);
							}
							//cout<<mpi->getName()<<" end send"<<len<<endl;
							break;
						default:
							cerr<<"Error:Unknown mess_type"<<endl;
							assert(false);
							break;
					}
				}
				return -1;
			}
		}
};
//int main(int argc,char *argv[]){
int old_main(int argc,char *argv[]){
	//srand(time(NULL));
	MPIHelper mpi(argc,argv);
	Save s(mpi.getName(),"x","y");
	//Save s("title","x","y");
	s.add(mpi.getID());
	//	s.add(100);
	//	ofstream out;
	//	out.open("abc.txt");
	//	out<<"hello from id:"<<mpi.getID()<<endl;
	//	out.close();
	Trace(mpi.getName());
	Trace(mpi.getID());
	if(mpi.isMaster()){
		for(int i=1;i<mpi.getNumProcesses();i++){
			mpi.send(i,i);
		}
	}else{
		int t;
		mpi.recv(t,0);
		Trace(t);
	}
}
///////////////////////
ParallelDE *de=0;
void IntHandler(int s){
	if(de!=NULL){
		de->endAEvaluate();
		delete de;
	}
	cout<<"My:Caught signal SIGINT"<<endl;
	exit(1);
}
void saveConfigData(int id,const char *f,const char *algorithm,int run,int MaxRun,int numOfProcesses,int MaxFEs,int PopSize,int numDim,double F,double CR,const char *state,double usedTime,double absError){
		ofstream runConfig;
		char buff[5000];
		sprintf(buff,"Run-configuration-%d.txt",id);
		runConfig<<"ID\tFunction\tAlgorithm\tRun\tMaxRun\tNumOfProcesses\tMaxFEs\tPopSize\tNumDim\tF-parameter\tCR-paramter\tState\tUsedTime\tAbsError"<<endl;
		sprintf(buff,"%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%g\t%s\t%g\t%g\n",id,f,algorithm,run,MaxRun,numOfProcesses,MaxFEs,
				PopSize,NumDim,F,CR,state,usedTime,absError);
		runConfig.flush();
		runConfig.close();
}

//int old_main(int argc,char *argv[]){
int main(int argc,char *argv[]){
	SignalHandleHelper::registerSignalHandler(IntHandler);
	//
	vector<double>x;
	double fx;
	//srand(time(NULL));
	MPIHelper mpi(argc,argv);
	//	SearchParam param("DE.json");
	SearchParam param("MyF1.json");
	de=new ParallelDE(&mpi);
	de->initParam(&param);
	//Function*f=new PECFunction(param.getInt("NumDim"));
	Function*f=new MyF1(param.getInt("NumDim"));
	/********Master data*************/
	if(mpi.isMaster()){
		cout<<"Runing "<<de->getName()<<" "<<endl;
		cout<<"FunName(MyBestF,Optima)"<<endl;
		//
		int configID=IDHelper::newID();
		int run=-1;//unused
		int MaxRun=-1;//unused
		/*
		runConfig<<"ID\tFunction\tAlgorithm\tRun\tMaxRun\tNumOfProcesses\tMaxFEs\tPopSize\tNumDim\tF-parameter\tCR-paramter\tState\tUsedTime\tAbsError"<<endl;
		sprintf(buff,"%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%g\t%g\t%s\t%g\t%g\n",configID,f->getName(),de->getName(),run,MaxRun,mpi->getNumProcesses(),param.getInt("MaxFEs"),
				param.getInt("PopSize"),param.getInt("NumDim"),F=param->getDouble("F"),CR=param->getDouble("CR"),state,usedTime,absError);


		   */

		runConfig.open(buff);
		char state[50]="start";
		double usedTime=-1;
		double absError=-1;
		Tic::tic("begin");
	}
	double res=de->getMax(f,param.getInt("MaxFEs"),x,fx);
	if(mpi.isMaster()){
		usedTime=Tic::tic("end");
		absError=fabs(fx-f->getFBest());
		strcpy(state,"end");

		printf("%s(%g,%g)",f->getName(),fx,f->getFBest());
	}
	delete de;
	return 0;
}

