#ifdef BASIC_DE
#else
#define BASIC_DE
#include<iostream>
#include<fstream>
#include "template.h"

using namespace std;


typedef SettingParser SearchParam ;

class BasicDE:public EA{
	private:
		//about function:f
		Function *f;
		double xmin;
		double xmax;
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
		//
	private:
		void updateX(){
			//main process
			RandomPermutation perm(PopSize);
			for(int i=0;i<PopSize;i++){
				perm.generate();
				int a=bestI; int b=perm.next(); int c=perm.next();
				int randDim=rand()%numDim;
				for(int j=0;j<numDim;j++){
					if(j==randDim||drand()<CR){
						tx[j]=x[a][j]+F*(x[b][j]-x[c][j]);
						if(tx[j]<xmin || tx[j]>xmax){
							tx[j]=drand(xmin,xmax);
						}
					}else{
						tx[j]=x[i][j];
					}
				}
				double ftx=f->evaluate(&tx[0]);
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
#undef SaveData
#define SaveData s.add(getBestFx());
			Save s(f->getName(),"Generation","F");
			SaveData;
			for(int g=1;g<=maxGeneration;g++){
				updateX();
				SaveData;
			}
#undef SaveData
		}
	public:
		virtual void initParam(SearchParam *param){
			PopSize=param->getInt("PopSize");
			F=param->getDouble("F");
			CR=param->getDouble("CR");
			setName(param->getString("Name"));
		}
		virtual double getMin(Function *f,int MaxFEs,vector<double>&out_x,double &out_fx){
			this->f=f;
			numDim=f->getNumDim();
			xmin=f->getRange(0);
			xmax=f->getRange(1);
			//population initializing....
			//allocate space.
			tx.resize(numDim);
			x.resize(PopSize);
			fx.resize(PopSize);
			for(int i=0;i<PopSize;i++){
				x[i].resize(numDim);
			}
			bestI=0;
			for(int i=0;i<PopSize;i++){
				for(int d=0;d<numDim;d++){
					x[i][d]=drand(xmin,xmax);
				} 
				fx[i]=f->evaluate(&x[i][0]);
				if(fx[i]<fx[bestI]){
					bestI=i;
				}
			}
			//update, main process.
			update(MaxFEs/PopSize-1);
			out_x=x[bestI];
			out_fx=fx[bestI];
			return out_fx;
		}
};
//read Range from param.
class BasicDE1:public EA{
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
		//
	private:
		void updateX(){
			//main process
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
					double ftx=f->evaluate(&tx[0]);
					if(ftx<fx[i]){
						x[i]=tx;
						fx[i]=ftx;
						if(ftx<fx[bestI]){
							bestI=i;
						}
					}
				}
			}
		}
	private:
		double getBestFx()const{
			return fx[bestI];
		}
		void update(int maxGeneration){
#undef SaveData
#define SaveData s.add(getBestFx());
			Save s(f->getName(),"Generation","F");
			SaveData;
			for(int g=1;g<=maxGeneration;g++){
				updateX();
				SaveData;
				//Tic::tic("update one X");
			}
#undef SaveData
		}
	public:
		void initParam(SearchParam *param){
			PopSize=param->getInt("PopSize");
			F=param->getDouble("F");
			CR=param->getDouble("CR");
			param->getBiVector("Range",range);
			setName(param->getString("Name"));
		}
		virtual double getMin(Function *f,int MaxFEs,vector<double>&out_x,double &out_fx){
			this->f=f;
			numDim=f->getNumDim();
			//population initializing....
			//allocate space.
			tx.resize(numDim);
			x.resize(PopSize);
			fx.resize(PopSize);
			for(int i=0;i<PopSize;i++){
				x[i].resize(numDim);
			}
			bestI=0;
			for(int i=0;i<PopSize;i++){
				for(int d=0;d<numDim;d++){
					x[i][d]=drand(range[d][0],range[d][1]);
				} 
				fx[i]=f->evaluate(&x[i][0]);
				if(fx[i]<fx[bestI]){
					bestI=i;
				}
			}
			//update, main process.
			update(MaxFEs/PopSize-1);
			out_x=x[bestI];
			out_fx=fx[bestI];
			return out_fx;
		}
};
#endif

