//Define USE_MPI to use mpi.
//without USE_MPI, mpi is not used.
#ifdef MPI_HELPER_H
#else
#define MPI_HELPER_H

#ifdef USE_MPI
//be careful of, this block should be the first code in main file.
#define OMPI_IMPORTS
#include "mpi.h"
#endif
#include "stdio.h"

class MPIHelper{

#define TAG   100
#define TAG1  101
#define TAG2  102
#define TAG3  103
#define TAG4  104

	int id;
	int numProcess;
#ifdef USE_MPI
	char name[MPI_MAX_PROCESSOR_NAME];
#else
	char name[100];
#endif

#ifdef USE_MPI
	MPI_Status status;
	MPI_Request request;
#endif
	public:
		MPIHelper(int &argc,char ** &argv){
			char buff[100]="localhost";
#ifdef USE_MPI
			int len;
			MPI_Init(&argc,&argv);
			MPI_Comm_rank(MPI_COMM_WORLD,&id);
			MPI_Comm_size(MPI_COMM_WORLD,&numProcess);
			MPI_Get_processor_name(buff,&len);
#else
		id=0;
		numProcess=1;
#endif
		sprintf(name,"%s(%d)",buff,id);
		}
		const char *getName()const{
			return name;
		}
		bool isMaster()const{return id==0;}
		int getID()const{return id;}
		int getNumProcesses()const{return numProcess;}


		void send(double *buf,int count,int dest,int tag=TAG){
#ifdef USE_MPI
			MPI_Send(buf,count,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
#endif
		}
		void send(int *buf,int count,int dest,int tag=TAG1){
#ifdef USE_MPI
			MPI_Send(buf,count,MPI_INT,dest,tag,MPI_COMM_WORLD);
#endif
		}

		void send(int buf,int dest,int tag=TAG2){
#ifdef USE_MPI
			MPI_Send(&buf,1,MPI_INT,dest,tag,MPI_COMM_WORLD);
#endif
		}

		void send(double buf,int dest,int tag=TAG2){
#ifdef USE_MPI
			MPI_Send(&buf,1,MPI_DOUBLE,dest,tag,MPI_COMM_WORLD);
#endif
		}

		void asend(int buf,int dest,int tag=TAG3){
#ifdef USE_MPI
			MPI_Isend(&buf,1,MPI_INT,dest,tag,MPI_COMM_WORLD,&request);
#endif
		}
		/*
		//asychronized.
		void asend(int buf,int dest,int tag=TAG3){
#ifdef USE_MPI
			MPI_Send(&buf,1,MPI_INT,dest,tag,MPI_COMM_WORLD,&status);
#endif
		}
		*/

		void recv(double *buf,int count,int src,int tag=TAG){
#ifdef USE_MPI
			MPI_Recv(buf,count,MPI_DOUBLE,src,tag,MPI_COMM_WORLD,&status);
#endif
		}
		void recv(int *buf,int count,int src,int tag=TAG1){
#ifdef USE_MPI
			MPI_Recv(buf,count,MPI_INT,src,tag,MPI_COMM_WORLD,&status);
#endif
		}
		void recv(int &buf,int src,int tag=TAG2){
#ifdef USE_MPI
			MPI_Recv(&buf,1,MPI_INT,src,tag,MPI_COMM_WORLD,&status);
#endif
		}
		void recv(double &buf,int src,int tag=TAG2){
#ifdef USE_MPI
			MPI_Recv(&buf,1,MPI_DOUBLE,src,tag,MPI_COMM_WORLD,&status);
#endif
		}
/*
		const MPI_Status *getMPIStatus()const{
			return &status;
		}
*/
		~MPIHelper(){
#ifdef USE_MPI
			MPI_Finalize();
#endif
		}
};
#endif

