#ifdef ID_HELPER_H
#else
#define ID_HELPER_H
#include<fstream>
using namespace std;
class IDHelper{
	public:
		static int newID(const char *file=0){
			if(file==0){
				file="ID.txt";
			}
			ifstream readid;
			readid.open(file);
			int id=-1;
			readid>>id;
			readid.close();
			id++;
			ofstream writeid;
			writeid.open(file);
			writeid<<id;
			writeid.close();
			return id;
		}
};


/*
   Example usage:
#include<iostream>
using namespace std;
int main(){
	cout<<"my id:"<<IDHelper::newID()<<endl;
	return 0;
}
*/
#endif
