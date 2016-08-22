#ifdef SIGNAL_HANDLE_HELPER
#else
#define SIGNAL_HANDLE_HELPER
#include<signal.h>
class SignalHandleHelper{
	private:
		SignalHandleHelper(){
		}
	public:
typedef void SignalHandler(int);
static void registerSignalHandler(SignalHandler h,int type=SIGINT){
			struct sigaction sigIntHandler;
			sigIntHandler.sa_handler=h;
			sigemptyset(&sigIntHandler.sa_mask);
			sigIntHandler.sa_flags=0;
			sigaction(type,&sigIntHandler,NULL);
		}
};
/****example handler:
void myHandler(int s){
	printf("Caught signal %d\n",s);
	exit(0);
}
****/
#endif

