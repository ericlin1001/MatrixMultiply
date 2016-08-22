#ifdef PEC_FUNCTION_H
#define PEC_FUNCTION_H
#else
#include "include/functions.h"
#include "cir.h"
#define UNUSED 1
DefFunction(PECFunction,UNUSED,UNUSED,0)
	return _objective(xs);
EndDef
#endif

