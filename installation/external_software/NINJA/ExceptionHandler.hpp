#ifndef ExceptionHandler_HPP
#define ExceptionHandler_HPP

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <time.h>
#include <errno.h>


namespace Exception{
	void critical();
	void criticalErrno(const char* arg);
}

#endif
