#pragma once
#include "simData.h"
class InitialFunc { public: Data * data; InitialFunc(Data * d) : data(d) {} virtual ~InitialFunc() {} };
