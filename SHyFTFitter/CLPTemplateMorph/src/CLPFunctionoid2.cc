// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>
#include "CLPTemplateMorph/interface/CLPFunctionoid2.h"

using namespace std;


CLPFunctionoid2::CLPFunctionoid2()
{
}

// friends
ostream& operator<< (ostream& o_stream, const CLPFunctionoid2 &rhs)
{
   return o_stream;
} 
