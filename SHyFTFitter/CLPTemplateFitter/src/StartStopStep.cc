// -*- C++ -*-
#include <iostream>
#include <fstream>
#include <iomanip>

#include "TSystem.h"

#include "CLPTemplateFitter/interface/StartStopStep.h"

using namespace std;

StartStopStep::StartStopStep () :
   m_start (0.), m_stop (0.), m_step (0.)
{
}

StartStopStep::StartStopStep (double start, double stop, double step) :
   m_start (start), m_stop (stop), m_step (step)
{
}

// friends
ostream& operator<< (ostream& o_stream, const StartStopStep &rhs)
{
   o_stream << Form ("[%6.2f, %6.2f] - %5.3f", 
                     rhs.start(), rhs.stop(), rhs.step());
                     
   return o_stream;
} 
