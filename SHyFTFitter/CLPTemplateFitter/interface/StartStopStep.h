// -*- C++ -*-

#if !defined(StartStopStep_HH)
#define StartStopStep_HH

#include <vector>

class StartStopStep
{
   public:
      //////////////////////
      // Public Constants //
      //////////////////////

      /////////////
      // friends //
      /////////////
      // tells particle data how to print itself out
      friend std::ostream& operator<< (std::ostream& o_stream, 
                                       const StartStopStep &rhs);

      typedef std::vector< StartStopStep > Vec;

      //////////////////////////
      //            _         //
      // |\/|      |_         //
      // |  |EMBER | UNCTIONS //
      //                      //
      //////////////////////////

      /////////////////////////////////
      // Constructors and Destructor //
      /////////////////////////////////
      StartStopStep ();
      StartStopStep (double start, double stop, double step);

       ////////////////
      // One Liners //
      ////////////////
      //////////////////////
      // Access Functions //
      //////////////////////
      double start() const { return m_start; }
      double stop()  const { return m_stop; }
      double step()  const { return m_step; }

      ///////////////////////
      // Setting Functions //
      ///////////////////////
      void setStart (double start) { m_start = start; }
      void setStop  (double stop)  { m_stop  = stop; }
      void setStep  (double step)  { m_step  = step; }
      void setStartStopStep (double start, double stop, double step)
      { m_start = start; m_stop = stop; m_step = step;}

      //////////////////////////////
      // Regular Member Functions //
      //////////////////////////////

      /////////////////////////////
      // Static Member Functions //
      /////////////////////////////

  protected:
      /////////////////////////
      // Protected Constants //
      /////////////////////////

      ///////////////////////////
      // Protected Member Data //
      ///////////////////////////

  private:
      ///////////////////////
      // Private Constants //
      ///////////////////////

      /////////////////////////
      // Private Member Data //
      /////////////////////////
      double m_start;
      double m_stop;
      double m_step;
};

#endif // StartStopStep_HH
