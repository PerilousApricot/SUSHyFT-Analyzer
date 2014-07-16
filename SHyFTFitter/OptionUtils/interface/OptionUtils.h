// -*- C++ -*-

#if !defined(OptionUtils_H)
#define OptionUtils_H

#include <map>
#include <vector>
#include <string>

namespace OptionUtils
{
   /////////////////////
   // /////////////// //
   // // Constants // //
   // /////////////// //
   /////////////////////

   // typedefs
   typedef std::vector< std::string >           SVec;
   typedef std::vector< int >                   IVec;
   typedef std::vector< double >                DVec;
   typedef std::map< std::string, void* >       OptionMap;
   typedef std::map< std::string, bool  >       StrBoolMap;
   typedef std::map< std::string, std::string > StrStrMap;
   typedef SVec::const_iterator                 SVecConstIter;
   typedef IVec::const_iterator                 IVecConstIter;
   typedef DVec::const_iterator                 DVecConstIter;
   typedef OptionMap::iterator                  OptionMapIter;
   typedef StrBoolMap::const_iterator           StrBoolMapConstIter;

   // constants
   const std::string kSpaces = " \t";

   ///////////////
   // Functions //
   ///////////////

   // parse the command line arguments.  If 'returnArgs' is true, then
   // any non-assignments and non-options will be returned.
   SVec parseArguments (int argc, char** argv, bool returnArgs = false);

   // set a usage string for '--help' option
   void setUsageString (const std::string &usage);

   // prints out '--help' screen, then exits.
   void help();

   // remove an ending (e.g., '.root') from a string
   std::string removeEnding (const std::string &input, 
                                    const std::string &ending);   

   // splits a line into words
   void split (SVec &retval, std::string line, std::string match = " \t",
               bool ignoreComments = true);

   // removes '# ....' comment
   void removeComment (std::string &line);

   // removes leading and trailing spaces
   void removeLeadingAndTrailingSpaces (std::string &line);

   // given a line, finds first non-space word and rest of line
   void findCommand (const std::string &line,
                     std::string &command,
                     std::string &rest);

   // print out all of the variables hooked up
   void printOptionValues();

   // converts a string to lower case characters
   void lowercaseString(std::string &arg); 

   // converts a single character to lower case
   char toLower (char &ch);

   // add variable to option maps
   void addOptionKey (std::string key, int         &variable, 
                      const std::string &description = "");
   void addOptionKey (std::string key, double      &variable, 
                      const std::string &description = "");
   void addOptionKey (std::string key, bool        &variable, 
                      const std::string &description = "");
   void addOptionKey (std::string key, std::string &variable, 
                      const std::string &description = "");
   void addOptionKey (std::string key, SVec        &variable, 
                      const std::string &description = "");
   void addOptionKey (std::string key, IVec        &variable, 
                      const std::string &description = "");
   void addOptionKey (std::string key, DVec        &variable, 
                      const std::string &description = "");

   // returns true if a variable has been modified from the command
   // line.
   bool valueHasBeenModified (const std::string &key);

   // Sets a variable 'varname' to a 'value' from a string
   // 'varname=value'.  If 'dontOverrideChange' is set true, then the
   // function will NOT set a variable that has been already set.
   // This will allow you to read in a file to set most variables and
   // still be allowed to make changes from the command line.  If you
   // want to ignore the first 'n' characters, simply set 'offset=n'.
   bool setVariableFromString (const std::string &arg,
                                      bool dontOverrideChange = false,
                                      int offset = 0);

   // sets variable options from file where lines are formatted as
   // - var=value
   bool setVariablesFromFile (const std::string &filename);

   /////////////////////
   // /////////////// //
   // // Variables // //
   // /////////////// //
   /////////////////////

   extern OptionMap    ns_intMap;
   extern OptionMap    ns_doubleMap;
   extern OptionMap    ns_stringMap;
   extern OptionMap    ns_svecMap;
   extern OptionMap    ns_ivecMap;
   extern OptionMap    ns_dvecMap;
   extern OptionMap    ns_boolMap;
   extern StrBoolMap   ns_variableModifiedMap;
   extern StrStrMap    ns_variableDescriptionMap;
   extern std::string  ns_usageString;  
   extern std::string  ns_argv0;  

};

#endif // OptionUtils_H
