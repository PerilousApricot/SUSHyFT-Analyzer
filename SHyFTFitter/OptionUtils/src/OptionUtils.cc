// -*- C++ -*-
#include <iostream>
#include <string>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <iomanip>

#include "OptionUtils/interface/OptionUtils.h"
#include "OptionUtils/interface/dout.h"
#include "OptionUtils/interface/dumpSTL.icc"

using namespace std;

// Variable definitions
OptionUtils::OptionMap  OptionUtils::ns_intMap;
OptionUtils::OptionMap  OptionUtils::ns_doubleMap;
OptionUtils::OptionMap  OptionUtils::ns_stringMap;
OptionUtils::OptionMap  OptionUtils::ns_svecMap;
OptionUtils::OptionMap  OptionUtils::ns_ivecMap;
OptionUtils::OptionMap  OptionUtils::ns_dvecMap;
OptionUtils::OptionMap  OptionUtils::ns_boolMap;
OptionUtils::StrBoolMap OptionUtils::ns_variableModifiedMap;
OptionUtils::StrStrMap  OptionUtils::ns_variableDescriptionMap;
string                  OptionUtils::ns_usageString;
string                  OptionUtils::ns_argv0;

OptionUtils::SVec
OptionUtils::parseArguments (int argc, char** argv, bool returnArgs)
{   
   bool printOptions = false;
   SVec argsVec;
   ns_argv0 = argv[0];
   for (int loop = 1; loop < argc; ++loop)
   {
      string arg = argv[loop];
      string::size_type where = arg.find_first_of("=");
      if (string::npos != where)
      {
         if (! setVariableFromString (arg) )
         {
            cerr << "Don't understand: " << arg << endl;
            exit(0);
         }
         continue;
      } // tag=value strings
      else if (arg.at(0) == '-')
      {
         string::size_type where = arg.find_first_not_of("-");
         if (string::npos == where)
         {
            // a poorly formed option
            cerr << "Don't understand: " << arg << endl;
            exit(0);
            continue;
         }
         lowercaseString (arg);
         char first = arg.at (where);
         // Print the values
         if ('p' == first)
         {
            printOptions = true;
            continue;
         }
         // Exit after printing values
         if ('e' == first || 'h' == first)
         {
            help();
         }
         // if we're still here, then we've got a problem.
         cerr << "Don't understand: " << arg << endl;
         exit(0);
      } // -arg strings
      if (returnArgs)
      {
         argsVec.push_back (arg);
      } else {
         cerr << "Don't understand: " << arg << endl;
         exit(0);
      }
   } // for loop
   if (printOptions)
   {
      printOptionValues();
   } // if printOptions
   return argsVec;
}

void 
OptionUtils::setUsageString (const std::string &usage) 
{ 
   ns_usageString = usage; 
}

void
OptionUtils::help()
{
   if (ns_usageString.length())
   {
      cout << "Usage: " << ns_argv0 << " " << ns_usageString << endl;
   }
   printOptionValues();
   exit (0);
}

void 
OptionUtils::split (SVec &retval, string line, string match, 
                    bool ignoreComments)
{
   if (ignoreComments)
   {
      removeComment (line);
   } // if ignoreComments
   retval.clear();
   // find the first non-space
   string::size_type start1 = line.find_first_not_of(kSpaces);
   // Is the first non-space character a '#'
   char firstCh = line[start1];
   if ('#' == firstCh)
   {
      // this line is a comment
      return;
   }

   line += match; // get last word of line
   string::size_type last = string::npos;
   string::size_type current = line.find_first_of(match);
   while (string::npos != current)
   {
      string::size_type pos;
      if (string::npos != last)
      {
         pos = last + 1;
      } else {
         pos = 0;
      }
      string part = line.substr( pos, current - last - 1);
      // don't bother adding 0 length strings
      if (part.length()) 
      {
         retval.push_back(part);
      }
      last = current;
      current = line.find_first_of(match, current + 1);
   } // while we're finding spaces
}

void
OptionUtils::removeComment (string &line)
{
   string::size_type location = line.find ("#");
   if (string::npos != location)
   {
      // we've got a comment.  Strip it out
      line = line.substr (0, location - 1);
   } // if found
}

void
OptionUtils::removeLeadingAndTrailingSpaces (std::string &line)
{
   string::size_type pos = line.find_first_not_of (kSpaces);
   if (string::npos == pos)
   {
      // we don't have anything here at all.  Just quit now
      return;
   }
   if (pos)
   {
      // We have spaces at the beginning.
      line = line.substr (pos);
   }
   pos = line.find_last_not_of (kSpaces);
   if (pos + 1 != line.length())
   {
      // we've got spaces at the end
      line = line.substr (0, pos + 1);
   }
}
 
string
OptionUtils::removeEnding (const string &input, const string &ending)
{
   string::size_type position = input.rfind(ending);
   if (input.length() - ending.length() == position)
   {
      // we've got it
      return input.substr(0, position);
   }
   // If we're still here, it wasn't there
   return input;
}

void 
OptionUtils::findCommand (const string &line,
                          string &command,
                          string &rest)
{
   command = rest = "";
   string::size_type nonspace = line.find_first_not_of (kSpaces);
   if (string::npos == nonspace)
   {
      // we don't have anything here at all.  Just quit now
      return;
   }
   string::size_type space = line.find_first_of (kSpaces, nonspace);
   if (string::npos == space)
   {
      // we only have a command and nothing else
      command = line.substr (nonspace);
      return;
   }
   command = line.substr (nonspace, space - nonspace);
   rest    = line.substr (space + 1);
   removeLeadingAndTrailingSpaces (rest);
}


void
OptionUtils::printOptionValues()
{
   cout << "------------------------------------------------------------------" 
        << endl << "Option Values:" << endl;
   // Print the bools first
   if (ns_boolMap.size())
   {
      cout << "  Bool options:" << endl;
   }
   for (OptionMapIter iter = ns_boolMap.begin(); 
       ns_boolMap.end() != iter; 
       ++iter) 
   {
      string &description = ns_variableDescriptionMap[ iter->first ];
      bool *valuePtr = (bool *) iter->second;
      cout << "    " << iter->first << " = ";
      if (*valuePtr)
      {
         cout << "true";
      } else {
         cout << "false";
      }
      if (description.length())
      {
         cout << " - " << description;
      }
      cout << endl;
   } // for iter

   // Print the ints next
   if (ns_intMap.size())
   {
      cout << "  Int options:" << endl;
   }
   for (OptionMapIter iter = ns_intMap.begin(); 
       ns_intMap.end() != iter; 
       ++iter) 
   {
      string &description = ns_variableDescriptionMap[ iter->first ];
      int *valuePtr = (int *)iter->second;
      cout << "    " << iter->first << " = " << *valuePtr;
      if (description.length())
      {
         cout << " - " << description;
      }
      cout << endl;
   } // for iter

   // Print the doubles next
   if (ns_doubleMap.size())
   {
      cout << "  Double options:" << endl;
   }
   for (OptionMapIter iter = ns_doubleMap.begin(); 
       ns_doubleMap.end() != iter; 
       ++iter) 
   {
      string &description = ns_variableDescriptionMap[ iter->first ];
      double *valuePtr = (double *)iter->second;
      cout << "    " << iter->first << " = " << *valuePtr;
      if (description.length())
      {
         cout << " - " << description;
      }
      cout << endl;
   } // for iter

   // Print the strings next
   if (ns_stringMap.size())
   {
      cout << "  String options:" << endl;
   }
   for (OptionMapIter iter = ns_stringMap.begin(); 
       ns_stringMap.end() != iter; 
       ++iter) 
   {
      string &description = ns_variableDescriptionMap[ iter->first ];
      string *valuePtr = (string *) iter->second;
      cout << "    " << iter->first << " = '" << *valuePtr << "'";
      if (description.length())
      {
         cout << " - " << description;
      }
      cout << endl;
   } // for iter
   // SVec
   if (ns_svecMap.size())
   {
      cout << "  SVec options:" << endl;
   }
   for (OptionMapIter iter = ns_svecMap.begin(); 
       ns_svecMap.end() != iter; 
       ++iter) 
   {
      string &description = ns_variableDescriptionMap[ iter->first ];
      SVec *valuePtr = (SVec *) iter->second;
      cout << "    " << iter->first << " = ";
      dumpSTL (*valuePtr); 
      cout << endl;
      if (description.length())
      {
         cout << "      - " << description;
      }
      cout << endl;
   } // for iter
   // IVec
   if (ns_ivecMap.size())
   {
      cout << "  IVec options:" << endl;
   }
   for (OptionMapIter iter = ns_ivecMap.begin(); 
       ns_ivecMap.end() != iter; 
       ++iter) 
   {
      string &description = ns_variableDescriptionMap[ iter->first ];
      IVec *valuePtr = (IVec *) iter->second;
      cout << "    " << iter->first << " = ";
      dumpSTL (*valuePtr);
      cout << endl;
      if (description.length())
      {
         cout << "      - " << description;
      }
      cout << endl;
   } // for iter
   // DVec
   if (ns_dvecMap.size())
   {
      cout << "  DVec options:" << endl;
   }
   for (OptionMapIter iter = ns_dvecMap.begin(); 
       ns_dvecMap.end() != iter; 
       ++iter) 
   {
      string &description = ns_variableDescriptionMap[ iter->first ];
      DVec *valuePtr = (DVec *) iter->second;
      cout << "    " << iter->first << " = ";
      dumpSTL (*valuePtr);
      cout << endl;
      if (description.length())
      {
         cout << "      - " << description;
      }
      cout << endl;
   } // for iter
   cout << "------------------------------------------------------------------" 
        << endl;
}

void 
OptionUtils::lowercaseString(string& arg)
{
   // assumes 'toLower(ch)' modifies ch
   std::for_each (arg.begin(), arg.end(), OptionUtils::toLower);
   // // assumes 'toLower(ch)' returns the lower case char
   // std::transform (arg.begin(), arg.end(), arg.begin(), 
   //                 OptionUtils::toLower);
}

char 
OptionUtils::toLower (char& ch)
{
   ch = tolower (ch);
   return ch;
}

void
OptionUtils::addOptionKey (string key, int &variable, 
                           const string &description)
{
   lowercaseString (key);
   ns_intMap[key] = (void*) &variable;
   ns_variableModifiedMap[key] = false;
   ns_variableDescriptionMap[key] = description;
}

void
OptionUtils::addOptionKey (string key, double &variable, 
                           const string &description)
{
   lowercaseString (key);
   ns_doubleMap[key] = (void*) &variable;
   ns_variableModifiedMap[key] = false;
   ns_variableDescriptionMap[key] = description;
}

void
OptionUtils::addOptionKey (string key, bool &variable, 
                           const string &description)
{
   lowercaseString (key);
   ns_boolMap[key] = (void*) &variable;
   ns_variableModifiedMap[key] = false;
   ns_variableDescriptionMap[key] = description;
}

void
OptionUtils::addOptionKey (string key, string &variable, 
                           const string &description)
{
   lowercaseString (key);
   ns_stringMap[key] = (void*) &variable;
   ns_variableModifiedMap[key] = false;
   ns_variableDescriptionMap[key] = description;
}

void
OptionUtils::addOptionKey (string key, SVec &variable, 
                           const string &description)
{
   lowercaseString (key);
   ns_svecMap[key] = (void*) &variable;
   // Although we will never mark this variable as modified, we still
   // need to put an entry into the map.
   ns_variableModifiedMap[key] = false;
   ns_variableDescriptionMap[key] = description;
}

void
OptionUtils::addOptionKey (string key, IVec &variable, 
                           const string &description)
{
   lowercaseString (key);
   ns_ivecMap[key] = (void*) &variable;
   // Although we will never mark this variable as modified, we still
   // need to put an entry into the map.
   ns_variableModifiedMap[key] = false;
   ns_variableDescriptionMap[key] = description;
}

void
OptionUtils::addOptionKey (string key, DVec &variable, 
                           const string &description)
{
   lowercaseString (key);
   ns_dvecMap[key] = (void*) &variable;
   // Although we will never mark this variable as modified, we still
   // need to put an entry into the map.
   ns_variableModifiedMap[key] = false;
   ns_variableDescriptionMap[key] = description;
}

bool
OptionUtils::valueHasBeenModified (const string &key)
{
   StrBoolMapConstIter iter = ns_variableModifiedMap.find (key);
   if (ns_variableModifiedMap.end() == iter)
   {
      // Not found.  Not a valid option
      cerr << "OptionUtils::valueHasBeenModfied () Error: '" 
           << key << "' is not a valid key." << endl;
      return false;
   }
   return iter->second;
}

bool 
OptionUtils::setVariableFromString (const string &arg,
                                    bool dontOverrideChange,
                                    int offset)
{
   string::size_type where = arg.find_first_of("=", offset + 1);
   string varname = arg.substr (offset, where - offset);
   string value   = arg.substr (where + 1);
   lowercaseString (varname);
   // check to make sure this is a valid option
   StrBoolMapConstIter sbiter = ns_variableModifiedMap.find (varname);
   if (ns_variableModifiedMap.end() == sbiter)
   {
      // Not found.  Not a valid option
      return false;
   }
  
   // if 'dontOverrideChange' is set, then we are being asked to NOT
   // change any variables that have already been changed.
   if (dontOverrideChange && valueHasBeenModified (varname) )
   {
      // don't go any further
      return true;
   }
   // doubles
   OptionMapIter iter = ns_doubleMap.find(varname);
   if (ns_doubleMap.end() != iter)
   {
      // we found it
      double *doublePtr = (double*) iter->second;
      *doublePtr = atof(value.c_str());
      ns_variableModifiedMap[varname] = true;
      return true;
   }
   // integers
   iter = ns_intMap.find(varname);
   if (ns_intMap.end() != iter)
   {
      // we found it
      int *intPtr = (int*) iter->second;
      // use 'atof' instead of 'atoi' to get scientific notation
      *intPtr = (int) atof( value.c_str() );
      ns_variableModifiedMap[varname] = true;
      return true;
   }
   // strings
   iter = ns_stringMap.find(varname);
   if (ns_stringMap.end() != iter)
   {
      // we found it
      string *stringPtr = (string*) iter->second;
      *stringPtr = value;
      ns_variableModifiedMap[varname] = true;
      return true;
   }
   // SVec
   iter = ns_svecMap.find(varname);
   if (ns_svecMap.end() != iter)
   {
      // we found it
      SVec *svecPtr = (SVec*) iter->second;
      SVec words;
      split (words, value, ",");
      for (SVecConstIter svecIter = words.begin();
           words.end() != svecIter;
           ++svecIter)
      {
         svecPtr->push_back( *svecIter );
      }
      // we don't want to mark this as modified because we can add
      // many values to this
      // ns_variableModifiedMap[varname] = true;
      return true;
   }
   // IVec
   iter = ns_ivecMap.find(varname);
   if (ns_ivecMap.end() != iter)
   {
      // we found it
      IVec *ivecPtr = (IVec*) iter->second;
      SVec words;
      split (words, value, ",");
      for (SVecConstIter svecIter = words.begin();
           words.end() != svecIter;
           ++svecIter)
      {
         ivecPtr->push_back( atoi( svecIter->c_str() ) );
      }
      // we don't want to mark this as modified because we can add
      // many values to this
      // ns_variableModifiedMap[varname] = true;
      return true;
   }
   // DVec
   iter = ns_dvecMap.find(varname);
   if (ns_dvecMap.end() != iter)
   {
      // we found it
      DVec *dvecPtr = (DVec*) iter->second;
      SVec words;
      split (words, value, ",");
      for (SVecConstIter svecIter = words.begin();
           words.end() != svecIter;
           ++svecIter)
      {
         dvecPtr->push_back( atof( svecIter->c_str() ) );
      }
      // we don't want to mark this as modified because we can add
      // many values to this
      // ns_variableModifiedMap[varname] = true;
      return true;
   }
   // bools
   iter = ns_boolMap.find(varname);
   if (ns_boolMap.end() != iter)
   {
      // we found it
      bool *boolPtr = (bool*) iter->second;
      int val = atoi(value.c_str());
      if (val)
      {
         *boolPtr = true;
      } else {
         *boolPtr = false;
      }               
      ns_variableModifiedMap[varname] = true;
      return true;
   }
   // We didn't find your variable.
   return false;
}

bool
OptionUtils::setVariablesFromFile (const string &filename)
{
   ifstream source(filename.c_str(), ios::in);
   if (! source)
   {
      cerr << "file " << filename << "could not be opened" << endl;
      return false;
   }
   string line;
   while (getline (source, line))
   {
      // find the first nonspace
      string::size_type where = line.find_first_not_of(kSpaces);
      if (string::npos == where)
      {
         // no non-spaces
         continue;
      } 
      char first = line.at (where);
      if ('-' != first)
      {
         continue;
      }
      where = line.find_first_not_of(kSpaces, where + 1);
      if (string::npos == where)
      {
         // no non-spaces
         continue;
      }
      // Get string starting at first nonspace after '-'.  Copy it to
      // another string without copying any spaces and stopping at the
      // first '#'.
      string withspaces = line.substr (where);
      string nospaces;
      for (int position = 0; 
           position < (int) withspaces.length(); 
           ++position)
      {
         char ch = withspaces[position];
         if ('#' == ch)
         {
            // start of a comment
            break;
         } else if (' ' == ch || '\t' == ch)
         {
            continue;
         }
         nospaces += ch;
      } // for position
      if (! setVariableFromString (nospaces, true) )
      {
         cerr << "Don't understand line" << endl << line << endl
              << "in options file '" << filename << "'.  Aborting."
              << endl;
         exit(0);
      } // if setting variable failed
   } // while getline
   return true;
}
