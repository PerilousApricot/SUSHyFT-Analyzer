#! /usr/bin/env python

import optparse
import ROOT
import re
import sys
import os
import pprint


def shlibDir():
    return "/".join( os.path.abspath(sys.argv[0]).split('/')[:-2] ) + '/shlib'


def envString (shlibDir):
    cmssw_base = os.environ.get('CMSSW_BASE')
    if not cmssw_base:
        raise RuntimeError, "You must have CMSSW setup first"
    #ld_lib_path = os.environ.get('LD_LIBRARY_PATH')
    return '- env = cd %s; . /uscmst1/prod/sw/cms/bashrc prod; eval `scramv1 runtime -sh`; export LD_LIBRARY_PATH=%s:$LD_LIBRARY_PATH; cd -' % \
           (cmssw_base, shlibDir)


def executableNameAndDir (options):
    command, directory = '',''
    if options.exe:
        command = os.path.abspath (options.exe)
        if not os.path.exists (command):
            command += '.exe'
    else:
        # looks like we'd better guess
        directory = os.environ['PWD']
        # we'll look in the current directory and its two parents
        for tries in range (3):
            base = os.path.basename (directory)
            command = "%s/%s.exe" % (directory, base)
            if os.path.exists (command):
                break
            # if we're still here, then we need move up a directory
            # and try again
            directory = '/'.join( directory.split('/')[:-1] )
    if not directory:
        directory = os.path.dirname (command)
    if not os.path.exists (command):
        raise RuntimeError, "Command '%s' not found." % command
    return command, directory
      

if __name__ == "__main__":
    # Setup options parser
    parser = optparse.OptionParser \
             ("usage: %prog [options] ")
    parser.add_option ('--numPEs', dest='numPEs', type="string", default='5e2',
                       help="Number of PEs per section (default: %default)")
    parser.add_option ('--sections', dest='sections', type="int", default=20,
                       help="Number of sections (default: %default)")
    parser.add_option ('--seed', dest='seed', type="string",
                       default='randomSeed',
                       help="Name of random seed option (default: %default)")
    parser.add_option ('--lum', dest='lum', action="append",
                       type="string", default=[],
                       help="Comma separated luminosity numbers")
    parser.add_option ('--output', dest='output', type="string",
                       default='',
                       help="What should be prepended to output (default: %default)")
    parser.add_option ('--other', dest='other', action="append",
                       type="string", default=[],
                       help="Other options to be passed to fitter")
    parser.add_option ('--exe', dest='exe', type="string",
                       default='',
                       help="Executable name (if different than current or parent directory name)")
    parser.add_option ('--commandFile', dest='commandFile', type="string",
                       default='jobs.cmd',
                       help="Command file name (default '%default'")
    parser.add_option ('--newCommand', dest='newCommand', action='store_true',
                       help="Force new command file (i.e., do not append)")
    parser.add_option ('--shlibDir', dest='shlibDir', type="string",
                       default=shlibDir(),
                       help="Shlib directory (default '%default'")
    parser.add_option ('--doubleColon', dest='doubleColon', action='store_true',
                       help="Use '::' as delimeter")
    options, args = parser.parse_args()
    commaRE = re.compile (r'\s*,\s*')
    colonRE = re.compile (r'\s*:\s*')
    if options.doubleColon:
        colonRE = re.compile (r'\s*::\s*')
    lumis = []
    for chunk in options.lum:
        lumis.extend( commaRE.split (chunk) )
    if not lumis:
        raise RuntimeError, "You must provide at least one luminosity number"
    configs = []
    nameSet = set()
    tagSet  = set()
    for chunk in args:
        chunk = chunk.strip()
        pieces = colonRE.split (chunk)
        if len(pieces) != 2:
            raise RuntimeError, "Do not understand arg '%s'" % chunk
        if pieces[0] in nameSet:
            raise RuntimeError, "You can't use config '%s' twice" % pieces[0]
        if pieces[1] in tagSet:
            raise RuntimeError, "You can't use tag '%s' twice" % pieces[1]
        nameSet.add (pieces[0])
        tagSet.add  (pieces[1])
        configs.append( tuple( pieces ) )
    baseCommand, directory = executableNameAndDir (options)
    # prepare command file
    if not os.path.exists (options.commandFile) or options.newCommand:
        # create a file from scratch
        target = open (options.commandFile, 'w')
        target.write ('# -*- sh -*-\n%s\n' % envString (options.shlibDir))
        target.write ('\n\n# Commands:\n')
    else:
        target = open (options.commandFile, 'a')
    for config in configs:
        for lum in lumis:            
            partTag = "%s_lum%s" % (config[1], lum)
            noLumiTag = "%s" % (config[1])
            for section in range (1, options.sections + 1):
                log = "%s_%03d" % (partTag, section)
                command = "%-30s   %s " % (log + '_$(JID)', baseCommand)
                command += config[0]
                command += " output=" + options.output + noLumiTag + "_$(JID)"
                command += " numPEs=%s intLumi=%s %s=%s " \
                           % (options.numPEs, lum, options.seed, section)
                command += ' inputPrefix=%s' % directory
                command += " ".join(options.other)
                target.write('%s\n' % command)
    target.close()
