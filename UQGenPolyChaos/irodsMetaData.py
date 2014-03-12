import sys
import os
import pickle
import ConfigParser

# setup command line parser to control execution
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--root_dir", 
                  action="store", dest="root_dir", default=None,
                  help="root directory for collections", metavar="DIR")
parser.add_option("--data_dir",
                  action="store", dest="data_dir", default=None,
                  help="data directory to add the meta data", metavar="DIR")
(options, args) = parser.parse_args()
if (options.data_dir != None and options.root_dir != None):
  # set a dictionary for passing to fem via Python kwargs
  # read in params previously stored in dictionary and written
  pkl_file = open('%s/CaseInfo.pkl' % options.data_dir , 'rb')
  fem_params = pickle.load(pkl_file)
  pkl_file.close()
  
  fem_options = fem_params['options'] 
  fem_config  = fem_params['config_parser']
  # copy all values from the input file
  for section in fem_config.sections():
    for name,value in  fem_config.items(section):
      os.system( "imeta add -C %s/%s '%s/%s' '%s' " % (options.root_dir,options.data_dir,section,name,value))
      #getpot.SetIniValue( "%s/%s" % (section,name) , value ) 
  
  for key,value in  fem_options.__dict__.items():
    os.system("imeta add -C %s/%s 'options/%s' '%s'" % (options.root_dir,options.data_dir,key,value))
else:
  parser.print_help()
  print options
