#!/usr/bin/python
#@Created by Jose Fernandez Navarrro <jose.fernandez.navarro@scilifelab.se>

import ConfigParser
import argparse
import os
import sys

class Config(object):
    """
    Reads the command line and the  config file.
    Merges options from both sources, giving priority to the command line.
    """
    
    class Args(argparse.Namespace):
        def __init__(self):
            self.properties = {}

    class SetProperty(argparse.Action):
        """
        Used with argparse to parse arguments setting property values.
        Creates an attribute 'property' in the results namespace containing
        all the property-value pairs read from the command line.
        """
        def __call__(self, parser, namespace, value, option_string=None):
            name, v = value.split('=', 1)
            namespace.properties[name] = v

    def __init__(self):
        
        self.cmd_parser = argparse.ArgumentParser(description='Distributed ST data analysis pipeline.')
        # make the parser print help whenever there's a usage error
        def error(message):
            sys.stderr.write('error: %s\n\n' % message)
            self.cmd_parser.print_help()
            sys.exit(2)

        self.cmd_parser.error = error
        ##############

        self.cmd_parser.add_argument('input', metavar='INPUT', help='input path')
        self.cmd_parser.add_argument('output', metavar='OUTPUT', help='output path')
        self.cmd_parser.add_argument('-r', '--num-reducers', metavar='INT', type=int, dest="num_reducers",
                help="Number of reduce tasks. Specify 0 to perform a single process analysis (default: 3 * num task trackers).")
        self.cmd_parser.add_argument('-cf', '--config-file', metavar='FILE', dest="config_file", default=os.path.join(os.path.expanduser('~'), '.cfg'),
                help='Override the default config file')
        self.cmd_parser.add_argument('-D', metavar="PROP=VALUE", action=type(self).SetProperty,
                help='Set a property value, such as -D mapred.compress.map.output=true')
        self.cmd_parser.add_argument('--reference', type=argparse.FileType('w'), dest="reference",
                help="Reference genome file to use to map the reads")
        
    def load_config_and_cmd_line(self, argv=sys.argv[1:]):
        # we scan the command line first in case the user wants to
        # override the default config file location
        args, left_over = self.cmd_parser.parse_known_args(args=argv, namespace=Config.Args())

        # load the config for this program, if the file exists
        config = ConfigParser.ConfigParser()

        # was a config file different from the default specified on the command line?
        try:
            if args.config_file != self.cmd_parser.get_default("config_file"):
                # in this case, make sure it exists and is readable
                if not os.path.exists(args.config_file):
                    print >>sys.stderr, "WARNING: The specified Seal config file %s doens't exist" % args.seal_config
                if not os.access(args.config_file, os.R_OK):
                    print >>sys.stderr, "WARNING: The specified Seal config file %s isn't readable" % args.seal_config
                config.read(args.config_file) # no problems.  Load the file.
            else:
                # check the default path.  If the file exists and is readable we'll load it
                if os.path.exists(args.config_file):
                    if os.access(args.config_file, os.R_OK):
                        config.read(args.config_file)
                    else:
                        print >>sys.stderr, "WARNING:  config file %s exists but isn't readable" % args.seal_config
                        
        except ConfigParser.Error as e: # catch errors from parsing the config file
            print >>sys.stderr, "Error in configuration file %s\n%s" % (args.seal_config, str(e))

        # override configuration properties from file with the ones
        # provided on the command line.
        for name, value in config.items(section="DEFAULT"):
            if not args.properties.has_key(name):
                args.properties[name] = value

        return args, left_over
