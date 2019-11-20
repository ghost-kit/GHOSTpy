import sys
import getopt

class lsArgs(object):

    def __init__(self, valid_args=None):
        assert isinstance(valid_args,dict)

        self.arg_dict = dict()
        self.command_called = sys.argv[0]
        self.args_called = sys.argv[1:]

        # build argument calls
        argstring = ""
        longargs = []
        argtuples = []

        for key in valid_args:
            argstring += "{}:".format(key)
            longargs.append("{}=".format(valid_args[key]))
            argtuples.append(tuple(["-{}".format(key), "--{}=".format(valid_args[key])]))

        # get arguments
        try:
            opts, args = getopt.getopt(self.args_called, argstring, longargs)
        except getopt.GetoptError:
            self.valid = False
            return

        # record valid arguments
        keys = valid_args.keys()
        for opt, arg in opts:
            ct = 0
            for cond in argtuples:
                if opt in cond:
                    self.arg_dict[valid_args[keys[ct]]] = arg
                ct += 1

        # record missed arguments
        for key in valid_args:
            if valid_args[key] not in self.arg_dict:
                self.arg_dict[valid_args[key]] = None





