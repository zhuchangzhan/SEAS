"""

Here is a construction site for building SEAS.

Will have to clean up once finished

Let's build the user input section first

"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.Common_Utils.configurable as config


import pprint
pp = pprint.PrettyPrinter(indent=4)




def main():
    
    user_input = config.Configuration("../../config/user_input.cfg")
    
    pp.pprint(user_input)
    
    
if __name__ == "__main__":
    main()
















