import unittest
from base_template import *

class InputLoadingTestCase(unittest.TestCase):
    
    def test_float_value(self):
        user_input = load_user_input()
        self.assertTrue(float(user_input["Telescope"]["Binning"]) == -1)

    def test_bool_value(self):
        user_input = load_user_input()
        self.assertTrue(bool(user_input["Xsec"]["Cloud"]["Enable"]))


if __name__ == '__main__':
    unittest.main()