#!/usr/bin/env python
#
# Copyright (C) 2017 - Massachusetts Institute of Technology (MIT)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
functions related to interpolation

"""

import sys
from scipy import interpolate




def interpolate1d(x,y,X,style="linear"):
    """
    Linear interpolation only for now.
    """
    
    if style == "linear":
        interp = interpolate.interp1d(x,y)
        try:
            Y = interp(X)
        except:
            print(X)
            print(x)
            sys.exit()
        return Y
    else:
        print "other interpolation methods not implemented yet"
        sys.exit()

def interpolate2d(x,y,X,Y,Z):
    """
    Used to need it but no longer necessary
    Here is just a placeholder for future implementation if needed again
    """
    return ""



