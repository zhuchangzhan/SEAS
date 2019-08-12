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
A timer code to help count time

"""
import time


class simple_timer():
    
    def __init__(self, precision=2):
        
        self.precision = str(precision)
        
        self.format = "".join(["%.",self.precision,"f"])
        
        
        self.current_time = time.time()
        self.initial_time = self.current_time
        self.counter_time = self.current_time
        
    def update(self): 
        
        self.counter_time = time.time()
        
    def elapse(self):
        "time for each lapse"
        
        self.current_time = time.time()
        result = self.format%(self.current_time - self.counter_time)
        self.update()
        
        return result

    def progress(self,count,total):
        "Estimated Time of completion"
        
        self.update()
        self.current_time = time.time()
        
        percent = float(count)/float(total)
        
        time_passed = self.current_time - self.initial_time
        time_passed_precision = ("%"+".%sf"%(int(self.precision)/2))%time_passed
        
        
        try:
            estimated_time_remaining = ((1.-percent)/percent)*time_passed
        except ZeroDivisionError:
            estimated_time_remaining = "???"
            
        
        
        percent_precision = ("%"+".%sf"%(int(self.precision)/2))%(percent*100)
        try:
            remaining_precision = ("%"+".%sf"%(int(self.precision)/2))%(estimated_time_remaining)
            hour_precision = ("%"+".%sf"%(int(self.precision)/2))%(estimated_time_remaining/3600.)
            
        except:
            remaining_precision = estimated_time_remaining
            hour_precision = estimated_time_remaining
            
        result = " %ss Passed, %s Completed, Estimated %ss (%shr) Remaining"%(time_passed_precision, percent_precision, remaining_precision,hour_precision)
        
        return result

    def total(self):
        "total time"
        
        self.current_time = time.time()
        return self.format%(self.current_time - self.initial_time)
    


        
        
        
        