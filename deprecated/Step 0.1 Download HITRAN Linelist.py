"""
This might be deprecated because it seems to slow down calculation significantly
"""

import os
import sys

DIR = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(DIR, '../..'))

import SEAS_Utils.Web_Utils.web_downloader as wd


def download_HITRAN_Linelists():

    ID = wd.get_HITRAN_ID()
    
    HITRAN_Lines = "../../SEAS_Input/Line_List/HITRAN_Line_List/"
    numin = 0
    numax = 50000
    
    for i in ID:
        wd.HITRAN_Line_List_downloader(HITRAN_Lines,i,numin,numax,True,False)
        




if __name__ == "__main__":
    download_HITRAN_Linelists()