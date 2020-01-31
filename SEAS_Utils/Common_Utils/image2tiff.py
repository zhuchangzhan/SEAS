"""

Utility module to convert png and jpeg image to tiff.

Astrobiology require image format to be TIF/EPS


# potentially add argparse for terminal friendly use
"""

import os
from PIL import Image



def image2tiff(srcname,srcdir="",outname="",outdir=""):
    
    if srcdir == "":
        img = Image.open(srcname)
        
        source_name = srcname.split("/")[-1].split(".")[0]       
        if outname == "":
            
            img.save(os.path.join(outdir,'%s.tiff'%source_name),dpi=(300, 300))


def check_img_dpi(srcname,srcdir = ""):
    
    if srcdir == "":
        img = Image.open(srcname)
        print(img.info['dpi'])