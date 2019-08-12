"""
This Module is the graphic user interface wrapper for accessing the all small molecules database 


Things to do,
1. clean up print statments
2. config files
3. comments
4. dependencies. reduce code blob.
5. more dynamic plotting options like add, remove specific molecule 
6. add spectra analysis tools like area under curve, etc
7. recreate spectra table in db_interact and remove local location hardlink

Code Ownership: Sara Seager, Zhuchang Zhan
"""

from rdkit import Chem
from rdkit.Chem import Draw

import numpy as np
import db_management as dbm 
import DIRfile as df
from Tkinter import *
from textwrap import wrap
from PIL import Image, ImageTk

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

IMPLEMENT = True    
VERBOSE = True
VERSION = 0.2

font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }

COLORS = ['r','b','g','c','m','y','k',"0.9","0.8","0.7","0.6","0.5","0.4","0.3","0.2","0.1"]

def data_process(filename,spliter=""):
    """
    process data files into numpy arrays that can be plotted.
    input files are arbitrary columns.
    returns a numpy array and the file title
    """
    ndata = []
    origin = []
    count = 0
    for line in open(filename, 'r'):
        try:
            if "#" in line:
                origin.append(line)
            if spliter == "":
                init = line.split()                           # acquire data from the txt file
            else:
                init = line.split(spliter)
            data = np.array(init, dtype=np.float,ndmin=2)  # convert data into numpy array           
            if ndata == []:     
                ndata = data
            else:
                ndata = np.concatenate((ndata,data))      # data stitching
            count+=1
        except ValueError:
            pass
    ndata = np.transpose(ndata)
    return ndata,origin


class StdoutRedirector(object):
    def __init__(self,text_widget):
        self.text_space = text_widget

    def write(self,string):
        self.text_space.insert('end', string)
        self.text_space.see('end')


class Example(Frame):
  
    def __init__(self, parent):
        Frame.__init__(self, parent)   
        
        db_name = df.DBDIR+"/molecule_db.db"
        self.c,self.conn = dbm.access_db(db_name)
        
        self.parent = parent    
        self.INCHIKEY = "VNWKTOKETHGBQD-UHFFFAOYSA-N"  # default set to CH4   
    
        data = self.c.execute("SELECT iupac_name,inchikey,inchi,formula from ID WHERE inchikey='%s'"%self.INCHIKEY)
        
        #self.INCHI = 
        self.FORMULA = "CH4"
        self.wave = "Wavenumber" 
        self.color = 0
        self.initUI()
        self.grid()
        
    def initUI(self):
        """
        creates the GUI and set all get widgets
        """
        
        # Top Center Software Title Frame
        self.frame0 = Frame(self).grid(row=0, column=0)
        self.beginlabel = Label(self.frame0, text="Spectra Search for All Small Molecules").grid(row=0, column=0)
        
        
        # Top Left search Frame
        self.frame1 = Frame(self)
        self.frame1.grid(row=1, column=0, rowspan=2,sticky=N)
        
        self.inchikeylable = Label(self.frame1, text="InChiKey").grid(row=0, column=0)
        self.Inchikey = StringVar()
        self.Inchikey.set(self.INCHIKEY)
        self.inchikey = Entry(self.frame1, textvariable=self.Inchikey).grid(row=0, column=1)
        
        self.inchilable = Label(self.frame1, text="InChi").grid(row=1, column=0)
        self.Inchi= StringVar()
        self.inchi= Entry(self.frame1, textvariable=self.Inchi).grid(row=1, column=1)
    
        self.formulalable = Label(self.frame1, text="Formula").grid(row=2, column=0)
        self.Formula= StringVar()
        self.formula= Entry(self.frame1, textvariable=self.Formula).grid(row=2, column=1)
    
        self.smilelable = Label(self.frame1, text="Smile").grid(row=3, column=0)
        self.Smile= StringVar()
        self.smile= Entry(self.frame1, textvariable=self.Smile).grid(row=3, column=1)
    
        self.commonlable = Label(self.frame1, text="Common Name").grid(row=4, column=0)
        self.Common= StringVar()
        self.common= Entry(self.frame1, textvariable=self.Common).grid(row=4, column=1)

        self.b0 = Button(self.frame1,text="Info",command=self.item_info).grid(row=5, column=0)           
        self.b00 = Button(self.frame1,text="Clear",command=self.clear_data).grid(row=5, column=1)      
      
        self.b1 = Button(self.frame1,text="Search",command=self.search_spectra).grid(row=6, column=0)
        self.b2 = Button(self.frame1,text="Display",command=self.display_spectra).grid(row=6, column=1)
        self.b3 = Button(self.frame1,text="Overlap",command=self.overlap_spectra).grid(row=6, column=2)  


        # Center Plot display frame
        self.frame2 = Frame(self)
        self.frame2.grid(row=2, column=1, rowspan = 3)
        
        fig = Figure(figsize=(7,5.5), dpi=100)
        fig.patch.set_facecolor('white')
        self.spectra_plot = fig.add_subplot(111)
       
        self.canvas = FigureCanvasTkAgg(fig, master=self.frame2)
        self.canvas.show()
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky=N) 
        
        """
        self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame2)
        self.toolbar.update()
        self.canvas._tkcanvas.grid(row=0, column=0, sticky=N)   
        """

        # Left Middle Top info messages and error messages Frame
        self.frame3 = Frame(self)
        self.frame3.grid(row=3, column=0, sticky=N)    
        self.signlable = Label(self.frame3, text="Search Results").grid(row=0, column=0, sticky=N)
        
        self.text_box = Text(self.frame3, wrap='word', height = 11, width=50)
        self.text_box.grid(row=1, column=0, sticky='NSWE', padx=4, pady=4)
        self.text_box.config(highlightbackground='black')
        #redirecting stdout and stderr to the result box
        sys.stdout = StdoutRedirector(self.text_box)
        sys.stderr = StdoutRedirector(self.text_box)


        # Left Middle Bottom Reference Frame
        self.frame4 = Frame(self)
        self.frame4.grid(row=4,column=0)
        self.signlable = Label(self.frame4, text="Molecule Reference").grid(row=0, column=0)
        self.select = Listbox(self.frame4, height= 11, width = 50)
        self.select.grid(row=1,column=0,columnspan=3)
        self.display_reference()
        self.ba = Button(self.frame4,text="Lookup",command=self.item_look_up).grid(row=2, column=0)
        self.bb = Button(self.frame4,text="Display",command=self.display_spectra).grid(row=2, column=1)
        self.bc = Button(self.frame4,text="Overlap",command=self.overlap_spectra).grid(row=2, column=2)  


        # Left Bottom signature and version frame
        self.frame5 = Frame(self)
        self.frame5.grid(row=5, column=0)    
        self.signlable = Label(self.frame5, text="Created and maintained by Zhuchang Zhan").grid(row=0, column=0)
        self.versionlable = Label(self.frame5, text="Version = %s"%VERSION).grid(row=1, column=0)


        # Center Top Plot function frame
        self.frame6 = Frame(self)
        self.frame6.grid(row=1,column=1)
        
        self.wavebutton = Button(self.frame6,text="Wavelength/Wave Number",comman=self.wave_switch).grid(row=0, column=0)
        self.clearbutton = Button(self.frame6,text="Clear Display",command=self.clear_display).grid(row=0, column=1) 
        self.quitbutton = Button(self.frame6,text="Quit Program",command=self._quit).grid(row=0, column=2)  

        # 2D structure frame
        
        """
        self.frame7 = Frame(self)
        self.frame7.grid(row=2, column=2)
        
        m = Chem.MolFromInchi('')
        Draw.MolToFile(m,"test.png")        
        
        self.img = Image.open('test.png')#df.DATADIR+'/%s/%s.png'%(self.INCHIKEY,self.INCHIKEY))
        self.tatras = ImageTk.PhotoImage(self.img)
        
        self.TwoDstructure = Canvas(self.frame7, width=self.img.size[0]+20, height=self.img.size[1]+20)
        self.TwoDstructure.create_rectangle(0, 0, self.img.size[0]+20, self.img.size[1]+20, outline="#f11", fill="#1f1", width=1)       
        self.image_on_canvas =  self.TwoDstructure.create_image(10, 10, anchor=NW, image=self.tatras,tags='2Dstruct')      
        self.TwoDstructure.pack()
        
        """

    def item_info(self):
        """
        Query the database for information regarding the designated molecule
        """

        if self.Inchikey.get() != "":
            key = "InChiKey"
            keyword = self.Inchikey.get()
        elif self.Inchi.get() != "":
            key = "Inchi"
            keyword = self.Inchi.get()
        elif self.Formula.get() != "":
            key = "Formula"
            keyword = self.Formula.get()
        elif self.Smile.get() != "":
            key = "Smiles"
            keyword = self.Smile.get()
        elif self.Common.get() !="":
            key = "Common_Name"
            keyword = self.Common.get()
        else:
            print "No Entry"     
            return   
        
        data = self.c.execute("SELECT InChiKey,InChi,Formula,Smiles,Common_Name,N,Molar Mass,Origin from ID where %s = '%s'"%(key,keyword))
        items = data.fetchall()
        for i in items:
            print "Keyword: ",keyword
            for j in i:
                print j
            print ""

    def clear_data(self):
        """
        Clears the search boxes
        """
        self.Inchikey.set("")
        self.Inchi.set("")
        self.Formula.set("")
        self.Smile.set("")
        self.Common.set("")


    def clear_display(self):
        """
        Clears the canvas and display an empty plot
        """
        self.spectra_plot.clear()
        self.spectra_plot.set_xlim([0,1])
        self.spectra_plot.set_ylim([0,1])
        self.canvas.draw()        

    def search_spectra(self):
        """
        Convert the data provided into inchikey and 
        look for the spectra in the database
        
        if entry is null, will not search
        but please note that currently only support one entry from box
        """
        self.noentry = False
        self.noinchi = False
        self.nospectra = False
        self.INCHIKEY = ""  #reset inchikey
        
        if self.Inchikey.get() != "":
            print self.Inchikey.get()
            command = "SELECT InChiKey,Formula from ID WHERE InchiKey = '%s'"%self.Inchikey.get()
        elif self.Inchi.get() != "":
            print self.Inchi.get()
            command = "SELECT InChiKey,Formula from ID WHERE Inchi = '%s'"%self.Inchi.get()
        elif self.Formula.get() != "":
            print self.Formula.get()
            command = "SELECT InChiKey,Formula from ID WHERE Formula = '%s'"%self.Formula.get()
        elif self.Smile.get() != "":
            print self.Smile.get()
            command = "SELECT InChiKey,Formula from ID WHERE Smile = '%s'"%self.Smile.get()
        elif self.Common.get() !="":
            print self.Common.get()
            command = "SELECT InChiKey,Formula from ID WHERE Common = '%s'"%self.Common.get()
        else:
            print "No Entry"
            self.noentry = True
            return
        
        inputdata = self.c.execute(command)
        results = inputdata.fetchone() #right now only dealing with single cases
        if results == None:
            print "Inchikey not exist for this input"
            self.noinchi = True
            return
        
        self.INCHIKEY = results[0]
        self.FORMULA = results[1]
        # search for the spectra
        command = "SELECT Local_Location,filename from Spectra WHERE InChiKey = '%s'"%self.INCHIKEY
        data = self.c.execute(command)
        item = data.fetchone() #right now only dealing with single cases
        if item == None:
            print "Spectra not exist for this Inchikey"
            self.nospectra = True
            return
        
        self.file_location = item[0]+"/"+item[1]
        print self.FORMULA,self.INCHIKEY,"is found."
        print "Location: ",self.file_location
        print "\n"*5
        
    def display_spectra(self):
        """
        display the spectra for given molecule, will reset the canvas so that 
        only one spectra is displayed.
        """
        self.titlename = []       # reset molecule
        self.spectra_plot.clear() # reset plot canvas
        self.color = 0            # reset color to default
        
        if self.noentry or self.noinchi or self.nospectra:
            return
        data,origin = data_process(self.file_location,"")
        for i in data[1:]:        
            #self.save_location = self.file_location.rstrip(".jdx")+"_%s"%count+".png"
            #simple_plot(data[0],i,False,True,self.save_location,"Wavenumber (1/CM)","Transmittance/Absorbance","Molecule Spectra Feature for %s"%self.INCHIKEY)
            x = data[0]
            y = i
            break #right now only dealing with single cases

        self.xdata = x
        self.ydata = y
        
        # reading what unit x and y is from the data file title
        XUNIT = ""
        YUNIT = ""
        for i in origin:
            if "XUNITS" in i:
                XUNIT = i
            if "YUNITS" in i:
                YUNIT = i
        if XUNIT == "" or YUNIT == "":
            print "Something is wrong with datafile, xunit or yunit is missing"
            
        # convert xdata to the right unit as set by self.wave
        if ("MICROMETERS" in XUNIT and self.wave == "Wavenumber") or ("CM" in XUNIT and self.wave == "Wavelength"):
            self.xdata = 10000.0/self.xdata    
        elif ("MICROMETERS" in XUNIT and self.wave == "Wavelength") or ("CM" in XUNIT and self.wave == "Wavenumber"):
            pass
        else:
            print "Something is wrong with xdata, please check file and settings"

        # converting absorbance to transmittance, not sure if method is correct. 
        # need some inputs for this. original file with absorbance range from 0-1000 ish.
        if "ABSORBANCE" in YUNIT:
            self.ydata = 1-self.ydata/float(max(self.ydata))

        # Configure plot title, x and y labels        
        self.titlename.append(self.FORMULA)
        self.spectra_plot.set_title("Molecule Spectra for %s"%self.FORMULA)
        if self.wave == "Wavenumber":
            self.spectra_plot.set_xlabel("Wavenumber (1/CM)")
        elif self.wave == "Wavelength":       
            self.spectra_plot.set_xlabel("WaveLength (micron)")
        else:
            print "Something is wrong with self.wave, please check settings"
        self.spectra_plot.set_ylabel("Transmittance") 
        self.spectra_plot.plot(self.xdata,self.ydata, linestyle="-", color=COLORS[self.color], marker="*",linewidth=1.0, markersize=2)        
        self.canvas.draw()
        self.display_2D_structure()

    def overlap_spectra(self):
        """
        plot new molecule onto the old ones
        """
        
        if self.noentry or self.noinchi or self.nospectra:
            return
        data,origin = data_process(self.file_location,"")
        count = 1
        for i in data[1:]:        
            x = data[0]
            y = i
            break #right now only dealing with single cases

        
        self.xdata = x
        self.ydata = y
        
        XUNIT = ""
        YUNIT = ""
        for i in origin:
            if "XUNITS" in i:
                XUNIT = i
            if "YUNITS" in i:
                YUNIT = i
        if XUNIT == "" or YUNIT == "":
            print "Something is wrong with datafile, xunit or yunit is missing"
            
        # convert xdata to the right unit as set by the self.wave
        if ("MICROMETERS" in XUNIT and self.wave == "Wavenumber") or ("CM" in XUNIT and self.wave == "Wavelength"):
            self.xdata = 10000.0/self.xdata    
        elif ("MICROMETERS" in XUNIT and self.wave == "Wavelength") or ("CM" in XUNIT and self.wave == "Wavenumber"):
            pass
        else:
            print "Something is wrong with xdata, please check file and settings"

        # converting absorbance to transmittance, not sure if method is correct. 
        # need some inputs for this. original file with absorbance range from 0-1000 ish.
        if "ABSORBANCE" in YUNIT:
            self.ydata = 1-self.ydata/float(max(self.ydata))
        
        # Configure plot title, x and y labels
        self.titlename.append(self.FORMULA)
        name = " "
        for i in self.titlename:
            name = "".join([name,i,","])
        self.spectra_plot.set_title("\n".join(wrap("Molecule Spectra for %s"%name[:-1], 60))) #no this is not a smiley face, it's to get rid of the last ","
        if self.wave == "Wavenumber":
            self.spectra_plot.set_xlabel("Wavenumber (1/CM)")
        elif self.wave == "Wavelength":       
            self.spectra_plot.set_xlabel("WaveLength (micron)")
        else:
            print "Something is wrong with self.wave, please check settings"
        self.spectra_plot.set_ylabel("Transmittance") 
        
        self.color += 1
        mark = "*"
        if self.color == len(COLORS):
            self.color = 0
            mark = "v"
        self.spectra_plot.plot(self.xdata,self.ydata, linestyle="-", color=COLORS[self.color], marker=mark,linewidth=1.0, markersize=3)        
        self.canvas.draw()
        self.display_2D_structure()
        

    def display_reference(self):
        """
        Query the database for molecules that have spectra and display them in the reference frame
        """
        data = self.c.execute("""SELECT ID.Formula, ID.InChiKey from ID,spectra 
                                 WHERE ID.InChiKey = spectra.InChikey AND spectra.Local_location NOT NULL 
                                 ORDER BY ID.Formula""")
        items = data.fetchall()
        self.inchikeylist = []
        for Formula, InChiKey in items:
            self.inchikeylist.append(InChiKey)
            format = Formula+" "*(20-len(Formula))
            
            self.select.insert (END, format+"\t"+InChiKey)

    def item_look_up(self):
        """
        lookup and search for a particular item
        """
        index = int(self.select.curselection()[0])
        print "InChiKey set to: ",self.inchikeylist[index]
        self.Inchikey.set(self.inchikeylist[index])
        self.search_spectra()

    def wave_switch(self):
        """
        Toggle between wavelength and wavenumber
        default is wavenumber since nist file are in wavenumbers.
        """
        if self.wave == "Wavenumber":
            self.wave = "Wavelength"
            self.display_spectra()
        elif self.wave == "Wavelength":
            self.wave = "Wavenumber"
            self.display_spectra()
        else:
            print "Something is wrong with the wave switch"

    def display_2D_structure(self):
        
        
        self.TwoDstructure.delete("2Dstruct")
        self.newimg = Image.open(df.DATADIR+'/%s/%s.png'%(self.INCHIKEY,self.INCHIKEY))
        self.newtatras = ImageTk.PhotoImage(self.newimg)      
        self.image_on_canvas =  self.TwoDstructure.create_image(10, 10, anchor=NW, image=self.newtatras,tags="2Dstruct")      
        self.TwoDstructure.grid(row=0, column=0)
    
        
    def _quit(self):
        """
        terminate program
        """
        self.quit()     # stops mainloop
        self.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate       



if __name__ == '__main__':
    root = Tk()
    root.wm_title("Spectra Search for All Small Molecules V%s"%VERSION)
    ex = Example(root)
    root.mainloop()  
    
    
    
    
"""
#Deprecated but may have future uses
 

        #update an image

"""