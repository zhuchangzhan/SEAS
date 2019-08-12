
"""
This Module is the graphic user interface wrapper for accessing the all small molecules database 


Current Problems:
    Compatibility issues and quick fix:
        Issues with upgraded version of matplotlib:
            for ubuntu:
                sudo apt-get install tk-dev libpng-dev libffi-dev dvipng texlive-latex-base
                pip uninstall matplotlib
                pip install matplotlib
            for mac:
                brew install tk-dev libpng-dev libffi-dev dvipng texlive-latex-base ? (need testing)
            or:
                Set MPL = False
        Issues with anaconda-Tkinter and importing the right PIL:
            Set PIL = False, will use MPL to plot.
        Issues with RDkit
            Set RDK = False, will use local generated images
    Please have either PIL or MPL set to true. Otherwise, fix your computer. 
        

0.3 update
    include dynamic plotting for spectra
    better looking GUI
    eliminate realpath 
    added Molecule class for cleaner data passing
    added Structure class to display molecular structure
    added jdx_file reader and eliminated redundancy


Things to do,
    legends for plot
    config files, comments
    add spectra analysis tools like area under curve, etc
    rescrape nist
    rework database, 
    recreate spectra table in db_interact 
    add a scroll bar to structure frame
    pip install packaging
    plot info
    button and method to save spectra image
    understand dpi
    image naming

Code Ownership: Sara Seager, Zhuchang Zhan
"""
VERSION = "0.3.3.1"

# System Default Packages
import os
import sys
import shutil
import numpy as np
from Tkinter import *
from textwrap import wrap



# Local Packages
import jdx
import db_management as dbm 
import db_interact as dbi

# External Package Controllers
RDK = False
PILL = False
MPL = False

# ploting packages
if PILL:
    import PIL
    from PIL import Image, ImageTk

import matplotlib
from matplotlib.figure import Figure

if MPL:
    matplotlib.use('TkAgg')
    from matplotlib import pylab as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
else:
    import matplotlib.pyplot as plt

# chemistry packages
if RDK:
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
        rdavl = True
    except:
        rdavl = False
        print "Warning: rdkit not properlly installed, structure acquired from local saves"
          



font = {'family' : 'serif',
        'color'  : 'darkred',
        'weight' : 'normal',
        'size'   : 16,
        }

COLORS = ['r','b','g','c','m','y','k',"0.9","0.8","0.7","0.6","0.5","0.4","0.3","0.2","0.1"]


def check_db(path):
    if os.path.exists("../database/molecule_db.db"):
        return True
    else:
        return dbi.import_data_local(path)

class StdoutRedirector(object):
    def __init__(self,text_widget):
        self.text_space = text_widget

    def write(self,string):
        self.text_space.insert('end', string)
        self.text_space.see('end')


class Molecule():
    
    def __init__(self,object,type,c):
        
        self.INCHIKEY   = ""  # default set to CH4   
        self.INCHI      = ""
        self.SMILES     = ""
        self.FORMULA    = ""
        self.IUPAC_NAME = ""
        self.CAS        = "" 
        self.has_spectra= ""
        self.path       = ""
        
        self.reference  = object
        self.ref_type   = type
        self.c          = c
        
        self.find_molecule()
    
    def __str__(self):
        cmd = "%s\n"*8
        return cmd%(self.INCHIKEY,self.INCHI,self.SMILES,self.FORMULA,self.IUPAC_NAME,self.CAS,self.has_spectra,self.path)
    
    def find_molecule(self):
        
        cmd = "SELECT ID.inchikey, ID.inchi, ID.smiles, ID.formula, ID.iupac_name, ID.cas,\
                      Spectra.has_spectra, Spectra.path \
               From ID,Spectra \
               WHERE ID.%s='%s' \
               AND ID.inchikey=Spectra.inchikey"%(self.ref_type,self.reference)
        result = self.c.execute(cmd)
        data = result.fetchall()[0]
        
        self.INCHIKEY    = data[0]
        self.INCHI       = data[1]
        self.SMILES      = data[2]
        self.FORMULA     = data[3]
        self.IUPAC_NAME  = data[4]
        self.CAS         = data[5]
        self.has_spectra = data[6]
        self.path        = data[7]
        
        
class Struct_GUI(Frame):

    def __init__(self, master, selected):
        Frame.__init__(self, master)
        
        self.selected_molecule = selected
        self.master = master
        self.initUI()


    def initUI(self):
        
        self.frame = Frame(self.master)
        self.frame.grid(row=0,column=0)
        self.quitButton = Button(self.frame, text = 'Quit', width = 25, command = self.close_windows)
        self.quitButton.grid(row=0,column=0)

        self.titleframe = Frame(self.master)
        self.titleframe.grid(row=1,column=0,columnspan=3)

        self.imageframe = Frame(self.master)
        self.imageframe.grid(row=2,column=0,columnspan=3)
        
        count, num = 0,0
        Formula = []
        for molecule in self.selected_molecule:
            
            smile = molecule.SMILES
            formula = molecule.FORMULA
            Formula.append(formula)
            
            if PILL:
                if RDK:
                    m = Chem.MolFromSmiles(smile)   #some mods here needed to config the image
                    Draw.MolToFile(m,"temp.png")  
                    back = Image.open("temp.png")
                    print "here"
                else:
                    try:
                        stored = "../2D_Structures/%s.png"%molecule.INCHIKEY
                    except:
                        stored = "temp.png"
                    
                    back = Image.open(stored)
                background = ImageTk.PhotoImage(back)
                label = Label(self.imageframe,anchor="center",image=background)
                label.image = background
                label.grid(column=count,row=num,sticky='EW')              
                
            elif MPL:
                if RDK:
                    m = Chem.MolFromSmiles(smile)   #some mods here needed to config the image
                    Draw.MolToFile(m,"temp.png") 
                    #image = mpimg.imread("temp.png") 
                    stored = "temp.png"
                else:
                    stored = "../2D_Structures/%s.png"%molecule.INCHIKEY

                image_load = plt.imread(stored)
                fig = plt.figure(figsize=(5,4))
                im = plt.imshow(image_load)
                
                # a tk.DrawingArea
                canvas = FigureCanvasTkAgg(fig, self.imageframe)
                canvas.show()
                canvas.get_tk_widget().grid(row=num,column=count)
            else:
                print "Do not understand plotting method"
                self.close_windows()
            
  
            
            label2 = Label(self.imageframe,text=formula,anchor="center")
            label2.grid(column=count,row=num+1,sticky='EW')   
          
            count +=1  
            if count%3 == 0:
                count = 0
                num += 2  
                
        self.titlelabel = Label(self.titleframe, text=",".join(Formula)).grid(row=0, column=0)

        
    def close_windows(self):
        self.master.destroy()

class Spectra_GUI(Frame):

    def __init__(self, master, selected):
        Frame.__init__(self, master)   
        
        self.selected_molecule = selected
        self.master = master
        self.initUI()
        
    def initUI(self):
        
        self.frame = Frame(self.master)
        self.frame.grid(row=0,column=0)
        self.quitButton = Button(self.frame, text = 'Quit', width = 25, command = self.close_windows)
        self.quitButton.grid(row=0,column=0)

        self.titleframe = Frame(self.master)
        self.titleframe.grid(row=1,column=0,columnspan=3)

        self.imageframe = Frame(self.master)
        self.imageframe.grid(row=2,column=0,columnspan=3)      
                
                    
class Spec_GUI(Frame):
  
    def __init__(self, master):
        Frame.__init__(self, master)   
        
        
        db_name = "../database/molecule_db.db"
        self.c,self.conn = dbm.access_db(db_name)
        
        reference = "VNWKTOKETHGBQD-UHFFFAOYSA-N"
        ref_type  = "inchikey"
        
        self.molecule = Molecule(reference,ref_type,self.c)
        
        self.yunit = "Transmittance"
        self.xunit = "Wavenumber" 
        
        
        self.spectra_in_db()
        self.selected = []
        
        self.initUI()
        self.grid()

    def initUI(self):

        
        # Top Left search Frame
        self.frame1 = Frame(self)
        self.frame1.grid(row=2, column=0, rowspan=2,sticky=N)
        
        self.inchikeylable = Label(self.frame1, text="InChiKey").grid(row=0, column=0)
        self.Inchikey = StringVar()
        self.Inchikey.set(self.molecule.INCHIKEY)
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
    
        self.iupaclable = Label(self.frame1, text="Iupac").grid(row=4, column=0)
        self.Iupac= StringVar()
        self.iupac= Entry(self.frame1, textvariable=self.Iupac).grid(row=4, column=1)        
        
        self.caslable = Label(self.frame1, text="cas").grid(row=5, column=0)
        self.Cas= StringVar()
        self.cas= Entry(self.frame1, textvariable=self.Cas).grid(row=5, column=1)  

        
        self.frame1a = Frame(self)
        self.frame1a.grid(row=3, column=0, rowspan=2,sticky=N)        
        self.ba = Button(self.frame1,text="Info ",command=self.item_info).grid(row=6, column=0)  
        self.bb = Button(self.frame1,text="Add  ",command=self.add_molecule).grid(row=6, column=1)  
        self.bc = Button(self.frame1,text="Clear",command=self.clear_data).grid(row=7, column=0)
  
        
        self.frame3 = Frame(self)
        self.frame3.grid(row=1, column=0,sticky="N")       # select of names
        self.select = Listbox(self.frame3, width=30, height=20)
        self.select.grid(row=1,column=0)   
        self.setSelect()


        self.frameb = Frame(self)
        self.frameb.grid(row=1,column=1)
        self.toplable = Button(self.frameb, text=" T ",command=self.top_molecule).grid(row=0, column=0)
        self.uplable = Button(self.frameb, text=" ^ ",command=self.up_molecule).grid(row=1, column=0)    
        self.downlable = Button(self.frameb, text=" v ",command=self.down_molecule).grid(row=2, column=0)
        self.bottomlable = Button(self.frameb, text=" B ",command=self.bottom_molecule).grid(row=3, column=0) 
        self.loadlable = Button(self.frameb, text="->",command=self.load_molecule).grid(row=4, column=0)    
        self.deletelable = Button(self.frameb, text="<-",command=self.delete_molecule).grid(row=5, column=0)    
        self.deletealllable = Button(self.frameb, text="<<",command=self.delete_all_molecule).grid(row=6, column=0)    
        

        self.frame3a = Frame(self)
        self.frame3a.grid(row=1, column=2,sticky="N")       # select of names
        self.select2 = Listbox(self.frame3a, width=30, height=20)
        self.select2.grid(row=1,column=0,columnspan=2)      
        self.spectrabutton   = Button(self.frame3a, text="Spectra",command=self.display_spectra).grid(row=2, column=0)    
        self.structurebutton = Button(self.frame3a, text="Structure",command=self.display_structure).grid(row=2, column=1)       
        self.spectrasbutton  = Button(self.frame3a, text="Spectras",command=self.displaying_spectra).grid(row=2, column=2)  


        self.frame2 = Frame(self)
        self.frame2.grid(row=1, column=3, rowspan=2,sticky="N")   
        
        #MPL is the default method for displaying spectra
        if MPL:     
            self.fig = Figure(figsize=(6.5,5.5), dpi=100)
            self.fig.patch.set_facecolor('white')
            self.spectra_plot = self.fig.add_subplot(111)
            
            self.canvas = FigureCanvasTkAgg(self.fig, master=self.frame2)
            self.canvas.show()
            self.canvas.get_tk_widget().grid(row=1, column=0, sticky=N)  
            
        elif PILL:
            self.fig = plt.figure(figsize=(6.5,5.5), dpi=100)
            self.fig.patch.set_facecolor('white')
            self.spectra_plot = self.fig.add_subplot(111)
            
            back = Image.open("blank.png")
            background = ImageTk.PhotoImage(back)
            self.spectralabel = Label(self.frame2,anchor="center",image=background)
            self.spectralabel.image = background
            self.spectralabel.grid(row=1, column=0, sticky=N)        
        else:
            print "Unknown plotting method"
            
        
        
        #self.toolbar = NavigationToolbar2TkAgg(self.canvas, self.frame2)
        #self.toolbar.update()
        #self.canvas._tkcanvas.grid(row=0, column=0, sticky=N)                 
        
          
          
        self.frame_debug = Frame(self)
        self.frame_debug.grid(row=2, column=2, sticky="N")    
        #self.signlable = Label(self.frame_debug, text="Search Results").grid(row=0, column=0, sticky="NW")
        self.text_box = Text(self.frame_debug, wrap='word', height = 12, width=30)
        self.text_box.grid(row=1, column=0, sticky='NW')
        self.text_box.config(highlightbackground='black')
        #redirecting stdout and stderr to the result box
        sys.stdout = StdoutRedirector(self.text_box)
        sys.stderr = StdoutRedirector(self.text_box)


        self.frame6 = Frame(self)
        self.frame6.grid(row=0,column=3)
        self.xbutton = Button(self.frame6,text="Wavelength/Wave Number",comman=self.x_switch).grid(row=0, column=0)
        self.ybutton = Button(self.frame6,text="Trans./Absorb.",comman=self.y_switch).grid(row=0, column=1)
        self.clearbutton = Button(self.frame6,text="Clear Display",command=self.clear_display).grid(row=0, column=2) 
        self.clearbutton = Button(self.frame6,text="Save Display",command=self.save_display).grid(row=0, column=3)
        self.quitbutton = Button(self.frame6,text="Quit Program",command=self._quit).grid(row=0, column=4)  



    def x_switch(self):
        """
        Toggle between wavelength and wavenumber
        default is wavenumber since nist file are in wavenumbers.
        """
        if self.xunit == "Wavenumber":
            self.xunit = "Wavelength"
            self.display_spectra()
        elif self.xunit == "Wavelength":
            self.xunit = "Wavenumber"
            self.display_spectra()
        else:
            print "Something is wrong with the wave switch"    

    def y_switch(self):
        """
        Toggle between transmittance and absorbance
        default is transmittance 
        """
        if self.yunit == "Transmittance":
            self.yunit = "Absorbance"
            self.display_spectra()
        elif self.yunit == "Absorbance":
            self.yunit = "Transmittance"
            self.display_spectra()
        else:
            print "Something is wrong with the wave switch"            

    def item_info(self):
        """
        Query the database for information regarding the designated molecule
        """

        if self.Inchikey.get() != "":
            type = "inchikey"
            keyword = self.Inchikey.get()
        elif self.Inchi.get() != "":
            type = "inchi"
            keyword = self.Inchi.get()
        elif self.Formula.get() != "":
            type = "formula"
            keyword = self.Formula.get()
        elif self.Smile.get() != "":
            type = "smiles"
            keyword = self.Smile.get()
        elif self.Iupac.get() !="":
            type = "iupac_name"
            keyword = self.Iupac.get()
        elif self.Cas.get() !="":
            type = "cas"
            keyword = self.Cas.get()
        else:
            print "No Entry"     
            return   
        
        self.molecule = Molecule(keyword,type,self.c)
        print self.molecule  
                
    def display_structure(self):
        """
        Open up the tab for structure
        """
        self.Structure_Window = Toplevel(self.master)
        self.app = Struct_GUI(self.Structure_Window,self.selected)  
    
    def displaying_spectra(self):
        """
        Open up the tab or spectra
        """
        self.Spectra_Window = Toplevel(self.master)
        self.app = Spectra_GUI(self.Spectra_Window,self.selected)  
        
    def spectra_in_db(self):
        
        temp= self.c.execute('SELECT ID.inchikey FROM ID,Spectra \
                         WHERE ID.inchikey=Spectra.inchikey AND Spectra.has_spectra="Y"')  
        inchikey_with_spectra = []
        for i in temp.fetchall():
            inchikey_with_spectra.append(i[0])
        
        self.result = []
        type = "inchikey"
        for inchikey in inchikey_with_spectra:
            molecule = Molecule(inchikey,type,self.c)
            self.result.append(molecule)
                
    def display_spectra(self):
        
        self.spectra_plot.clear() # reset plot canvas
        self.xdata = []
        self.ydata = []
        filename = []
        titlename = []
        
        # deeper look into here needed
        for i in self.selected:
            path = i.path
            for j in os.listdir(path):
                if "jdx" in j:
                    print path
                    titlename.append(j)
                    filename.append(path+"/"+j)  #how to handle multiple file?
        
        for file in filename:
            data = jdx.JdxFile(file)
            
            if self.xunit == "Wavenumber":
                x = data.wn()
            elif self.xunit == "Wavelength":
                x = data.wl()
                
            if self.yunit == "Transmittance":
                y = data.trans()
            elif self.yunit == "Absorbance":
                y = data.absorb()

            self.xdata.append(x)
            self.ydata.append(y)  
            
            
        name = " "
        for i in titlename:
            name = "".join([name,i,","])
        self.spectra_plot.set_title("\n".join(wrap("Molecule Spectra for %s"%name[:-1], 60))) #no this is not a smiley face, it's to get rid of the last ","


        if self.xunit == "Wavenumber":
            self.spectra_plot.set_xlabel("Wavenumber (1/CM)")
        elif self.xunit == "Wavelength":       
            self.spectra_plot.set_xlabel("WaveLength (micron)")
        else:
            print "Something is wrong with self.xunit, please check settings"

        if self.yunit == "Transmittance":
            self.spectra_plot.set_ylabel("Transmittance") 
        elif self.yunit == "Absorbance":
            self.spectra_plot.set_ylabel("Absorbance") 
        
        
        if MPL:
            for i in range(len(self.selected)):
                self.spectra_plot.plot(self.xdata[i],self.ydata[i], linestyle="-", color=COLORS[i], marker="*",linewidth=1.0, markersize=1)
                self.canvas.draw()
        elif PILL:
            for i in range(len(self.selected)):
                self.spectra_plot.plot(self.xdata[i],self.ydata[i], linestyle="-", color=COLORS[i], marker="*",linewidth=1.0, markersize=1)
            self.fig.savefig("spectra_temp.png")
            img2 = ImageTk.PhotoImage(Image.open("spectra_temp.png"))
            self.spectralabel.configure(image = img2)
            self.spectralabel.image = img2           
        
       

    def save_display(self):
        savename = "plots/"+str(self.spectra_plot.title)+".png"
        print "Spectra saved to %s"%savename
        
        if MPL:
            self.fig.savefig(savename)
        elif PILL:
            shutil.copy("spectra_temp.png",savename)
            

    def clear_data(self):
        """
        Clears the search boxes
        """
        self.Inchikey.set("")
        self.Inchi.set("")
        self.Formula.set("")
        self.Smile.set("")
        self.Iupac.set("")
        self.Cas.set("")        
  
    def clear_display(self):
        """
        Clears the canvas and display an empty plot
        """
        self.spectra_plot.clear()
        self.spectra_plot.set_xlim([0,1])
        self.spectra_plot.set_ylim([0,1])
        self.canvas.draw()               
            
    def top_molecule(self):
        selection = self.whichSelected2()
        
        temp = [self.selected[selection]]
        del self.selected[selection]
        temp.extend(self.selected)
        self.selected = temp
        
        self.setSelect2()
        
    def up_molecule(self):
        selection = self.whichSelected2()
        if selection == 0:
            return
        else:
            self.selected[selection],self.selected[selection-1] = self.selected[selection-1],self.selected[selection]
        self.setSelect2()
    
    def down_molecule(self):
        selection = self.whichSelected2()
        if selection == len(self.selected):
            return
        else:
            self.selected[selection],self.selected[selection+1] = self.selected[selection+1],self.selected[selection]
        self.setSelect2()
    
    def bottom_molecule(self):
        selection = self.whichSelected2()
        
        temp = [self.selected[selection]]
        del self.selected[selection]
        self.selected.extend(temp)
        
        self.setSelect2()

    def load_molecule(self):
        i = self.result[self.whichSelected()]
        self.selected.append(i)
        self.setSelect2()
        #self.selected.sort()
    
    def delete_molecule(self):
        
        selection = self.whichSelected2()
        del self.selected[selection]
        self.setSelect2()
    
    def delete_all_molecule(self):
        self.selected = []
        self.setSelect2()

    def add_molecule(self):
        
        self.item_info()
        
        if self.molecule.has_spectra == "Y":
            #item = [self.FORMULA,self.INCHIKEY,path,self.SMILES]
            self.selected.append(self.molecule)
            self.setSelect2()
        else:
            print "molecule does not have spectra, not added"
    
    def setSelect(self):
        for i in self.result:
            name = i.FORMULA+" "*(20-len(i.FORMULA))+"\t"+i.INCHIKEY
            self.select.insert(END, name)

    def setSelect2(self):
        self.select2.delete(0,END)
        for i in self.selected:
            name = i.FORMULA+" "*(20-len(i.FORMULA))+"\t"+i.INCHIKEY
            self.select2.insert (END, name)
    
    def whichSelected(self):
        try:
            return int(self.select.curselection()[0])
        except IndexError:
            print "None Selected"
            return
         
    def whichSelected2(self):
        return int(self.select2.curselection()[0])

    def _quit(self):
        """
        terminate program
        """
        self.quit()     # stops mainloop
        self.destroy()  # this is necessary on Windows to prevent
                        # Fatal Python Error: PyEval_RestoreThread: NULL tstate

        
if __name__ == '__main__':
    
    path = "../data"
    if check_db(path) == False:
        sys.exit()
    
    root = Tk()
    root.wm_title("Spectra Search for All Small Molecules V%s"%VERSION)
    ex = Spec_GUI(root)
    root.mainloop() 







