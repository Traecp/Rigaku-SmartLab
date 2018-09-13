# -*- coding: utf-8 -*-
"""How to use this script:
Open a command prompt window
Run: ipython --pylab
Type: %run SmartLab.py IP-RSM/IPRSM*.ras
"""
from __future__ import division
import os, sys, time, glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#
#mpl.rcParams['font.size'] = 18.0
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.handletextpad'] = 0.5
mpl.rcParams['legend.fontsize'] = 'medium'
mpl.rcParams['figure.subplot.bottom'] = 0.13
mpl.rcParams['figure.subplot.top'] = 0.93
mpl.rcParams['figure.subplot.left'] = 0.14
mpl.rcParams['figure.subplot.right'] = 0.915
mpl.rcParams['image.cmap'] = "jet"
mpl.rcParams['savefig.dpi'] = 300

class RasFile(object):
    RAS_HEADER_START = "*RAS_HEADER_START"
    RAS_HEADER_END   = "*RAS_HEADER_END"
    RAS_INT_START    = "*RAS_INT_START"
    RAS_INT_END      = "*RAS_INT_END"
    HEADER_SPLIT     = "\""
    DATA_SPLIT       = " "
    def __init__(self, filePath):
        self.filename = filePath
        self.header = dict()
        self.scan_data = []
        self.axis_names  = dict() #{0: "TwoTheta", 1: "Omega", ...} use internal name
        self.axis_values = dict() # {0:[1,2,..], 1:[45,20,000],...}
        self.MEAS_COND_AXIS_NAME = dict() #Array
        self.MEAS_COND_AXIS_NAME_INTERNAL = dict() #Array
        self.MEAS_COND_AXIS_OFFSET = dict() #Array
        self.MEAS_COND_AXIS_POSITION = dict() #Array
        self.MEAS_COND_AXIS_UNIT = dict() #Once
        
        self.scan_axis = ""
        self.scan_axis_internal = ""
        self.scan_angles = []
        self.number_of_scan = 0
        self.points_per_scan = 0
        self.parse()
    def parse(self):
        with open(self.filename, encoding="Latin-1", mode="r") as f:
            scan_start   = False
            scan_end     = False
            header_start = False
            scan_data    = []
            scan_angle   = []
            header_initialized = False
            for line in f:
                if line.strip():
                    line = line.strip()
                    # print("Scan start: ", scan_start)
                    if line.startswith(self.RAS_HEADER_START):
                        header_start = True
                        # print(line)
                        continue
                    if line.startswith(self.RAS_HEADER_END):
                        header_start = False
                        header_initialized = True
                        # print(line)
                        continue
                    if line.startswith(self.RAS_INT_START):
                        scan_start = True
                        # scan_end   = False
                        # print(line)
                        continue
                    if line.startswith(self.RAS_INT_END):
                        scan_start = False
                        # scan_end   = True
                        # print(line)
                        pad_points = self.points_per_scan - len(scan_data)
                        if pad_points >0:
                            print("Data not complete. Number of data point missing for this scan: ", pad_points)
                            pad_data   = [0]*pad_points
                            scan_data.extend(pad_data)
                        self.scan_data.append(scan_data)
                        self.number_of_scan +=1
                        scan_data = []
                        scan_angle= []
                        # continue
                        
                    if scan_start:
                        ls = line.split(self.DATA_SPLIT)
                        # print(ls)
                    elif header_start:
                        ls = line.split(self.HEADER_SPLIT)
                    else:
                        continue
                        
                        
                    if header_start:
                        key = ls[0][1:].strip()
                        val = ls[1].strip()
                        if not header_initialized: #If the header is read for the first time, we need to fill different metadata information (basically all)
                            self.header[key] = val #We collect all metadata in the header - done only Once.
                            if "MEAS_COND_AXIS_NAME-" in key:
                                tmp = key.split("-")
                                order = int(tmp[1].strip())
                                self.MEAS_COND_AXIS_NAME[order] = val
                            if "MEAS_COND_AXIS_NAME_INTERNAL-" in key:
                                tmp = key.split("-")
                                order = int(tmp[1].strip())
                                self.MEAS_COND_AXIS_NAME_INTERNAL[order] = val
                            if "MEAS_COND_AXIS_OFFSET-" in key:
                                tmp = key.split("-")
                                order = int(tmp[1].strip())
                                try:
                                    val = float(val)
                                except:
                                    val = 0
                                self.MEAS_COND_AXIS_OFFSET[order] = val
                            if "MEAS_COND_AXIS_POSITION-" in key:
                                tmp = key.split("-")
                                order = int(tmp[1].strip())
                                try:
                                    val = float(val)
                                    self.MEAS_COND_AXIS_POSITION[order] = [val]
                                except:
                                    self.MEAS_COND_AXIS_POSITION[order] = val
                            if "MEAS_COND_AXIS_UNIT-" in key:
                                tmp = key.split("-")
                                order = int(tmp[1].strip())
                                self.MEAS_COND_AXIS_UNIT[order] = val
                            if "MEAS_DATA_COUNT" in key:
                                self.points_per_scan = int(float(val))
                            if key == "MEAS_SCAN_AXIS_X":
                                self.scan_axis = val
                            if key == "MEAS_SCAN_AXIS_X_INTERNAL":
                                self.scan_axis_internal = val
                            if key == "MEAS_SCAN_START":
                                self.scan_angle_start = float(val)
                            if key == "MEAS_SCAN_STEP":
                                self.scan_angle_step = float(val)
                            if key == "MEAS_SCAN_STOP":
                                self.scan_angle_stop = float(val)
                                
                                
                        else: #Header already initialized, we add new position to the axis, if they are number and not string.
                            if "MEAS_COND_AXIS_POSITION-" in key:
                                tmp = key.split("-")
                                order = int(tmp[1].strip())
                                try:
                                    val = float(val)
                                    self.MEAS_COND_AXIS_POSITION[order].append(val)
                                except:
                                    continue
                                        
                    if scan_start:
                        a = float(ls[0].strip())
                        v = float(ls[1].strip())
                        scan_angle.append(a)
                        scan_data.append(v)
                        # print("Angle {:.2f} Intensity: {:.2f}".format(a,v))
                    
            self.scan_data = np.asarray(self.scan_data)
            if self.number_of_scan == 1:
                self.scan_data = self.scan_data[0]
            self.scan_angles = np.linspace(self.scan_angle_start, self.scan_angle_stop, self.points_per_scan)
            self.axis_names = self.MEAS_COND_AXIS_NAME
            self.axis_values= self.MEAS_COND_AXIS_POSITION
            print ("Parsing done. ", self.filename)
    
class PoleFigure:
    def __init__(self, filename):
        self.filename = filename
        self.ras      = RasFile(self.filename)
        
    def plot(self, vmin=None, vmax=None, log=True, save=False):
        self.data  = self.ras.scan_data    #2D grid
        self.phi   = self.ras.scan_angles  #2D grid
        for k in self.ras.axis_names.keys():
            if self.ras.axis_names[k]=="Chi":
                chi = self.ras.axis_values[k]
                self.chi = np.array(chi)
            if self.ras.axis_names[k]=="2-Theta":
                self.tth = self.ras.axis_values[k][0]
                
        self.phi   = np.radians(self.phi)
        if vmin==None:
            vmin = self.data.min()
            if log:
                vmin  = np.log10(vmin)
        if vmax==None:
            # vmax  = np.sort(self.data.ravel())[int(0.99*self.data.size)-1]
            vmax = self.data.max()*0.9
            if log:
                vmax = np.log10(vmax)
        self.fig = plt.figure()
        self.fig.subplots_adjust(top=0.9, bottom=0.2, left=0.2,right=0.95,hspace=0.2,wspace=0.2)
        self.ax  = self.fig.add_subplot(111, projection="polar")
        if log:
            self.img  = self.ax.contourf(self.phi, self.chi, np.log10(self.data), 50, vmin=vmin, vmax=vmax)
            self.cb = plt.colorbar(self.img, label="log10(Intensity)", format="%.2f")
        else:
            self.img  = self.ax.contourf(self.phi, self.chi, self.data, 50, vmin=vmin, vmax=vmax)
            self.cb = plt.colorbar(self.img, label="Intensity", format="%.2f")
        self.ax.set_title("2 Theta = %.4f deg."%self.tth)
        
        self.out_fig = os.path.splitext(self.ras.filename)[0]+".png"
        
        plt.show()
        if save:
            self.fig.savefig(self.out_fig)
        
   
class RSM:
    def __init__(self, filenames):
        self.filenames = filenames
        if len(self.filenames)>1:
            self.is_multiple_file_RSM = True
        else:
            self.is_multiple_file_RSM = False
            self.filenames = self.filenames[0]
            
        self.get_data()
        
    def get_data(self):
        if self.is_multiple_file_RSM:
            self.filenames.sort()
            data = []
            
            fn = self.filenames[0]
            self.filename_template = fn
            self.ras= RasFile(fn)
            self.all_axis_position = self.ras.MEAS_COND_AXIS_POSITION
            data.append(self.ras.scan_data)
            
            for f in self.filenames[1:]:
                raw = RasFile(f)
                data.append(raw.scan_data)
                all_pos = raw.MEAS_COND_AXIS_POSITION
                for k in all_pos.keys():
                    if type(all_pos[k])==list:
                        self.all_axis_position[k].extend(all_pos[k])
            self.first_axis_angle = self.ras.scan_angles
            self.first_axis_name  = self.ras.scan_axis
            self.data = np.asarray(data)
            # Which axis is the second scanning axis? Let's check
            for k in self.all_axis_position.keys():
                if type(self.all_axis_position[k])==list:
                    self.second_axis_angle = np.array(self.all_axis_position[k])
                    tmp = self.second_axis_angle - self.second_axis_angle[0]
                    if not np.allclose(tmp, 0):
                        self.second_axis_name = self.ras.MEAS_COND_AXIS_NAME[k]
                        break
                        
        else:
            self.filename_template = self.filenames
            self.ras = RasFile(self.filenames)
            self.data = self.ras.scan_data
            self.first_axis_name = self.ras.scan_axis
            self.first_axis_angle= self.ras.scan_angles
            self.all_axis_position = self.ras.MEAS_COND_AXIS_POSITION
            # Which axis is the second scanning axis? Let's check
            for k in self.all_axis_position.keys():
                if type(self.all_axis_position[k])==list:
                    self.second_axis_angle = np.array(self.all_axis_position[k])
                    tmp = self.second_axis_angle - self.second_axis_angle[0]
                    if not np.allclose(tmp, 0):
                        self.second_axis_name = self.ras.MEAS_COND_AXIS_NAME[k]
                        break
                        
    def plot_2D(self, vmin=None, vmax=None, log=True, save=False):
        self.fig,self.ax = plt.subplots(1,1)
        self.fig.subplots_adjust(top=0.9, bottom=0.2, left=0.2,right=0.95,hspace=0.2,wspace=0.2)
        if vmin==None:
            vmin=0
        if vmax==None:
            # vmax = np.sort(self.data.ravel())[int(0.99*self.data.size)-1]
            vmax = self.data.max()*0.9
            if log:
                vmax = np.log10(vmax)
        if log:
            self.img = self.ax.contourf(self.first_axis_angle, self.second_axis_angle, np.log10(self.data+1e-6), 50, vmin=vmin, vmax=vmax)
            self.cb  = plt.colorbar(self.img, label="Log$_{10}$ (Intensity)")
        else:
            self.img = self.ax.contourf(self.first_axis_angle, self.second_axis_angle, self.data, 50, vmin=vmin, vmax=vmax)
            self.cb  = plt.colorbar(self.img, label="Intensity")
        self.ax.set_xlabel(self.first_axis_name)
        self.ax.set_ylabel(self.second_axis_name)
        self.out_fig = os.path.splitext(self.filename_template)[0]+".png"
        
        plt.show()
        if save:
            self.fig.savefig(self.out_fig)
            
    def plot_1D(self, log=False, save=False):
        self.fig,self.ax = plt.subplots(1,1)
        self.fig.subplots_adjust(top=0.9, bottom=0.2, left=0.2,right=0.95,hspace=0.2,wspace=0.2)
        self.num_curves = self.second_axis_angle.size
        for i in range(self.num_curves):
            self.ax.plot(self.first_axis_angle, self.data[i], lw=2, label="{} = {:.2f}".format(self.second_axis_name, self.second_axis_angle[i]))
        if log:
            self.ax.set_yscale("log")
        self.ax.set_xlabel(self.first_axis_name)
        self.ax.set_ylabel("Intensity (cps)")
        self.ax.legend()
        self.out_fig = os.path.splitext(self.filename_template)[0]+".png"
        plt.show()
        if save:
            self.fig.savefig(self.out_fig)
            
def plot_one_curve(ras, log=False, save=False):
    fig,ax = plt.subplots(1,1)
    fig.subplots_adjust(top=0.9, bottom=0.2, left=0.2,right=0.95,hspace=0.2,wspace=0.2)
    ax.plot(ras.scan_angles, ras.scan_data, lw=2)
    if log:
        ax.set_yscale("log")
    ax.set_xlabel(ras.scan_axis)
    ax.set_ylabel("Intensity (cps)")
    out_fig = os.path.splitext(ras.filename)[0]+".png"
    plt.show()
    if save:
        fig.savefig(out_fig)
        
if __name__ == "__main__":
    
    if len(sys.argv)<2:
        print("Usage: python SmartLab.py  path_to_the_RAS_file(s)")
    else:
        allfiles= sys.argv[1:]
        first_file = allfiles[0]
        ras = RasFile(first_file)
        if "FILE_DATA_TYPE" in ras.header.keys():
            tmp = ras.header["FILE_DATA_TYPE"]
            if tmp == "RAS_3DE_RSM":
                mode = "RSM"
            elif tmp == "RAS_3DE_POLEFIG":
                mode = "PF"
        else:
            mode = "1D"
            
        nbr_files = len(allfiles)
        
        if mode.upper() == "PF":
            for fn in allfiles:
                pf=PoleFigure(fn)
                pf.plot()
        if mode.upper() == "RSM":
            for fn in allfiles:
                rsm = RSM(fn)
                rsm.plot_2D()
        if mode.upper() == "1D" and nbr_files>1:
            rsm = RSM(allfiles)
            rsm.plot_1D()
            rsm.plot_2D()
        if mode.upper() == "1D" and nbr_files==1:
            plot_one_curve(ras)