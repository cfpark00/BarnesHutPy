import numpy as np
import matplotlib.pyplot as plt
import glob
from vpython import *
import os
import glob




scene = canvas(width=900,height=800,background=vector(1,1,1))
full=0
size=1000.
objt=[]
r=10
dm=False
seeboxes=True
bound=None
trail=False
colcode={int(size/2):vector(101/255,67/255,33/255),int(size/4):vector(124/255,10/255,2/255),int(size/8):vector(220/255,20/255,5/255),
int(size/16):vector(1,36/255,0),int(size/32):vector(0,128/255,1),int(size/64):vector(115/255,194/255,251/255),int(size/128):vector(1,1,1),
int(size/256):vector(1,1,1),int(size/512):vector(1,1,1),int(size/1024):vector(1,1,1),int(size/2048):vector(1,1,1),int(size/4096):vector(1,1,1),
int(size/8192):vector(1,1,1),int(size/16384):vector(1,1,1),int(size/32768):vector(1,1,1)}
"""
colcode={int(size/2):vector(0,0,0),int(size/4):vector(0,0,0.5),int(size/8):vector(0,0,1),
int(size/16):vector(0,1,0.5),int(size/32):vector(0,1,1),int(size/64):vector(0,1,0),int(size/128):vector(1,1,0),
int(size/256):vector(1,0,0),int(size/512):vector(1,1,1),int(size/1024):vector(1,1,1),int(size/2048):vector(1,1,1),int(size/4096):vector(1,1,1),
int(size/8192):vector(1,1,1),int(size/16384):vector(1,1,1)}
"""
objtemp={}
obj={}
boxes={}
delkey=[]
while True:
    inputtt=input("FILEPATH: ")
    for key in obj:
        delkey.append(key)
    for key in delkey:
        obj[key].visible=False
        del obj[key]
        if seeboxes:
            boxes[key].visible=False
            del boxes[key]
        

    delkey=[]
    fl=glob.glob(inputtt)
    f=open(fl[0],"r")
    print("Reading file...")

    lincount=0
    for line in f:
        if lincount==0:
            lincount+=1
            continue
        a=line.strip().split(" ")
        if dm:
            if (int(a[8])==2):
                continue
        objtemp[int(a[0])]=[vector(float(a[1])-size/2,float(a[2])-size/2,float(a[3])-size/2),colcode[int(float(a[4]))],float(a[4]),float(a[5])-size/2,float(a[6])-size/2,float(a[7])-size/2]

    print("Initializing")

    
    

    def makebox(ox,oy,oz,s):
        box=curve(pos=[(ox,oy,oz), (ox+s,oy,oz), (ox+s,oy+s,oz),(ox,oy+s,oz),(ox,oy,oz)],color=vector(0,1,0), radius=3)
        box.append(pos=(ox,oy,oz+s))
        box.append(pos=(ox+s,oy,oz+s))
        box.append(pos=(ox+s,oy+s,oz+s))
        box.append(pos=(ox,oy+s,oz+s))
        box.append(pos=(ox,oy,oz+s))
        box.append(pos=(ox+s,oy,oz+s))
        box.append(pos=(ox+s,oy,oz))
        box.append(pos=(ox+s,oy+s,oz))
        box.append(pos=(ox+s,oy+s,oz+s))
        box.append(pos=(ox,oy+s,oz+s))
        box.append(pos=(ox,oy+s,oz))
        return box


    def drawuniversebound(ox,oy,oz,s):
        global bound
        bound=curve(pos=[(ox,oy,oz), (ox+s,oy,oz), (ox+s,oy+s,oz),(ox,oy+s,oz),(ox,oy,oz)],color=vector(1,1,1), radius=5)
        bound.append(pos=(ox,oy,oz+s))
        bound.append(pos=(ox+s,oy,oz+s))
        bound.append(pos=(ox+s,oy+s,oz+s))
        bound.append(pos=(ox,oy+s,oz+s))
        bound.append(pos=(ox,oy,oz+s))
        bound.append(pos=(ox+s,oy,oz+s))
        bound.append(pos=(ox+s,oy,oz))
        bound.append(pos=(ox+s,oy+s,oz))
        bound.append(pos=(ox+s,oy+s,oz+s))
        bound.append(pos=(ox,oy+s,oz+s))
        bound.append(pos=(ox,oy+s,oz))

    #drawuniversebound(-size/2,-size/2,-size/2,size)

    for key in objtemp:
        obj[key]=sphere(pos=objtemp[key][0],radius=r,color=vector(1,0,0),make_trail=trail,retain=30)
        if seeboxes:
            boxes[key]=makebox(objtemp[key][3],objtemp[key][4],objtemp[key][5],objtemp[key][2])
    sleep(0.0000001)