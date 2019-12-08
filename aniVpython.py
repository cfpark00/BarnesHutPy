import numpy as np
import matplotlib.pyplot as plt
import glob
from vpython import *
import os
from PIL import ImageGrab


import glob
def sortby(a):
    return int(a.split("_")[1])
fl=glob.glob("dataBH/t*")
fl.sort(key=sortby)

full=0
size=30000.
objt=[]
trail=False

scene = canvas(width=900,height=800)

print("Reading file...")
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
numfiles=len(fl)
count=0
skipping=0
endding=numfiles
itercalctime=[]
dm=False
for n in fl:
    if count<skipping:
        count+=1
        continue
    if count>endding:
        break
    f=open(n,"r")
    objtemp={}
    lincount=0
    for line in f:
        if lincount==0:
            itercalctime.append(float(line.strip().split(" ")[1]))
            lincount+=1
            continue
        a=line.strip().split(" ")
        #print(n)
        #print(line)
        if dm:
            if int(a[8])==0:
                continue
        objtemp[int(a[0])]=[vector(float(a[1])-size/2,float(a[2])-size/2,float(a[3])-size/2),colcode[int(float(a[4]))],float(a[4]),float(a[5])-size/2,float(a[6])-size/2,float(a[7])-size/2]
    objt.append(objtemp)
    full+=1
    count+=1
    print(str(count*100./numfiles)+" % done")

del objtemp

ff=open("calctimehist","w")
for el in itercalctime:
    ff.write(str(el)+" ")
ff.close()

print("Initializing")

r=70#size of particle visual
col=vector(1,0,0)

obj={}
seeboxes=False

def makebox(ox,oy,oz,s):
    box=curve(pos=[(ox,oy,oz), (ox+s,oy,oz), (ox+s,oy+s,oz),(ox,oy+s,oz),(ox,oy,oz)],color=vector(57/255,1,20/255), radius=10)
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

    """
    box.append(curve(pos=[(ox,oy,oz+s), (ox+s,oy,oz+s), (ox+s,oy+s,oz+s),(ox,oy+s,oz+s),(ox,oy,oz+s)], radius=30))
    box.append(curve(pos=[(ox,oy,oz), (ox,oy,oz+s)], radius=30))
    box.append(curve(pos=[(ox+s,oy,oz), (ox+s,oy,oz+s)], radius=30))
    box.append(curve(pos=[(ox,oy+s,oz), (ox,oy+s,oz+s)], radius=30))
    box.append(curve(pos=[(ox+s,oy+s,oz), (ox+s,oy+s,oz+s)], radius=30))
    """
    return box

boxes={}

initialids=[]
for key in objt[0]:
    obj[key]=sphere(pos=objt[0][key][0],radius=r,color=objt[0][key][1],make_trail=trail,retain=30)
    if seeboxes:
        boxes[key]=makebox(objt[0][key][3],objt[0][key][4],objt[0][key][5],objt[0][key][2])
    initialids.append(key)

bound=None


def drawuniversebound(ox,oy,oz,s):
    global bound
    bound=curve(pos=[(ox,oy,oz), (ox+s,oy,oz), (ox+s,oy+s,oz),(ox,oy+s,oz),(ox,oy,oz)],color=vector(1,1,1), radius=30)
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

drawuniversebound(-size/2,-size/2,-size/2,size)

def move():
    global iterations
    checked=[]
    for key in objt[iterations]:
        obj[key].pos=objt[iterations][key][0]
        obj[key].color=objt[iterations][key][1]
        if seeboxes:
            boxes[key]=makebox(objt[iterations][key][3],objt[iterations][key][4],objt[iterations][key][5],objt[iterations][key][2])
        checked.append(key)

    for oneid in initialids:
        if oneid in checked:
            continue
        obj[oneid].visible=False
        if seeboxes:
            boxes[oneid].visible=False
        initialids.remove(oneid)
        del obj[oneid]
        if seeboxes:
            del boxes[oneid]



iterations=0
save=False
while iterations<full:
    if iterations==0:
        pass
    move()
    if save:
        im = ImageGrab.grab((0,0,500,500))
        im.save("IMAGETEMP/imgtemp"+str(iterations)+".png")
    sleep(0.005)
    iterations+=1
    print(iterations)
    if False and mag(objt[iterations][1][0]-objt[iterations-1][1][0])>2000: #for debugging
        print(iterations,objt[iterations][1][0],objt[iterations-1][1][0],objt[iterations+1][1][0])
        raise Exception("What?")

