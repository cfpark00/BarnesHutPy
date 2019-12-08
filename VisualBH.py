import numpy as np
import time as time

deg=0
combideg=11
savedirectory="dataBH"
class node:
	def __init__(self,idnum,m,cm,leaf,size,origin,v=np.full(3,0.),rpast=None):
		self.id=idnum
		self.leaf=leaf
		self.children=None
		self.size=size
		self.m=m
		self.cm=cm
		self.origin=origin
		self.acc=np.full(3,0.)
		self.v=v
		self.delete=False
		self.rpast=rpast
	def __str__(self):
		a="{} with cm={} and mass {}".format("Leaf" if self.leaf else "Tree with {} children and depth {}".format(sum(x is not None for x in self.children),self.getdepth()),str(self.cm),str(self.m))
		if self.children==None:
			return a
		a+="\n ["
		for child in self.children:
			if child==None:
				a+=" ,"
				continue
			a+=str(child.m)
			a+=","
		a+="]"
		return a	
	def add_particle(self,idnum,m,r,v=np.full(3,0.),rpast=None):
		global deg,combideg
		if deg>combideg:
			self.mergeparticles(m,r,v)
			if idnum==-1:
				self.id=-1
			deg=0
			return 0
		octant=self.find_octant(r-self.origin)
		indice=self.octanttoindice(octant)
		if self.leaf:
			self.selfdeepen()
		if self.children[indice]==None:
			self.children[indice]=node(idnum,m,r,True,self.size/2,self.origin+octant*(self.size/2),v,rpast)
			self.cm=self.cm*self.m
			self.cm=self.cm+m*r
			self.m=self.m+m
			self.cm=self.cm/self.m
			deg=0
			return 0
		deg+=1
		self.children[indice].add_particle(idnum,m,r,v,rpast)
		self.cm=self.cm*self.m
		self.cm=self.cm+m*r
		self.m=self.m+m
		self.cm=self.cm/self.m
	def mergeparticles(self,m,r,v):
		self.v=(self.m*self.v+m*v)
		self.m=self.m+m
		self.v=self.v/self.m
	def writedata(self,n):
		f=open(savedirectory+"/t_"+str(n),"a")
		if self.leaf:
			self.writeone(f)
			return 0
		for child in self.children:
			if child==None:
				continue
			child.writeone(f)
		f.close()
	def writeone(self,f):
		if self.leaf:
			f.write(str(self.id)+" "+str(self.cm[0])+" "+str(self.cm[1])+" "+str(self.cm[2])+" "+str(self.size)+" "+str(self.origin[0])+" "+str(self.origin[1])+" "+str(self.origin[2]))
			f.write("\n")
		else:
			for child in self.children:
				if child==None:
					continue
				child.writeone(f)
	def getdepth(self):
		maxdepth=0
		if self.leaf:
			return 0
		for child in self.children:
			if child==None:
				continue
			a=child.getdepth()
			if a>maxdepth:
				maxdepth=a
		return maxdepth+1
	def selfdeepen(self):
		if not self.leaf:
			raise Exception("Not leaf, can\'t deepen")
		self.leaf=False
		self.children=[None for i in range(8)]
		octant=self.find_octant(self.cm-self.origin)
		indice=self.octanttoindice(octant)
		self.children[indice]=node(self.id,self.m,self.cm,True,self.size/2,self.origin+octant*(self.size/2),self.v,self.rpast)
		self.id=0
	def find_octant(self,r):#r should be given relative to position
		x=0
		y=0
		z=0
		if r[0]>self.size/2:
			x=1
		if r[1]>self.size/2:
			y=1
		if r[2]>self.size/2:
			z=1
		return np.array([x,y,z])
	def octanttoindice(self,v):

		return 4*v[0]+2*v[1]+v[2]
	def distance(self,othernode):

		return np.sqrt(np.sum(np.square(self.cm-othernode.cm)))
	def criterion(self,othernode):

		return othernode.size/self.distance(othernode)
	def calcvone(self):
		global dt
		self.v=np.add(self.v,self.acc*dt)
	def calcaccfrom(self,othernode):
		if not self.leaf:
			raise Exception('Force calc for non-leaf')
		if othernode.leaf:
			if not othernode==Universe:
				raise Exception('Bad access for acc calc')
			print("One particle universe")
			return 0
		for child in othernode.children:
			if child==None:
				continue
			if self==child:
				continue
			if child.leaf:
				physicsacc(self,child)
				continue
			if self.criterion(child)<threshold:
				physicsacc(self,child)
			else:
				self.calcaccfrom(child)


def physicsacc(ori,b):#all force-physics is in here
	global remove_singularity, presslen, eps, r_eq,a,n,singlemass, externalpressure,meshsize,gf
	if remove_singularity:
		if thmodel=="simple":
			ori.acc=np.add(ori.acc,(b.cm-ori.cm)*G*b.m/(ori.distance(b)**2+presslen**2)/ori.distance(b))
			return 0
		if thmodel=="vdw":
			ori.acc=np.add(ori.acc,((b.cm-ori.cm)/ori.distance(b))*(12*eps)*r_eq*((r_eq/ori.distance(b))**(-13)-2*(r_eq/ori.distance(b))**(-7)))
			ori.acc=np.add(ori.acc,(b.cm-ori.cm)*G*b.m/(ori.distance(b)**2+presslen**2)/ori.distance(b))
			return 0
	ori.acc=np.add(ori.acc,(b.cm-ori.cm)*G*b.m/ori.distance(b)**3)

def physicsadv(a,fromfnc):
	global dt,method, t
	if method=="euler" and fromfnc==1:
		a.calcvone()
		return 0
	if fromfnc==1:
		return 0
	if method=="euler":
		a.cm=a.cm+a.v*dt
		return 0
	if t==0:
		a.rpast=a.cm
		a.cm=a.cm+a.v*dt+0.5*a.acc*(dt**2)
	else:
		temp=a.cm
		a.cm=2*a.cm-a.rpast+a.acc*(dt**2)
		a.rpast=temp



def calcacc(self):
	global threshold,Universe,G,meshsize,a
	if self.leaf:
		if self.id==-1:
			return 0
		self.acc=np.full(3,0.)
		self.calcaccfrom(Universe)
		if externalpressure=="normal":
			#print(self,self.id,self.cm,self.rpast,meshsize,a)
			self.acc=np.add(self.acc,gf[int(self.cm[0]*meshsize/a),int(self.cm[1]*meshsize/a),int(self.cm[2]*meshsize/a)])
		physicsadv(self,1)
		return 0
	for child in self.children:
		if child==None:
			continue
		calcacc(child)

def advance(self):
	global Universe,singlemass
	if self.leaf:
		if self.id==-1:
			return 0
		physicsadv(self,0)
		return 0
	for child in self.children:
		if child==None:
			continue
		advance(child)

def feed(datanode):
	global nodes,method
	if datanode.leaf:
		nodes.append([datanode.id,datanode.m,datanode.cm,datanode.v,None if method=="euler" else datanode.rpast])
		return 0
	for child in datanode.children:
		if child==None:
			continue
		feed(child)

def formatr(rin):
	global a
	if bc=="torus":
		return rin%a
	if bc=="box":
		if rin[0]<0:
			rin[0]*=-1
		elif rin[0]>a:
			rin[0]=2*a-rin[0]
		if rin[1]<0:
			rin[1]*=-1
		elif rin[1]>a:
			rin[1]=2*a-rin[1]
		if rin[2]<0:
			rin[2]*=-1
		elif rin[2]>a:
			rin[2]=2*a-rin[2]
		return rin


def reconstruct():
	global a,Universe,nodes,method
	Universe=node(nodes[0][0],nodes[0][1],formatr(nodes[0][2]),True,a,np.full(3,0.),nodes[0][3],None if method=="euler" or nodes[0][0]==-1 else formatr(nodes[0][4]))
	for i in range(1,len(nodes)):
		Universe.add_particle(nodes[i][0],nodes[i][1],formatr(nodes[i][2]),nodes[i][3],None if method=="euler" or nodes[i][0]==-1 else formatr(nodes[i][4]))
	nodes=[]

nodes=[]

##################### START ############################

#############
#Parameters
#############
G=2000
singlemass=100
a=1000  #oneside length
n=100 #number of particles
threshold=1 #the Barnes-Hut threshold
dt=0.5
totaliternum=1000



#############
#options
#############

###
# Numerical integration
###
method="verlet"#Choose between "euler","verlet","RK4"

###
# Interaction physics
###

coagulation=True
#degree at top of code

remove_singularity=True #Effect of thermal pressure
#Themal model
thmodel="simple"#choose between "ideal", "vdw", "simple"
#force parameter eps and equilibrium distance r_eq for "vdw"
eps=1000.
r_eq=100.
#Pressure length for "simple"
presslen=100. #effective length of pressure


###
# Boundary condition
###
bc="box" #choose between "torus", "box", "qinfinite"
externalpressure="normal"
usepotdata=False
writepotdata=True
potdatafilename="potdata16"

### Hubble expansion
expansion=True
H=100


###
# Systems
###
central=False
centralmass=10000
vtangential="norma"#"normal"#"normal" for disc rotation
vrand=False
vmag=20#default not real magnitude
# Disc initial condition
disc=False


####
# Calculate balancing potential
####
meshsize=16
if usepotdata and externalpressure=="normal":
	gf=np.zeros(shape=(meshsize,meshsize,meshsize), dtype=(float,3))
	potdatfile=open(potdatafilename,"r")
	print("Reading potential "+potdatafilename)
	for line in potdatfile:
		l=line.replace(":"," ").split(" ")
		gf[int(l[0])][int(l[1])][int(l[2])][0]=float(l[3])
		gf[int(l[0])][int(l[1])][int(l[2])][1]=float(l[4])
		gf[int(l[0])][int(l[1])][int(l[2])][2]=float(l[5])
	potdatfile.close()
elif externalpressure=="normal":
	fac=G*n*singlemass/meshsize**3
	rat=a/meshsize
	gf=np.zeros(shape=(meshsize,meshsize,meshsize), dtype=(float,3))
	#Intense calculation
	count=0.
	for x in range(meshsize):
		for y in range(meshsize):
			print(str(count*100/meshsize**2)+" % of potential calculation done")
			for z in range(meshsize):
				for xx in range(meshsize):
					for yy in range(meshsize):
						for zz in range(meshsize):
							rmag3=(rat**(3))*((x-xx)**2+(y-yy)**2+(z-zz)**2)**(3/2)
							if rmag3==0:
								continue
							gf[x][y][z][0]+=(x-xx)*(rat)*fac/rmag3
							gf[x][y][z][1]+=(y-yy)*(rat)*fac/rmag3
							gf[x][y][z][2]+=(z-zz)*(rat)*fac/rmag3
			count+=1
if writepotdata and not usepotdata and externalpressure=="normal":
	fpot=open(savedirectory+"/potdata"+str(meshsize)+",GM"+str(fac),"w")
	for x in range(meshsize):
		for y in range(meshsize):
			for z in range(meshsize):
				fpot.write(str(x)+" "+str(y)+" "+str(z)+":"+str(gf[x][y][z][0])+" "+str(gf[x][y][z][1])+" "+str(gf[x][y][z][2])+"\n")
	fpot.close()

print("potdata read")
#############
#Initialize
np.random.seed()

if central:
	Universe=node(-1,centralmass,np.full(3,a/2.),True,a,np.full(3,0.),np.full(3,0.))
else:
	r=a*np.array([np.random.random(),np.random.random(),np.random.random()])
	v=vmag*np.array([np.random.random(),np.random.random(),np.random.random()])
	if not vrand:
		v=np.full(3,0.)#v is zero
	Universe=node(1,singlemass,r,True,a,np.full(3,0.),v)
print(Universe)
for i in range(n-1):
	r=a*np.array([np.random.random(),np.random.random(),np.random.random()])
	rcent=r-np.array([a/2.,a/2.,a/2.])
	if disc:
		rcent[2]=rcent[2]/10.
		r[2]=rcent[2]+a/2.
	rcentmag=np.sqrt(np.sum(np.square(rcent)))
	v=vmag*np.array([np.random.random(),np.random.random(),np.random.random()])
	if not vrand:
		v=np.full(3,0.)#v is zero
	if vtangential=="normal":
		vmag=np.sqrt(G*centralmass/rcentmag)
		v=vmag*np.array([-rcent[1]/rcentmag,rcent[0]/rcentmag,np.random.random()/10.])
	Universe.add_particle(i+2,singlemass,r,v)
	print(Universe)

t=0
tlog=[0.,0.,0.,0.]
t00=time.time()
for i in range(totaliternum):
	t0=time.time()
	calcacc(Universe)
	t1=time.time()
	tlog[0]+=t1-t0
	print("Time for force calc: "+str(t1-t0))
	advance(Universe)
	t2=time.time()
	tlog[1]+=t2-t1
	print("Time for advance: "+str(t2-t1))
	feed(Universe)
	reconstruct()
	t3=time.time()
	tlog[2]+=t3-t2
	print("Time for reconstruct: "+str(t3-t2))
	f=open(savedirectory+"/t_"+str(int(10*t)),"w")
	f.write("StepTime "+str(t3-t0)+"\n")
	f.close()
	Universe.writedata(int(10*t))
	t4=time.time()
	tlog[3]+=t4-t3
	print("Time for data writing: "+str(t4-t3))
	print(Universe)
	t=t+dt
	print("Progress "+str(i*100./totaliternum)+"%")
	if not i==0:
		print("Estimated time left: "+str((totaliternum-i)*(time.time()-t00)/(1.*i))+" seconds.")
	print("\n")
	
print("end")

print(thmodel+" thermal model computation completed with dt= "+ str(dt)+" for "+ str(totaliternum*dt)+" total time.")
print(method+" method used for numerical integration.")
print(" ")
print("Time for force calc: "+str(tlog[0])+". This is "+str(tlog[0]*100./np.sum(tlog))+" %")
print("Time for advance: "+str(tlog[1])+". This is "+str(tlog[1]*100./np.sum(tlog))+" %")
print("Time for reconstruct: "+str(tlog[2])+". This is "+str(tlog[2]*100./np.sum(tlog))+" %")
print("Time for data writing: "+str(tlog[3])+". This is "+str(tlog[3]*100./np.sum(tlog))+" %")

logfile=open(savedirectory+"/RUNLOG","w")
logfile.write("Data run ended at "+time.ctime(int(time.time()))+"\n")
logfile.write("Time for force calc: "+str(tlog[0])+". This is "+str(tlog[0]*100./np.sum(tlog))+" %\n")
logfile.write("Time for advance: "+str(tlog[1])+". This is "+str(tlog[1]*100./np.sum(tlog))+" %\n")
logfile.write("Time for reconstruct: "+str(tlog[2])+". This is "+str(tlog[2]*100./np.sum(tlog))+" %\n")
logfile.write("Time for data writing: "+str(tlog[3])+". This is "+str(tlog[3]*100./np.sum(tlog))+" %\n")
logfile.write("G= "+str(G)+"\n")
logfile.write("n= "+str(n)+"\n")
logfile.write("singlemass= "+str(singlemass)+"\n")
logfile.write("a(total size)= "+str(a)+"\n")
logfile.write("Barnes-Hut Threshold= "+str(threshold)+"\n")
logfile.write("dt= "+str(dt)+"\n")
logfile.write("totaliternum= "+str(totaliternum)+"\n")
logfile.write("combideg= "+str(combideg)+"\n")
logfile.write("method: "+str(method)+"\n")
logfile.write("coagulation= "+str(coagulation)+"\n")
logfile.write("remove_singularity= "+str(remove_singularity)+"\n")
logfile.write("thermal model= "+str(thmodel)+"\n")
logfile.write("eps= "+str(eps)+"\n")
logfile.write("r_eq= "+str(r_eq)+"\n")
logfile.write("presslen= "+str(presslen)+"\n")
logfile.write("boundary condition= "+str(bc)+"\n")
logfile.write("externalpressure= "+str(externalpressure)+"\n")
logfile.write("central= "+str(central)+"\n")
logfile.write("centralmass= "+str(centralmass)+"\n")
logfile.write("vtangential= "+str(vtangential)+"\n")
logfile.write("vrand= "+str(vrand)+"\n")
logfile.write("vmag= "+str(vmag)+"\n")
logfile.write("disc= "+str(disc)+"\n")

logfile.close()


