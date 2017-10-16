#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
import sys
import os

from math import cos,sin,radians,sqrt,fabs,acos,pi,degrees,floor,ceil
import struct
summArray=[]

patToFIdir = "/home/egor/dose/FI/"
class GeoGrid:
	def __init__(self,countx,county,lonmin,latmin,dlon,dlat,lonmax,latmax):
		self.lonmin = lonmin
		self.lonmax = lonmax
		self.latmin = latmin
		self.latmax = latmax
		self.countx = countx
		self.county = county
		self.fmin   = 0.0
		self.fmax   = 0.0		
		self.dlon = dlon
		self.dlat = dlat
		self.data   = [0.0 for o in range(countx*county)]
	
	def findMin(self):
		if(len(self.data) == 0):
			return -1
		ret = self.data[0]
		for v in self.data:
			if(ret > v):
				ret = v
		return ret
		
	def findMax(self):
		if(len(self.data) == 0):
			return -1
		ret = self.data[0]
		for v in self.data:
			if(ret < v):
				ret = v
		return ret			

	def printASCIIGRDFile(self,filename):
		f = open(filename, 'wt')
		f.write("DSAA" + '\r\n')
		f.write(str(self.countx)+"\t"+str(self.county) + '\r\n')
		f.write('%0.12f' % (self.lonmin)+"\t"+'%0.12f' % (self.lonmax) + '\r\n')
		f.write('%0.12f' % (self.latmin)+"\t"+'%0.12f' % (self.latmax) + '\r\n')
		f.write('%0.12e' % self.findMin()+"\t"+'%0.12e' % self.findMax() + '\r\n')
		for j in range(0,self.county):
			for i in range(0,self.countx):
				f.write('%0.12e' % self.data[i*self.county+j]+'\t')
			f.write('\r\n')
		
		f.close()
		return

	def getValue(self,lon,lat):
		if(lon<self.lonmin or lon>self.lonmax):
			return 0.
		if(lat<self.latmin or lat>self.latmax):
			return 0.
		ci = int(floor((lon-self.lonmin)/self.dlon))
		cj = int(floor((lat-self.latmin)/self.dlat))
		ci1 = 0
		ci2 = 0
                cj1 = 0
		cj2 = 0
		
		if ci == self.countx-1:
			ci2 = ci
			ci1 = ci-1
		else:
			ci1 = ci
			ci2 = ci+1
	    
		if cj == self.county-1:
			cj2 = cj
			cj1 = cj-1
		else:
			cj1 = cj
			cj2 = cj+1
		
		f11 = self.data[ci1*self.county+cj1] #+ -
		f21 = self.data[ci2*self.county+cj1] #+ -

		f12 = self.data[ci1*self.county+cj2] #- +
		f22 = self.data[ci2*self.county+cj2] #- +
		
		
		x1 = self.lonmin+self.dlon*ci1
		x2 = self.lonmin+self.dlon*ci2
		y1 = self.latmin+self.dlat*cj1
		y2 = self.latmin+self.dlat*cj2
		
		fy1=(f21-f11)/(x2-x1)*lon+(f11-(f21-f11)/(x2-x1)*x1)

		fy2=(f22-f12)/(x2-x1)*lon+(f12-(f22-f12)/(x2-x1)*x1)

		v = (fy2-fy1)/(y2-y1)*lat+(fy1-(fy2-fy1)/(y2-y1)*y1)
		return  v
	
	def angleOf(self,dX,dY):
		dFi = 0.0
		dR = sqrt(fabs(dX*dX+dY*dY));
		if (fabs(dR) > 0):
			dFi  = acos(dX/dR);
		if (dX <= 0 and dY < 0):
			dFi  = 2.0*pi - dFi;
		if (dX >  0 and dY < 0):
			dFi  = 2.0*pi - dFi;
		return dFi;
		
def prepareArray(num, path):
		 
		global summArray
		
		myWorkPath = str(path+"/res/")
		d = patToFIdir
		fileName = str(d+"f"+num+".bin")
		Name = ("f"+num+".bin")
		inp = open(str(fileName),"rb")
		datastr = inp.read(64)
		offsetStart = 64
		#float(fil.countx),float(fil.county),fil.lonmin,fil.latmin,fil.dlon,fil.dlat,fil.lonmax,fil.latmax
		headerStr = list(struct.unpack('<%dd' % (len(datastr)/8), datastr))
		print headerStr
		grid95 = GeoGrid(int(headerStr[0]),int(headerStr[1]),headerStr[2],headerStr[3],headerStr[4],headerStr[5],headerStr[6],headerStr[7])
		gridMax = GeoGrid(int(headerStr[0]),int(headerStr[1]),headerStr[2],headerStr[3],headerStr[4],headerStr[5],headerStr[6],headerStr[7])
		countx = int(headerStr[0])
		county = int(headerStr[1])
		print grid95.lonmin, grid95.latmin
		count = countx*county
		offsetForNextGrid = countx*county*4
		gridsCount = (os.stat(fileName).st_size-offsetStart)/offsetForNextGrid-1
		posOf95Pbase = int(ceil(float(gridsCount)*1.00))-1 #номер сетки превращаем в позицию сетки т.е. -1
		print gridsCount,offsetForNextGrid,posOf95Pbase
		#отличие от реализации Асфандиярова - работаем с данными точности float, чтобы ограничивать объемы дискового пространства
		summArray = [[0.0]*gridsCount]*count
		for pp in range(count):
			arr = []
			for grid in range(gridsCount):
				offset = offsetStart + grid*offsetForNextGrid+4*pp
				inp.seek(offset)
				binDat = inp.read(4)
				arr.append(struct.unpack('<f', binDat)[0])
			arr.sort()
			summArray[pp]=summArray[pp]+arr
			#print arr
			grid95.data[pp] = arr[posOf95Pbase]
		#print grid95.data
		print summArray[pp]
		inp.close()
		#grid95.printASCIIGRDFile(myWorkPath+Name.split(".")[0]+".grd")
		print "ok!"

def main():
	global summArray
	mypath=os.path.dirname(os.path.realpath( __file__ ))
	os.chdir(mypath)
	summArray
	print "Vvedite cherez probel nomera funktsionalov"
	string = str(input())
	lst = []
	if string.find(" "):
		lst=string.split(' ')
	else:
		lst = int(string)
	for number in lst:
		prepareArray(number, mypath)
	return

if __name__ == "__main__":
	sys.exit(main())

