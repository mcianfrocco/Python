#!/usr/bin/env python

#This program will take an input text file of traces from the tracking program
#'Fiesta' (https://www.bcube-dresden.de/fiesta/wiki/FIESTA) and then curate the traces. 
#Curating involves discarding traces that are too short, move in the wrong direction, etc.   
#
#To run
#./curate_particles_calc_MSD.py
#
#And the help options will be display on the command line.
#
#Dependencies: python, numpy
#

import math
import glob
import linecache
import optparse
from sys import *
import os,sys,re
from optparse import OptionParser
import numpy as np

#=========================
def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog -i <input.txt> --brightLim=[brightnessLimit] --minRun=[minimum] --numFrames=[frames] --stdCutoff=[cutoff] --pixMovt=[pixels] -n")
        parser.add_option("-i",dest="input",type="string",metavar="FILE",
                help="Master particle list, output from Particle Tracker in ImageJ")
        parser.add_option("--brightLim",dest="brightLim",type="float", metavar="FLOAT",
                help="Brightness limit for first moment (ImageJ) or Intensity limit (Fiesta)")
	parser.add_option("--minRun",dest="minRun",type="float", metavar="FLOAT",
                help="Minimum run length (pixels)")
	parser.add_option("--numFrames",dest="numFrames",type="float", metavar="FLOAT",
                help="Number of frames to use")
	parser.add_option("--stdCutoff",dest="std",type="float", metavar="FLOAT",
                help="Standard deviation cutoff to remove traces X*stdDev away from average")
	parser.add_option("--pixMovt",dest="pixMovt",type="float", metavar="FLOAT",
                help="Minimum distance that a paricle must move during 'diffusion' (subpixel values OK)")
	parser.add_option("-n", action="store_true",dest="Flag",default=False,
                help="Flag to NOT consider standard deviation cutoff")
        parser.add_option("-f", action="store_true",dest="Fiesta",default=False,
                help="Flag if data file is from Fiesta")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 2:
                parser.print_help()
                sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

#==============
def checkConflicts(params):

	listRuns = glob.glob("%s_*"%(params['input'][:-4]))

	if len(listRuns) != 0:
		print 'Runs have already been parsed, exiting'
		sys.exit()
#===============
def getMaxTrajectory(params):

        numLine = len(open(params['input'],'r').readlines())
        lastLine = linecache.getline(params['input'],numLine)
        lastLineSplit = lastLine.split()
        maxTrajectory = float(lastLineSplit[1])
        if params['debug'] is True:
                print 'MaxTrajectory = %01d' %(maxTrajectory)
	
	return maxTrajectory

#================
def excludeBrightParticles(params):

	#To create a dictionary of bad traces due to brightness

	f1 = open(params['input'],'r')

	d = dict()
	traj = 1
	maxTrajectory = getMaxTrajectory(params)
	
	for line in f1: 

		l = line.split()
		trajectory = float(l[1])
		intensity = float(l[6])

		if intensity < params['brightLim']:
			if trajectory in d:
				continue
			else:
				d[trajectory]=intensity
	if params['debug'] is True:
		print d
		print len(d)	
	return d

#================
def excludeBrightParticles_Fiesta(params):

        #To create a dictionary of bad traces due to brightness

        f1 = open(params['input'],'r')

        d = dict()
	
        for line in f1:

                l = line.split()
		
		if len(line) == 1:
			if params['debug'] is True:
				print 'blank line'
			continue
	
                if l[0] == 'Molecule': 
			trajectory = float(l[1])
			if params['debug'] is True:
				print 'Trajectory %s == %s' %(l[1],line)
			continue

		if l[0] == 'frame': 
			if params['debug'] is True:
				print 'frame'
			continue
                intensity = float(l[6])

                if intensity < params['brightLim']:
                        if trajectory in d:
                                continue
                        else:
                                d[trajectory]=intensity
        if params['debug'] is True:
                print d
                print len(d)
        return d,trajectory

#================
def parseMasterParticleList(params,goodParticles):

	maxTrajectory = getMaxTrajectory(params)
	i = 1

	while i <= maxTrajectory:
		if i in goodParticles:
			f1 = open(params['input'],'r')
			o1 = open('%s_%01d.txt'%(params['input'][:-4],i),'w')
			
			counter = 1

			for line in f1:

				l = line.split()

				traj = l[1]
				intensity = float(l[6])
	
				if counter <= params['numFrames']:

					if float(traj) == i: 
						o1.write(line)
						counter = 1 + counter
			o1.close()
			f1.close()
		else:
			print "Particle %01d excluded, too bright" %(i)
	
		i = i + 1

#================
def parseMasterParticleList_Fiesta(params,goodParticles,maxTrajectory):

        i = 1

        while i <= maxTrajectory:
                if i in goodParticles:
                        f1 = open(params['input'],'r')
                        o1 = open('%s_%01d.txt'%(params['input'][:-4],i),'w')

                        counter = 1

                        for line in f1:

                                l = line.split()

                		if len(line) == 1:
                        		continue

                		if l[0] == 'Molecule':
                        		traj = float(l[1])
                        		continue

                		if l[0] == 'frame':
	                        	continue
	
                                if counter <= params['numFrames']:

                                        if traj == i:
                                                o1.write(line)
                                                counter = 1 + counter
                        o1.close()
                        f1.close()
                else:
                        print "Particle %01d excluded, too bright" %(i)

                i = i + 1

#================
def removeRunsTooShort(params):

	if params['minRun'] < params['numFrames']:
		limit = params['numFrames']

	listRuns = glob.glob("%s_*"%(params['input'][:-4]))
	
	for run in listRuns:

		length = len(open(run,'r').readlines())

		if length < limit:
			os.remove(run)
			print '%s removed, too few frames: %01d' %(run,length)	

#================
def mergeData(params,extension):

	listRuns = glob.glob("%s_*%s"%(params['input'][:-4],extension))

	count = 1

	for run in listRuns:
		if count == 1:
			merged = np.loadtxt(run)		
			count = 1 + count
		else:
			array = np.loadtxt(run)  
			merged = np.vstack((merged,array))
			count = 1 + count 

	merged = merged.T

	np.savetxt('%s_merge_%s'%(params['input'][:-4],extension),merged,fmt='%f')

#=================
def calcDistanceMSD(params):
	
	#pixel size (nm)
	pix = 159

	listRuns = glob.glob("%s_*"%(params['input'][:-4]))

	for run in listRuns:

		numLines = len(open(run,'r').readlines())
		i = 2
		total=0
		total2=0
	
		kymoy = open('%s_kymoY.txt'%(run[:-4]),'w')
		kymox = open('%s_kymoX.txt'%(run[:-4]),'w')
		outdist = open('%s_distY.txt'%(run[:-4]),'w')
		outcum = open('%s_cumulDistY.txt'%(run[:-4]),'w')
		MSD = open('%s_MSDY.txt'%(run[:-4]),'w')
		outdist2D = open('%s_dist2D.txt'%(run[:-4]),'w')
                outcum2D = open('%s_cumulDist2D.txt'%(run[:-4]),'w')
                MSD2D = open('%s_MSD2D.txt'%(run[:-4]),'w')

		while i <= numLines:

			l1 = linecache.getline(run,i-1).split()
			l2 = linecache.getline(run,i).split()

			if params['Fiesta'] is False:
	
				X1 = float(l1[3])
				X2 = float(l2[3])
				Y1 = float(l1[4])
				Y2 = float(l2[4])

			if params['Fiesta'] is True:
				X1 = float(l1[2])
                                X2 = float(l2[2])
                                Y1 = float(l1[3])
                                Y2 = float(l2[3])
				
			DIST2 = math.sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1))
			total2 = total2 + DIST2
			meansq2 = total2*total2

			DIST = math.fabs(Y2-Y1)*pix
			total = total + DIST
			meansq = total*total

			kymoy.write('%s\n'%(str(Y1/pix)))
			kymox.write('%s\n'%(str(X1/pix)))
			outdist.write('%s\n'%(str(DIST)))
			outcum.write('%s\n'%(str(total)))
			MSD.write('%s\n'%(str(meansq)))
			outdist2D.write('%s\n'%(str(DIST2)))
                        outcum2D.write('%s\n'%(str(total2)))
                        MSD2D.write('%s\n'%(str(meansq2)))
			
			i = i + 1

#=======================
def removeRunsNotOnMT(params,std_dev_cutoff):

	listRuns = glob.glob("%s_*"%(params['input'][:-4]))
	
	#Calculate average Y value
	counter=1
	average=0
	for run in listRuns:

		f1 = open(run,'r')

		for line in f1:

			l = line.split()
			Y = float(l[4])
			average = Y + average
			counter = 1 + counter	

	averageY = average/counter
	if params['debug'] is True:
		print 'Average Y value = %s' %(str(averageY))


	#Calculate standard deviation of Y
	
	for run in listRuns:
		
		f1 = open(run,'r')

		standardsum=0

		for line in f1:
			l = line.split()
			Y = float(l[4])

			std1 = (averageY-Y)*(averageY-Y)
			standardsum = std1 + standardsum
	stddev = math.sqrt(standardsum/counter)
	if params['debug'] is True:
		print 'Standard deviation of Y = %s' %(str(stddev))

	#remove files that are outside of range

	for run in listRuns:

		testline = linecache.getline(run,10)

                l = testline.split()
                Y = float(l[4])

		if Y > (averageY+(std_dev_cutoff*stddev)) or Y < (averageY-(std_dev_cutoff*stddev)):
			os.remove(run)
			print '%s removed. Y=%s, cutoff=%s-%s' %(run,str(Y),str(averageY+(std_dev_cutoff*stddev)),str(averageY-(std_dev_cutoff*stddev)))

#========================
def removeStaticRuns(params):

	listRuns = glob.glob("%s_*"%(params['input'][:-4]))

	for run in listRuns:
		#get first position

		if params['Fiesta'] is False:
			firstY = float(linecache.getline(run,1).split()[4])
			firstX = float(linecache.getline(run,1).split()[3])

		if params['Fiesta'] is True:
			firstY = float(linecache.getline(run,1).split()[3])
                        firstX = float(linecache.getline(run,1).split()[2])

		maxMovt = 0
                f1 = open(run,'r')
		for line in f1:
                        l = line.split()

			if params['Fiesta'] is False:
				Y = float(l[4])
				X = float(l[3])

			if params['Fiesta'] is True:
                                Y = float(l[3])
                                X = float(l[2])

			testMovt = math.sqrt((firstY-Y)*(firstY-Y)+(firstX-X)*(firstX-X))

			if testMovt > maxMovt:
				maxMovt = testMovt

		if maxMovt < params['pixMovt']:
			os.remove(run)
			print '%s removed. Total movt=%s; Movt cutoff=%s' %(run,str(maxMovt),str(params['pixMovt']))

#==============================
if __name__ == "__main__":

	params=setupParserOptions()
	checkConflicts(params)
	
	if params['Fiesta'] is False:
		goodParticles = excludeBrightParticles(params)	
		parseMasterParticleList(params,goodParticles)

	if params['Fiesta'] is True:
		goodParticles,maxTrajectory = excludeBrightParticles_Fiesta(params)
		parseMasterParticleList_Fiesta(params,goodParticles,maxTrajectory)

        removeRunsTooShort(params)
        if params['Flag'] is False:
        	removeRunsNotOnMT(params,params['std'])
        removeStaticRuns(params)
        calcDistanceMSD(params)
        mergeData(params,'MSDY.txt')
        mergeData(params,'MSD2D.txt')
        mergeData(params,'kymoY.txt')
        mergeData(params,'kymoX.txt')
