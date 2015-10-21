#!usr/bin/env python

import os, pdb, sys, time, getopt, decimal, math
import argparse, logging, logging.handlers, random
from astropy.units import cds
import numpy as np
import manga_utils.dateObs2HA as ha

try:
	from platedb.APODatabaseConnection import db
	import platedb.ModelClasses as platedb
	import mangadb.ModelClasses as mangadb
	import sqlalchemy
except:
	print '\n Make sure to set up platedb and mangadb before running'
	print '\n\n'
	raise


### Initialize the main logger
def initLogger():
	'''Initialize the main DOS logger'''
	
	logfile='loadMockData'
	
	log = logging.getLogger(logfile)
	fmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s','%Y-%m-%d %H:%M:%S')
	log.setLevel(logging.DEBUG)
    	
	# main log
	handler = logging.handlers.TimedRotatingFileHandler(logfile+'.log', 'midnight', 1, 5)
	handler.setFormatter(fmt)
	handler.setLevel(logging.DEBUG)
	log.addHandler(handler)

	# to console
	confmt = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s','%Y-%m-%d %H:%M:%S')	
	console = logging.StreamHandler()
	console.setLevel(logging.INFO)
	console.setFormatter(confmt)
	log.addHandler(console)
	
	return log

### Parse command line arguments
def parseArgs():
	'''Parse the command line arguments'''

	parser = argparse.ArgumentParser(prog='loadMockData', usage='%(prog)s [options]', 
			formatter_class=argparse.RawDescriptionHelpFormatter, 
		    description='loads mock data into platedb and mangadb')

	parser.add_argument('-s', '--set', type=int, dest='set', help='number of exposures in one set, also splits up the MJDs', default=9)		   
	parser.add_argument('-n', '--number', type=int, dest='number', help='number of entries to load into the database', default=2)
	parser.add_argument('-p', '--plate', action='store_false', dest='plate', help='turn off loading into platedb, default is to load', default=True)
	parser.add_argument('-m', '--manga', action='store_false', dest='manga', help='turn off loading into mangadb, default is to load', default=True)
	parser.add_argument('-b', '--database', action='store_true', dest='dboff', help='turn on loading into datadb, default is to load', default=False)
	parser.add_argument('-d', '--delete', action='store_true', dest='delete', help='delete all entries in dbs', default=False)
	parser.add_argument('-l', '--list', action='store_true', dest='list', help='list all entries instead of loading into dbs', default=False)

	opts = parser.parse_args()
					
	return opts

### get dither offsets
def getDOffs(dpos):
	'''get dither offset'''
	
	if dpos == 'N':
		return (-0.417, +0.721)
	if dpos == 'S':
		return (-0.417, -0.721)
	if dpos == 'E':
		return (+0.833, +0.000)	

	 		
### random object class
class Object:
	'''random object class'''
	
	def __init__(self, number=None, mjd=56741):
		self.number = number
		self.mjd = mjd
		self.plate = 7443
		self.expNum = 00177375+self.number
		self.startTime = (self.mjd*cds.MJD).to(cds.s).value
		self.expTime = 900.1
		self.camera = random.choice(['b1','b2','r1','r2'])
		self.filename = 'mgscisky-{0}-{1}-{2}.fits'.format(self.plate,self.camera,self.exposure)
		self.survey = 'MaNGA'
		self.flavor = 'Science'
		self.sn2 = random.uniform(1,9)
		self.dpos = (number % 3 == 0) and 'N' or (number % 3 == 1) and 'S' or (number % 3 == 2) and 'E'
		self.dra = getDOffs(self.dpos)[0]
		self.ddec = getDOffs(self.dpos)[1]
		self.ha = random.uniform(0,360)
		self.seeing = random.uniform(0.8,2.6)
		self.transparency = random.uniform(0.5,1.0)
	
	def __repr__(self):
		return self.__str__()
		
	def __str__(self):
		return ('Number: {0} \n'.format(self.number) + 'MJD: {0} \n'.format(self.mjd)+
		'expNum: {0} \n'.format(self.expNum) + 'startTime {0} \n'.format(self.startTime)+
		'expTime: {0} \n'.format(self.expTime) + 'camera {0} \n'.format(self.camera)+
		'survey: {0} \n'.format(self.survey) + 'flavor {0} \n'.format(self.flavor)+
		'sn2: {0} \n'.format(self.sn2) + 'dpos: {0} \n'.format(self.dpos)+
		'dra: {0} \n'.format(self.dra) + 'ddec: {0} \n'.format(self.ddec)+
		'ha: {0} \n'.format(self.ha) + 'seeing: {0} \n'.format(self.seeing)+
		'transparency: {0} \n'.format(self.transparency));

	def getSN2(self, num):
		if num == 1: return random.gauss(3.6,2)
		if num == 2: return random.gauss(7.5,2)
	
### Load the Mock Data into Platedb
def loadPlateDB(session, obj,log):
	'''load the mock data'''
	
	# Plate DB loading
	log.info('Loading into platedb')
	
	# create plate here
	log.info('Creating plate {0}'.format(obj.number))
	plate=None			

	# create plugging here
	log.info('Creating plugging {0}'.format(obj.number))
	plugging=None
						
	# create plate pointing here
	log.info('Creating plate pointing {0}'.format(obj.number))
	platePointing=None
			
	# create observation
	log.info('Creating observation {0}'.format(obj.number))
	observation = platedb.Observation()
	session.add(observation)
	observation.plate_pointing = platePointing
	observation.plugging = plugging
	observation.mjd = obj.mjd
	observation.comment = ''

	# create exposure
	log.info('Creating exposure {0}'.format(obj.number))
	camera  = session.query(platedb.Camera).filter_by(label=obj.camera).one()
	survey  = session.query(platedb.Survey).filter_by(label='MaNGA').one()
	flavor  = session.query(platedb.ExposureFlavor).filter_by(label='Science').one()
	exposure = platedb.Exposure()
	session.add(exposure)
	exposure.exposure_no = obj.expNum
	exposure.observation = observation
	exposure.start_time = obj.startTime
	exposure.exposure_time = obj.expTime
	exposure.camera = camera
	exposure.survey = survey
	exposure.flavor = flavor
	exposure.comment = ' '
	exposure.exposure_status_pk = 1 #excellent
			
	# create camera frame
	log.info('Creating camera frame {0}'.format(obj.number))
	cframe = platedb.CameraFrame()
	session.add(cframe)
	cframe.exposure = exposure
	cframe.camera = camera
	cframe.survey = survey
	cframe.comment = ''
	cframe.sn2 = obj.sn2
			
	# flush the plate db session
	if observation and exposure and cframe:
		log.info('Flushing the platedb session.')
		session.flush()

###
### Load the Mock Data into Mangadb
def loadMangaDB(session, obj,log):
	'''load the mock data'''

	# get plate db exposure
	try:
		pdbexposure = session.query(platedb.Exposure).with_lockmode('update').filter_by(exposure_no=obj.expNum).one()
	except:
		pdbexposure=None
		
	# create exposure
	log.info('Creating mangadb exposure {0}'.format(obj.number))
	exposure = mangadb.Exposure()
	session.add(exposure)
	exposure.platedb_exposure_pk = None if pdbexposure==None else pdbexposure.pk
	exposure.dither_position = obj.dpos
	exposure.dither_ra = obj.dra
	exposure.dither_dec = obj.ddec
	exposure.ha = obj.ha
	exposure.seeing = obj.seeing
	exposure.transparency = obj.transparency
	exposure.comment = ' '
	exposure.exposure_status_pk = 1 #excellent

	# create sn2_values
	cam = obj.camera
	pipename = session.query(mangadb.PipelineName).filter(mangadb.PipelineName.label == 'DOS').one()
	sn2val = mangadb.SN2Values()
	session.add(sn2val)
	sn2val.b1_sn2 = obj.getSN2(1)
	sn2val.r1_sn2 = obj.getSN2(2)
	sn2val.b2_sn2 = obj.getSN2(1)
	sn2val.r2_sn2 = obj.getSN2(2)
	sn2val.exposure = exposure
	sn2val.pipelinename = pipename		
			
	# flush the plate db session
	if exposure:
		log.info('Flushing the mangadb session.')
		session.flush()

### load mock data
def loadMockData(session, opts, log):
	'''Load mock data'''
	
	log.info('Loading {0} number of entries into platedb, mangadb'.format(opts.number))

	# set up MJDs
	mjdlist = [[i]*opts.set for i in np.arange(int(math.ceil(opts.number/float(opts.set))))+56741]
	mjd=[inner for outer in mjdlist for inner in outer] # nested comprehension loop

	# loop over all objects
	for entry in range(opts.number):
		# generate random object
		obj = Object(entry+1, mjd=mjd[entry])
				
		#pdb.set_trace()
		if opts.list:
			log.info('Object: \n {0}'.format(obj))
		else:	
			# load into platedb
			if opts.plate == True : loadPlateDB(session, obj, log)
			# load into mangadb
			if opts.manga == True : loadMangaDB(session, obj, log)

### delete all 
def deleteAll(session, opts, log):
	'''delete everything'''
	
	log.info('Deleting everything from platedb and mangadb')
	manga = session.query(mangadb.Exposure).delete()
	obs = session.query(platedb.Observation).delete()
	exp = session.query(platedb.Exposure).delete()
	cam = session.query(platedb.CameraFrame).delete()
		
### Main
def main(args):
	'''Main method to load mock data into platedb and mangadb'''

	# parse command line
	opts = parseArgs()

	# initialize logger to console
	log=initLogger()
	
	# start DB session and load mock data
	session=None
	try:
		session=db.Session()
		session.begin(subtransactions=True)
		if opts.delete == False:
			log.info('loading mock data')
			loadMockData(session, opts, log)
		else:	
			log.info('deleting all entries')
			deleteAll(session, opts, log)
			
		session.commit()
	except:
		if session != None:
			log.error('Error occurred while loading mock data.  Rolling back session changes.')
			session.rollback()
		raise
	finally:
		if session != None:
			session.close()
			db.engine.dispose()
						
###
if __name__ == '__main__':
	main(sys.argv[1:])

