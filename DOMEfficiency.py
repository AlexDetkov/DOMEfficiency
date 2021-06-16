#!/usr/bin/python

import numpy as np
from icecube.dataio import I3File
from icecube.dataclasses import I3Position
from icecube.phys_services import Cylinder, I3Calculator
from icecube import simclasses, recclasses

# Variables
i3GeoLocation = '/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE.i3.gz'
detectorPadding = 0

trackBinSize = 50
stochasticLength = 200

class CherenkovHit():
	def __init__(self, charge, cherenkovDistance, trackDistance):
		self.charge = charge
		self.cherenkovDistance = cherenkovDistance
		self.trackDistance = trackDistance
	def intensity(self):
		return self.charge * np.exp(self.cherenkovDistance/stochasticLength) * self.cherenkovDistance
	def __repr__(self):
		return "CherenkovHit({}, {}, {})".format(self.charge, self.cherenkovDistance, self.trackDistance)

class CherenkovHitsBins(dict):
	def __init__(self, *args, **kw):
		super(CherenkovHitsBins, self).__init__(*args, **kw)
		self.itemlist = super(CherenkovHitsBins ,self).keys()
	def __setitem__(self, key, value):
		if key not in self:
			self.itemlist.append(key)
			super(CherenkovHitsBins, self).__setitem__(key, [value])
		else:
			self[key].append(value)
	def __iter__(self):
		return iter(self.itemlist)
	def keys(self):
		return self.itemlist
	def values(self):
		 return [self[key] for key in self]
	def itervalues(self):
		return (self[key] for key in self)

# Read Geometry
i3geoFile = I3File(i3GeoLocation, 'r')
i3geo = None
for frame in i3geoFile:
	# geometry frame
	if frame.Stop.id is 'G':
		i3geo = frame['I3Geometry']
i3geoFile.close()
if i3geo == None:
	print("Geometry File Not Loaded...\nExiting....")
	exit()


# InIce Detector Cylinder
i3DOMRange =  np.array([[0,0],[0,0],[0,0]])
for omgeo in i3geo.omgeo:
	
	# Ignore Non InIce DOMS
	if omgeo.second().omtype.name != 'IceCube':
		continue

	i3DOMPosition = omgeo.second().position
	i3DOMRange[0][0] = min(i3DOMRange[0][0], i3DOMPosition.x)
	i3DOMRange[0][1] = max(i3DOMRange[0][1], i3DOMPosition.x)
	i3DOMRange[1][0] = min(i3DOMRange[1][0], i3DOMPosition.y)
	i3DOMRange[1][1] = max(i3DOMRange[1][1], i3DOMPosition.y)
	i3DOMRange[2][0] = min(i3DOMRange[2][0], i3DOMPosition.z)
	i3DOMRange[2][1] = max(i3DOMRange[2][1], i3DOMPosition.z)

i3cylinderCenter = I3Position(np.mean(i3DOMRange[0]), np.mean(i3DOMRange[1]), np.mean(i3DOMRange[2]))
cylinderLength = (i3DOMRange[2][1] - i3DOMRange[2][0]) + 2 * detectorPadding
cylinderRadius = np.max(np.abs(np.concatenate([i3DOMRange[0] - i3cylinderCenter.x, i3DOMRange[1] - i3cylinderCenter.y]))) + detectorPadding
detectorCylinder = Cylinder(cylinderLength, cylinderRadius, i3cylinderCenter)


intensities  = []
for fileID in range(0, 1):

	print(fileID)
	i3fileLocation = '/data/icecube/atmmuon/eff100/0000000-0000999/Level2_eff100_IC86-2020_corsika_21269_0000{}.i3.zst'.format(fileID)
	
	# Read I3File
	i3file = I3File(i3fileLocation, 'r')
	
	for n, frame in enumerate(i3file):

		# Physics Frame
		if frame.Stop.id is 'P':
		
			# Only Anaylze Second Physics Frame
			if not frame.Has('I3TriggerHierarchy'):
				continue
			
			# Only One Track Passes Through Detector Cylinder	
			MMCTrackList = frame['MMCTrackList']
			i3particleList = [MMCTrack.particle for MMCTrack in MMCTrackList]
			particleMask = [not np.isnan(detectorCylinder.intersection(i3particle.pos, i3particle.dir).first) for i3particle in i3particleList]		
			if sum(particleMask) != 1:
					continue
			i3particle = i3particleList[[i for i, x in enumerate(particleMask) if x][0]]
			
			print(n + 4, len(i3particleList))
			continue
			
			# Get DOM Hits
			i3pulseSeriesMap = frame['InIcePulses'].apply(frame)
			
			cherenkovHitsBins = CherenkovHitsBins()
			for omgeo in i3geo.omgeo:
				
				# DOM Can Be Hit By Cherenkov Radiation
				i3DOMPosition = omgeo.second().position
				cherenkovDistance = I3Calculator.cherenkov_distance(i3particle, i3DOMPosition)
				if cherenkovDistance != cherenkovDistance:
					continue
					
				# Limit Cherenkov Distance maybe no lower limit
				if cherenkovDistance > 100 or cherenkovDistance < 60:
					continue
				
				# Get Track Bin
				i3cherenkovPosition = I3Calculator.cherenkov_position(i3particle, i3DOMPosition)
				trackDistance = (i3cherenkovPosition - i3particle.pos).magnitude			
				trackBin = int(trackDistance/trackBinSize)			
				
				# Get DOM Pulses
				i3pulseSeries = i3pulseSeriesMap.get(omgeo.first())
				
				# No DOM Pulses
				if i3pulseSeries == None:
					cherenkovHitsBins[trackBin] = CherenkovHit(0, cherenkovDistance, trackDistance)
					continue
				
				chargeSum = 0
				for i3pulse in i3pulseSeries:
					# Limit Cherenkov Time Residual (not nessesary)
					#timeResidual = I3Calculator.time_residual(i3particle, i3DOMPosition, i3pulse.time)
					#if timeResidual < -1000 or timeResidual > 1000:
					#	continue
					chargeSum += i3pulse.charge
				cherenkovHitsBins[trackBin] = CherenkovHit(chargeSum, cherenkovDistance, trackDistance)
			
			for bins, hits in cherenkovHitsBins.items():
				meanIntensity = np.mean([hit.intensity() for hit in hits])
				if meanIntensity != 0:
					intensities.append(meanIntensity)

				
	i3file.close()


import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

plt.title('Intensity Histogram')
plt.xlabel("Intensity")
plt.ylabel("Frequency")
plt.yscale('log',base=10) 
plt.hist(intensities, bins='auto')
plt.savefig('intensity_histogram.png')
plt.clf()
