#!/usr/bin/python
import numpy as np
from icecube.dataio import I3File
from icecube.dataclasses import I3Position
from icecube.phys_services import Cylinder
from icecube.topsimulator import intersect_cylinder

#################################################################################################################
# Track														#
#################################################################################################################

class Track():
	def __init__(self, length):
		self.length = length
		self.emissions = []
		self.losses = []
		
class TrackEvent():
	def __init__(self, distanceOnTrack):
		self.distanceOnTrack = distanceOnTrack

class CherenkovEmission(TrackEvent):
	def __init__(self, distanceOnTrack, intensity):
		TrackEvent.__init__(self, distanceOnTrack)
		self.intensity = intensity

class StochasticLoss(TrackEvent):
	def __init__(self, distanceOnTrack, energy):
		TrackEvent.__init__(self, distanceOnTrack)
		self.energy = energy

#################################################################################################################
# Geometry													#
#################################################################################################################

class I3GeometryFile(I3File):
	def __init__(self, fileLocation):
		I3File.__init__(self, fileLocation, 'r')
		self.i3geo = None
		for frame in self:
			# geometry frame
			if frame.Stop.id is 'G':
				self.i3geo = frame['I3Geometry']
		
		
	def getI3geometry(self):
		return self.i3geo
		
	def getDetectorCylinder(self):
		i3DOMRange =  np.array([[0,0],[0,0],[0,0]])
		for omgeo in self.i3geo.omgeo:
			
			# ignore non InIce DOMS
			if omgeo.second().omtype.name != 'IceCube':
				continue

			i3DOMPosition = omgeo.second().position
			i3DOMRange[0][0] = min(i3DOMRange[0][0], i3DOMPosition.x)
			i3DOMRange[0][1] = max(i3DOMRange[0][1], i3DOMPosition.x)
			i3DOMRange[1][0] = min(i3DOMRange[1][0], i3DOMPosition.y)
			i3DOMRange[1][1] = max(i3DOMRange[1][1], i3DOMPosition.y)
			i3DOMRange[2][0] = min(i3DOMRange[2][0], i3DOMPosition.z)
			i3DOMRange[2][1] = max(i3DOMRange[2][1], i3DOMPosition.z)

		i3center = I3Position(np.mean(i3DOMRange[0]), np.mean(i3DOMRange[1]), np.mean(i3DOMRange[2]))
		length = (i3DOMRange[2][1] - i3DOMRange[2][0])
		radius = np.max(np.abs(np.concatenate([i3DOMRange[0] - i3center.x, i3DOMRange[1] - i3center.y])))
		return Cylinder(length, radius, i3center)
	
detectorPadding = 100
def I3ParticleIntersectsCylinder(i3particle, cylinder):
	return intersect_cylinder(cylinder.center, cylinder.length + 2 * detectorPadding, cylinder.radius + detectorPadding, i3particle, i3particle)
	
dustLayerZRange = [-150, -50]
def CherenkovPassesThroughDustLayer(i3cherenkovEmissionPosition, i3DOMPosition):
	if (dustLayerZRange[0] < i3cherenkovEmissionPosition.z < dustLayerZRange[1]) or \
	(dustLayerZRange[0] < i3DOMPosition.z < dustLayerZRange[1]) or \
	((min(i3cherenkovEmissionPosition.z, i3DOMPosition.z) < dustLayerZRange[0]) and (dustLayerZRange[1] < max(i3cherenkovEmissionPosition.z, i3DOMPosition.z))):
		return True
	return False

#################################################################################################################
# DOM														#
#################################################################################################################

def DOMAngleAcceptance(x): # angle in degrees
	return 1.0 - 3.59e-3 * x + 5.11e-5 * x**2 - 4.27e-6 * x**3 + 5.56e-8 * x**4 -2.73e-10 * x**5 + 4.76e-13 * x**6
	
stochasticLength = 100
def Intensity(charge, cherenkovDistance):
	return charge * np.exp(cherenkovDistance/stochasticLength) * cherenkovDistance


	
