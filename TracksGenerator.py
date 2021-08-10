#!/usr/bin/python
import sys, pickle
from os import listdir, path
from icepy import *
from icecube.dataio import I3File
from icecube import simclasses, recclasses
from icecube.phys_services import I3Calculator

#################################################################################################################
# Variables													#
#################################################################################################################

trackFile = 'tracks100.pickle'
i3filesDirectory = '/data/icecube/domeff_analysis/reco_sim_nominal/allMMCTree'
i3geometryFile = '/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE.i3.gz'

trackEndPadding = -50

#################################################################################################################
# Geometry													#
#################################################################################################################

i3geoFile = I3GeometryFile(i3geometryFile)
i3geo = i3geoFile.getI3geometry()
detectorCylinder = i3geoFile.getDetectorCylinder()
i3geoFile.close()

#################################################################################################################
# Generate Tracks												#
#################################################################################################################
tracks = []

totalFileCount = len(listdir(i3filesDirectory))
for fileCount, i3fileName in enumerate(listdir(i3filesDirectory)):
	print 'File Progress:{:03d}%\r'.format((fileCount * 100)/totalFileCount),
	sys.stdout.flush()
		
	i3fileLocation = path.join(i3filesDirectory, i3fileName)
	i3file = I3File(i3fileLocation, 'r')
	
	for frameCount, frame in enumerate(i3file):
		# physics frame
		if frame.Stop.id is 'P':

			# only anaylze second physics frame
			if not frame.Has('MMCTrackList') or not frame.Has('I3MCTree') or not frame.Has('SplitInIcePulses'):
				continue
				
			MMCTrackList = frame['MMCTrackList']
			particleMask = [I3ParticleIntersectsCylinder(MMCTrack.particle, detectorCylinder) for MMCTrack in MMCTrackList]
			
			# only one track passes detector cylinder
			if sum(particleMask) != 1:
				continue

			MMCTrack = np.array(MMCTrackList)[particleMask][0]
			i3particle = MMCTrack.particle
			
			# begining, end, and length of track
			i3trackStartPosition = I3Position(MMCTrack.xi, MMCTrack.yi, MMCTrack.zi)
			i3trackEndPosition = I3Position(MMCTrack.xf, MMCTrack.yf, MMCTrack.zf)		
			trackLength = (i3trackEndPosition - i3trackStartPosition).magnitude
			
			# get daughters and determine high stochastic loss bins
			I3MCTree = frame['I3MCTree']
			i3daughters = I3MCTree.get_daughters(i3particle)
			
			track = Track(trackLength)
			
			for i3daughterParticle in i3daughters:
			
				# find energy losses
				if i3daughterParticle.type != i3particle.type:

					daughterTrackDistance = (i3daughterParticle.pos - i3trackStartPosition).magnitude
					
					if daughterTrackDistance > (trackLength + trackEndPadding):
						continue
					
					track.losses.append(StochasticLoss(daughterTrackDistance, i3daughterParticle.energy))
			
			# get DOM hits
			i3pulseSeriesMap = frame['SplitInIcePulses'].apply(frame)
			
			for omgeo in i3geo.omgeo:
				
				i3DOMPosition = omgeo.second().position
				cherenkovDistance = I3Calculator.cherenkov_distance(i3particle, i3DOMPosition)
				
				# DOM can be hit by cherenkov radiation
				if cherenkovDistance != cherenkovDistance:
					continue
					
				# limit cherenkov distance
				if cherenkovDistance > 100 or cherenkovDistance < 20:
					continue
				
				# get emission track distance
				i3cherenkovEmissionPosition = I3Calculator.cherenkov_position(i3particle, i3DOMPosition)
				emissionTrackDistance = (i3cherenkovEmissionPosition - i3trackStartPosition).magnitude
				
				# cherenkov distance is within track
				if emissionTrackDistance > (trackLength + trackEndPadding):
					continue
				
				# light doesnt pass through dust layer
				if CherenkovPassesThroughDustLayer(i3cherenkovEmissionPosition, i3DOMPosition):
					continue
				
				# get DOM impact angle and angle acceptance
				DOMImpactAngleDeg = np.rad2deg(I3Calculator.cherenkov_approach_angle(i3particle, i3DOMPosition, omgeo.second().direction))
				
				# filter impact angle
				if DOMImpactAngleDeg > 135:
					continue
				
				# get DOM pulses
				i3pulseSeries = i3pulseSeriesMap.get(omgeo.first())
				
				# no DOM Pulses
				if i3pulseSeries == None:
					track.emissions.append(CherenkovEmission(emissionTrackDistance, 0)) 
					continue
				
				# sum DOM charges
				chargeSum = 0
				for i3pulse in i3pulseSeries:
					timeResidual = I3Calculator.time_residual(i3particle, i3DOMPosition, i3pulse.time)
					if -100 < timeResidual < 100:
						chargeSum += i3pulse.charge			
				
				# fix charge based on acceptance
				DOMacceptance = DOMAngleAcceptance(DOMImpactAngleDeg)
				chargeSum /= DOMacceptance
				
				track.emissions.append(CherenkovEmission(emissionTrackDistance, Intensity(chargeSum, cherenkovDistance)))
				
			tracks.append(track)
	i3file.close()

print 'File Progress:100%'

#################################################################################################################
# Save Tracks													#
#################################################################################################################

print 'Saving Tracks in {}'.format(trackFile)
with open(trackFile, 'wb') as pickleFile:
    pickle.dump(tracks, pickleFile) 
