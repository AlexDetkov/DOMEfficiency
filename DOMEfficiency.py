#!/usr/bin/python
import pickle
import operator
import numpy as np

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

#################################################################################################################
# Variables													#
#################################################################################################################

trackFile = 'tracks100.pickle'

#################################################################################################################
# Load Tracks													#
#################################################################################################################

tracks = None
with open(trackFile, 'rb') as trackFile:
    tracks = pickle.load(trackFile) 

#################################################################################################################
# Analyze Tracks												#
#################################################################################################################

DOMsPerBin = 5
lossPerUnitLengthCutoff = 0.015

binOverAvgCutoffs = []
noLossIntensities = []
editedIntensities = []

for binOverAvgIntensityCutoff in np.arange(1, 10, 0.2):

	print("Bin Cutoff: {}".format(binOverAvgIntensityCutoff))
	noLossIntensityPerUnitLengthList = []
	editedIntensityPerUnitLengthList = []

	for track in tracks:
		intensitySum = sum([emission.intensity for emission in track.emissions])
		intensityPerUnitLength = intensitySum/track.length
		
		# no stocastic losses tracks
		#if len(track.losses) == 0:
		#	noLossIntensityPerUnitLengthList.append(intensityPerUnitLength)
			
		# small stocastic loss
		lossPerUnitLength = sum([loss.energy for loss in track.losses])/track.length	
		if lossPerUnitLength <= lossPerUnitLengthCutoff:
			noLossIntensityPerUnitLengthList.append(intensityPerUnitLength)
			continue
			
		if len(track.emissions) < 4 * DOMsPerBin or intensityPerUnitLength == 0:
			continue
	
			
		# sort emissions by distanceOnTrack
		sortedEmissions = sorted(track.emissions, key=operator.attrgetter('distanceOnTrack'))
		
		intensityPerUnitLengthInMaxBin = 0
		intensitySumInMaxBin = None
		binLengthInMaxBin = None		

		startIndex = 0
		lastIndex = len(sortedEmissions) - DOMsPerBin
		while startIndex <= lastIndex:
			endIndex = startIndex + DOMsPerBin - 1
			
			# get emissions in bin
			emissionsInBin = sortedEmissions[startIndex:(endIndex+1)]
			
			# get intensity in bin
			intensityInBin = sum([emission.intensity for emission in emissionsInBin])
			binStart = emissionsInBin[0].distanceOnTrack
			binEnd = emissionsInBin[-1].distanceOnTrack
			binLength = binEnd-binStart
			intensityInBinPerUnitLength = intensityInBin/binLength
				
			# possibly high loss bin
			if intensityInBinPerUnitLength > intensityPerUnitLengthInMaxBin:
				intensityPerUnitLengthInMaxBin = intensityInBinPerUnitLength
				intensitySumInMaxBin = intensityInBin
				binLengthInMaxBin = binLength					
				
			startIndex += 1
			
		editedTrackLength = track.length
		editedIntensitySum = intensitySum
		if intensityPerUnitLengthInMaxBin/intensityPerUnitLength > binOverAvgIntensityCutoff:
			editedIntensitySum -= intensitySumInMaxBin
			editedTrackLength -= binLengthInMaxBin
		
		editedIntensityPerUnitLength = editedIntensitySum/editedTrackLength

		if 0 < editedIntensityPerUnitLength:
			editedIntensityPerUnitLengthList.append(editedIntensityPerUnitLength)
		

	plt.title('Intensity Histogram')
	plt.xlabel("Intensity Per Unit Length")
	plt.ylabel("Frequency")
	plt.hist(editedIntensityPerUnitLengthList, bins=40, alpha=0.5, label="Edited Intensities")
	plt.hist(noLossIntensityPerUnitLengthList, bins=40, alpha=0.5, label="Intensities")
	plt.legend()
	plt.yscale('log', base=10)
	plt.savefig('Histogram/Log/{}.png'.format(binOverAvgIntensityCutoff))
	plt.clf()

	plt.title('Intensity Histogram')
	plt.xlabel("Intensity Per Unit Length")
	plt.ylabel("Frequency")
	plt.hist(editedIntensityPerUnitLengthList, bins=40, alpha=0.5, label="Edited Intensities")
	plt.hist(noLossIntensityPerUnitLengthList, bins=40, alpha=0.5, label="Intensities")
	plt.legend()
	plt.savefig('Histogram/Normal/{}.png'.format(binOverAvgIntensityCutoff))
	plt.clf()
		
	plt.title('Intensity Histogram')
	plt.xlabel("Intensity Per Unit Length")
	plt.ylabel("Frequency")
	plt.hist(editedIntensityPerUnitLengthList, bins=40, normed=True, alpha=0.5, label="Edited Intensities")
	plt.hist(noLossIntensityPerUnitLengthList, bins=40, normed=True, alpha=0.5, label="Intensities")
	plt.legend()
	plt.savefig('Histogram/Normed/{}.png'.format(binOverAvgIntensityCutoff))
	plt.clf()
	
	binOverAvgCutoffs.append(binOverAvgIntensityCutoff)
	editedIntensities.append(np.mean(editedIntensityPerUnitLengthList))
	noLossIntensities.append(np.mean(noLossIntensityPerUnitLengthList))

plt.title('Mean Intensity Per Unit Length vs Bin Over Avgerage Intensity Cutoff')
plt.xlabel("Bin Over Avgerage Intensity Cutoff")
plt.ylabel("Mean Intensity Per Unit Length")
plt.plot(binOverAvgCutoffs, editedIntensities, label="Mean Intensities of Edited Tracks")
plt.plot(binOverAvgCutoffs, noLossIntensities, label="Low Stochastic Loss Mean Intensities (<{:.3f})".format(lossPerUnitLengthCutoff))
plt.legend()
plt.savefig('Scatter/Intensity vs Cutoff.png'.format(binOverAvgIntensityCutoff), dpi=600)
