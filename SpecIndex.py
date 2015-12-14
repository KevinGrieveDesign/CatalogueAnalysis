#!/usr/bin/python
#misc packages
import os
import time
import csv
import argparse
import subprocess
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, FileTransferSpeed, FormatLabel, Percentage, ProgressBar, ReverseBar, RotatingMarker, SimpleProgress, Timer  
from datetime import datetime

#import the plotting tools
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


#import the fitting tools
import numpy as np
from numpy.polynomial import polynomial as P
from kapteyn import kmpfit

startTime = datetime.now()

def alpha(p, x):
	a, b = p
	return pow( 10, a*np.log10(x) + b )

def Curve(p, x):
	a,b,c = p
   	return a * np.sin(b*x+c)

def obj_printToScreen(fitobj, DictItem):
	print "\nFit status of " + DictItem + ":"
	print "===================="
	print "Best-fit parameters:        ", fitobj.params
	print "Asymptotic error:           ", fitobj.xerror
	print "Error assuming red.chi^2=1: ", fitobj.stderr
	print "Chi^2 min:                  ", fitobj.chi2_min
	print "Reduced Chi^2:              ", fitobj.rchi2_min
	print "Iterations:                 ", fitobj.niter
	print "Number of free pars.:       ", fitobj.nfree
	print "Degrees of freedom:         ", fitobj.dof

def obj_printToFile(fitobj, SourceName, Frequencies, Row, DestinationPath, LogFileName):
	Row = str(Row).replace('[','')
	Row = str(Row).replace(']','')
	Row = str(Row).replace('\'','')

	Frequencies = str(Frequencies).replace('[','')
	Frequencies = str(Frequencies).replace(']','')
	Frequencies = str(Frequencies).replace(',',' -')

	CSV_Line = ""
	CSV_Line = CSV_Line + str(SourceName) + ","
	CSV_Line = CSV_Line + str(Frequencies) + ","
	CSV_Line = CSV_Line + str(fitobj.params[0]) + "," + str(fitobj.params[1]) + ","
	CSV_Line = CSV_Line + str(fitobj.xerror[0]) + "," + str(fitobj.xerror[1]) + ","
	CSV_Line = CSV_Line + str(fitobj.stderr[0]) + "," + str(fitobj.stderr[1]) + ","
	CSV_Line = CSV_Line + str(fitobj.chi2_min) + ","
	CSV_Line = CSV_Line + str(fitobj.rchi2_min) + ","
	CSV_Line = CSV_Line + str(fitobj.niter) + ","
	CSV_Line = CSV_Line + str(fitobj.nfree) + ","
	CSV_Line = CSV_Line + str(fitobj.dof) 

	with open(DestinationPath + "/" + LogFileName + "StatsSimple.csv", "a") as myfile:
		myfile.write(CSV_Line + "\n")

	CSV_Line = CSV_Line + "," + str(Row) + "\n"

	with open(DestinationPath + "/" + LogFileName + "Stats.csv", "a") as myfile:
		myfile.write(CSV_Line)

def Headings_printToFile( Headings, DestinationPath, LogFileName):
	Headings = str(Headings).replace('[','')
	Headings = str(Headings).replace(']','')
	Headings = str(Headings).replace('\'','')

	CSV_Line = ""
	CSV_Line = CSV_Line + "Name,"
	CSV_Line = CSV_Line + "Frequencies,"
	CSV_Line = CSV_Line + "params 0,params 1,"
	CSV_Line = CSV_Line + "xerror 0,xerror 1,"
	CSV_Line = CSV_Line + "stderr 0,stderr 1,"
	CSV_Line = CSV_Line + "chi2_min,"
	CSV_Line = CSV_Line + "rchi2_min,"
	CSV_Line = CSV_Line + "niter,"
	CSV_Line = CSV_Line + "nfree,"
	CSV_Line = CSV_Line + "dof"

	with open(DestinationPath + "/" + LogFileName + "StatsSimple.csv", "w") as myfile:
		myfile.write(CSV_Line + "\n")

	CSV_Line = CSV_Line + "," + str(Headings) + "\n"

	with open(DestinationPath + "/" + LogFileName + "Stats.csv", "w") as myfile:
		myfile.write(CSV_Line)


def obj_printForTex(fitobj, SourceName, Frequencies, Row, DestinationPath, LogFileName, count, firstTime):
	Row = str(Row).replace('[','')
	Row = str(Row).replace(']','')
	Row = str(Row).replace('\'','')

	Frequencies = str(Frequencies).replace('[','')
	Frequencies = str(Frequencies).replace(']','')
	Frequencies = str(Frequencies).replace(',',' -')

	TexLine = ""

	if count == 1 or count == 3:
		#TexLine = TexLine + "\\begin{multicols}{2} "
		TexLine = TexLine + "\\begin{figure}[H] " 

		if firstTime == True and count == 1:
			TexLine = TexLine + "\\section{Spectral Index} "
			
	TexLine = TexLine + "\\includegraphics[scale=0.50]{Images/" + LogFileName + "/" + str(SourceName) + ".eps} "

	if count == 4 and firstTime == True:
		TexLine = TexLine + "\\caption{Spectral Index Plots for all sources that contain three or more detections.} "
 
	# TexLine = TexLine + "\\begingroup" 
	# TexLine = TexLine + "\\centering"
	# TexLine = TexLine + "\\includegraphics[scale=0.2]{Images/" + LogFileName + "/" + str(SourceName) + ".eps}"
	# TexLine = TexLine + "\\captionof{Spectral Index of source " + str(SourceName) + ", using " + str(Frequencies) + " with an \\ac{SI} of " + str(fitobj.params[0]) + "}"
	# TexLine = TexLine + "\\endgroup"

	if count % 2 == 0:
		TexLine = TexLine + "\\end{figure} "
		#TexLine = TexLine + "\\end{multicols} "

	if count == 4:
		if firstTime == False:
			TexLine = TexLine + "\\centering "	
			TexLine = TexLine + "Figure 4.1: \\textit{Continued}. "

		TexLine = TexLine + "\\clearpage "


	# TexLine = "\\begin{figure}[bht]" 
	# TexLine = TexLine + "\\centering"
	# TexLine = TexLine + "\\includegraphics[scale=0.2]{Images/" + LogFileName + "/" + str(SourceName) + ".eps}"
	# TexLine = TexLine + "\\caption{Spectral Index of source " + str(SourceName) + ", using " + str(Frequencies) + " with an \\ac{SI} of " + str(fitobj.params[0]) + "}"
	# TexLine = TexLine + "\\end{figure}"

	# CSV_Line = ""
	# CSV_Line = CSV_Line + str(SourceName) + ","
	# CSV_Line = CSV_Line + str(Frequencies) + ","
	# CSV_Line = CSV_Line + str(fitobj.params[0]) + "," + str(fitobj.params[1]) + ","
	# CSV_Line = CSV_Line + str(fitobj.xerror[0]) + "," + str(fitobj.xerror[1]) + ","
	# CSV_Line = CSV_Line + str(fitobj.stderr[0]) + "," + str(fitobj.stderr[1]) + ","
	# CSV_Line = CSV_Line + str(fitobj.chi2_min) + ","
	# CSV_Line = CSV_Line + str(fitobj.rchi2_min) + ","
	# CSV_Line = CSV_Line + str(fitobj.niter) + ","
	# CSV_Line = CSV_Line + str(fitobj.nfree) + ","
	# CSV_Line = CSV_Line + str(fitobj.dof) 

	with open(DestinationPath + "/" + LogFileName + "SpecIndex.tex", "a") as myfile:
		myfile.write(TexLine + "\n")

# The notifier function
def notify(title, subtitle, message):
    t = '-title {!r}'.format(title)
    s = '-subtitle {!r}'.format(subtitle)
    m = '-message {!r}'.format(message)
    a = '-sound \'Glass\''
    os.system('terminal-notifier {}'.format(' '.join([m, t, s, a])))

#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================

parser = argparse.ArgumentParser(description='This will generate spectral index plots from the catalogue')
parser.add_argument('-F' , '-f', '--file', help='This is the name of the Catalogue you wish to read from. Default: LMC.point.csv')
parser.add_argument('-L' , '-l', '--log', help='Write to the Stats.csv and SimpleStats.csv files. Default is true.')
args = parser.parse_args()

FileName = ""
# DestinationPath = "Plots/ExtendedOptical_SpecIndex/"
# DestinationPath = "Plots/CompactOptical_SpecIndex/"
DestinationPath = "Plots/LMC42/SpecIndex/"

Data = {} 
Style = {}
firstTime = True

ConditionalColumns = {}
ConditionalColumns['Good Source'] = "true"
ConditionalColumns['Extended Optical Source (MCELS)'] = "false"
#ConditionalColumns['Compact Optical Source at Centre'] = "true"
ConditionalColumns['Blended'] = "false"
#ConditionalColumns['KnownObject'] = ""

#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================

Style['AT20G Integ_20GHz'] = {
'Name'   : 'ATCA, AT20G - 20 GHz',
'xUpperLimit' : 35000,
'Colour' : 'orange',
'Shape'  : 'o'}

Style['AT20G Integ_8GHz'] = {
'Name'   : 'ATCA, AT20G - 8 GHz',
'xUpperLimit' : 35000,
'Colour' : 'orange',
'Shape'  : 'o'}

Style['AT20G Integ_5GHz'] = {
'Name'   : 'ATCA, AT20G - 5 GHz',
'xUpperLimit' : 35000,
'Colour' : 'orange',
'Shape'  : 'o'}

Style['3cm Integ']  = {
'Name'   : 'ATCA, 8.64 GHz',
'xUpperLimit' : 13000,
'Colour' : 'yellow',
'Shape'  : 'o'}

Style['PMN Flux']  = {
'Name'   : 'PMN, 4.85 GHz',
'xLowerLimit' : 1000,
'xUpperLimit' : 11000,
'Colour' : 'blue',
'Shape'  : 'o'}

Style['6cm Integ']  = {
'Name'   : 'ATCA, 4.80 GHz',
'xLowerLimit' : 1000,
'xUpperLimit' : 11000,
'Colour' : 'yellow',
'Shape'  : 'o'}

Style['Marx 1.4GHz Integ']  = {
'Name'   : 'ATCA, 1.42 GHz',
'xLowerLimit' : 600,
'xUpperLimit' : 5000,
'Colour' : 'red',
'Shape'  : 'o'}

Style['20cm Integ']  = {
'Name'   : 'ATCA, 1.42 GHz',
'xLowerLimit' : 600,
'xUpperLimit' : 5000,
'Colour' : 'green',
'Shape'  : 'o'}

Style['Marx 2.4GHz Integ']  = {
'Name'   : 'ATCA, 1.42 GHz',
'xLowerLimit' : 600,
'xUpperLimit' : 5000,
'Colour' : 'red',
'Shape'  : 'o'}

Style['36cm Integ']   = {
'Name'   : 'UTMOST, 0.843 GHz',
'xLowerLimit' : 550,
'xUpperLimit' : 5000,
'Colour' : 'green',
'Shape'  : 'o'}

Style['SUMSS Integ']   = {
'Name'   : 'SUMSS - 0.843 GHz',
'xLowerLimit' : 550,
'xUpperLimit' : 5000,
'Colour' : 'blue',
'Shape'  : 'o'}

Style['MRC 408 Flux']   = {
'Name'   : 'MRC - 0.408 GHz',
'xLowerLimit' : 250,
'Colour' : 'blue',
'Shape'  : 'o'}

Style['S408']   = {
'Name'   : 'MRC - 0.408 GHz',
'xLowerLimit' : 250,
'Colour' : 'blue',
'Shape'  : 'o'}





#'xCol' : ['19904', '8640', '4800-4850', '1384', '843-843', '408'],
#'xColLabel' : ['20000', '8640', '4800-4850', '1384', '843-843', '408'],

#'yCol' : ['AT20G Integ_20GHz', '3cm Integ', '6cm Integ-PMN Flux', '20cm Integ', '36cm Integ-SUMSS Integ', 'MRC 408 Flux'],
#'yColErr' : ['AT20G Integ_20GHz_Err', '3cm Integ_Err', '6cm Integ_Err-PMN Flux_Err', '20cm Integ_Err', '36cm Integ_Err-SUMSS Integ_Err', 'MRC 408 Flux_Err'],




#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================


# Data['Public1.1-6Sources-EPS'] = {
# 'xCol' : ['19904', '8640', '8640', '4800-4850', '4800', '1384', '843-843', '408'],
# 'xColLabel' : ['20000', '8640', '8640', '4800-4850', '4800', '1384', '843-843', '408'],

# 'yCol' : ['AT20G Integ_20GHz', '3cm Integ', 'AT20G Integ_8GHz', '6cm Integ-PMN Flux', 'AT20G Integ_5GHz', '20cm Integ', '36cm Integ-SUMSS Integ', 'MRC 408 Flux'],
# 'yColErr' : ['AT20G Integ_20GHz_Err', '3cm Integ_Err', 'AT20G Integ_8GHz_Err', '6cm Integ_Err-PMN Flux_Err', 'AT20G Integ_5GHz_Err', '20cm Integ_Err', '36cm Integ_Err-SUMSS Integ_Err', 'MRC 408 Flux_Err'],

# 'yColErrTol' : 10,

# 'xScale' : 'log',
# 'yScale' : 'log',

# 'Aspect' : 'equal',

# 'MinimumPoints' : 6,
# 'MaximumPoints' : 8,
# 'Type'   : 'scatter'}

Data['MarxTest2'] = {
# 'xCol' : ['19904', '8640', '8640', '4800-4850', '4800', '1384', '843-843', '408'],
# 'xColLabel' : ['20000', '8640', '8640', '4800-4850', '4800', '1384', '843-843', '408'],

'xCol' : ['19904', '8640', '8640', '4800-4850', '4800', '2378', '1384', '1380', '843-843', '408'],
'xColLabel' : ['20000', '8640', '8640', '4800-4850', '4800', '2378', '1384', '1380', '843-843', '408'],

'yCol' : ['AT20G Integ_20GHz', '3cm Integ', 'AT20G Integ_8GHz', '6cm Integ-PMN Flux', 'AT20G Integ_5GHz', 'Marx 2.4GHz Integ', '20cm Integ', 'Marx 1.4GHz Integ', '36cm Integ-SUMSS Integ', 'MRC 408 Flux'],
'yColErr' : ['AT20G Integ_20GHz_Err', '3cm Integ_Err', 'AT20G Integ_8GHz_Err', '6cm Integ_Err-PMN Flux_Err', 'AT20G Integ_5GHz_Err', 'Marx 2.4GHz Integ_Err', '20cm Integ_Err', 'Marx 1.4GHz Integ_Err', '36cm Integ_Err-SUMSS Integ_Err', 'MRC 408 Flux_Err'],

'yColErrTol' : 10,

'xScale' : 'log',
'yScale' : 'log',

'Aspect' : 'equal',

'MinimumPoints' : 3,
'MaximumPoints' : 10,
'Type'   : 'scatter'}

						
#==========================================================================================================================================
#==========================================================================================================================================
#==========================================================================================================================================

if args.file is not None:
	FileName  = args.file
else:
	FileName = "LMC.point.csv"

if args.log is not None:
	LogStats = args.log
else:
	LogStats = True

os.system('mkdir ' + DestinationPath)

#get the size of the csv file 
proc = subprocess.Popen("wc -l < " + FileName, shell=True, stdout=subprocess.PIPE)
LinesInFile = proc.communicate()[0]

maxval = int(LinesInFile) * len(Data)
Count = 0 
ExcludedCount = -1
testcol = 1

#create the progress bar from the size calculated
widgets = ['Creating Spectral Index Plots:', Percentage(), ' (', Counter(), '/' + str(maxval) + ') ', Bar()]
pbar = ProgressBar(widgets = widgets, maxval = maxval, poll=0.1).start()

#loop through the different plots that are desired.
for DictItem in list(Data.keys()):
	os.system('mkdir ' + DestinationPath + DictItem + '/')
	os.system('rm ' + DestinationPath + "/" + DictItem + "SpecIndex.tex")
	ExcludedCount = 0

	#open csv
	with open(FileName,'rb') as Catalogue:
		reader = csv.reader(Catalogue)
		Headings = next(reader)

		if Count == 0 and LogStats == True or str(LogStats) == "true":
			Headings_printToFile(Headings, DestinationPath, DictItem)

		Count += 1
		pbar.update(Count)

		#get the two important columns 
		GoodSource = Headings.index('Good Source')
		Blended = Headings.index('Blended')
	
		#loop through the csv
		for RowNum, Row in enumerate(reader):
			Data['x'] = []
			Data['y'] = [] 
			Data['yErr'] = [] 

			Data['Colour'] = [] 
			Data['Shape'] = [] 
			Data['Name'] = [] 

			ExcludeRow = False
			Count += 1
			pbar.update(Count)

			#check to see if i should exclude the line from the plot
			for Conditions in list(ConditionalColumns):
	 			if Row[Headings.index(Conditions)] != ConditionalColumns[Conditions]:
	 				ExcludeRow = True
	 				ExcludedCount += 1

			if ExcludeRow == False:			
				#loop through the different columns and get the flux data as well as the error on the flux
				for yColNum, yCol in enumerate(Data[DictItem]['yCol']):
					#if there is more than one survey at the same frequency and i only want to use one 
					#then i use a - in one list item to separate them
					if yCol.find('-') > -1:
						ColsFound = 0 

						#split the different surveys and loop through them. the one at the start has preference 
						for yColNumDecision, yColDecision in enumerate(yCol.split('-')):
							#if we have data in that column and havent found something before, then use it
							if bool(Row[Headings.index(yColDecision)]) and ColsFound == 0:
								ColsFound += 1

								#work out what the column name of the x and y_error Columns are. We know the y column as we are looping on it
								xColDecision = Data[DictItem]['xCol'][yColNum].split('-')
								yColDecisionErr = Data[DictItem]['yColErr'][yColNum].split('-')

								Data['x'].append(float(xColDecision[yColNumDecision]))     #x column
								Data['y'].append(float(Row[Headings.index(yColDecision)])) #y column

								#determine if the y error is greater than a pre-defined amount (see top of script), if yes, use predefined val
								if bool(Row[Headings.index(yColDecisionErr[yColNumDecision])]) and float(Row[Headings.index(yColDecisionErr[yColNumDecision])]) > 0 and ( float(Row[Headings.index(yColDecision)]) * float(Data[DictItem]['yColErrTol'])) > float(Row[Headings.index(yColDecisionErr[yColNumDecision])]):
								 	Data['yErr'].append( float(Row[Headings.index(yColDecisionErr[yColNumDecision])]))
								else:
								 	Data['yErr'].append( float(Row[Headings.index(yColDecision)])/float(Data[DictItem]['yColErrTol'])) 
								 	
								Data['Colour'].append(str(Style[yColDecision]['Colour']))
								Data['Shape'].append(str(Style[yColDecision]['Shape']))
								Data['Name'].append(str(yColDecision))

					elif bool(Row[Headings.index(yCol)]):
						#this does the same as above, checks if there is data in the column and uses it if there is. 
						Data['x'].append(float(Data[DictItem]['xCol'][yColNum]))  #x column
						Data['y'].append(float(Row[Headings.index(yCol)]))	 	  #y column

						if bool(Row[Headings.index(Data[DictItem]['yColErr'][yColNum])]) and float(Row[Headings.index(Data[DictItem]['yColErr'][yColNum])]) > 0 and ( float(Row[Headings.index(yCol)]) * float(Data[DictItem]['yColErrTol'])) > float(Row[Headings.index(Data[DictItem]['yColErr'][yColNum])]):
							Data['yErr'].append(float(Row[Headings.index(Data[DictItem]['yColErr'][yColNum])]))
						else:
							Data['yErr'].append(float(Row[Headings.index(yCol)])/float(Data[DictItem]['yColErrTol']))
							
						Data['Colour'].append(str(Style[yCol]['Colour']))
						Data['Shape'].append(str(Style[yCol]['Shape']))
						Data['Name'].append(str(yCol))
							
						#make sure that we have more than 1 data point on the object
				if len(Data['x']) >= Data[DictItem]['MinimumPoints'] and len(Data['x']) <= Data[DictItem]['MaximumPoints'] :
					#===============================================================================
					#============================Setup Plot=========================================
					#===============================================================================
					fontsize=18

					#generate the plot and make it look nice 
					fig = plt.figure()
					fig.set_size_inches(6.5 , 7)
					
					ax = fig.add_subplot(111)

					if 'Aspect' in Data[DictItem]:
						ax.set_aspect(Data[DictItem]['Aspect'])

					#set x the plot scale
					if bool(Data[DictItem]['xScale']):
						ax.set_xscale(Data[DictItem]['xScale'], fontsize=fontsize)

					#set y the plot scale
					if bool(Data[DictItem]['yScale']):
						ax.set_yscale(Data[DictItem]['yScale'], fontsize=fontsize)
					
					# xLowerLimit = Style[Data['Name'][len(Data['Name']) - 1]]['xLowerLimit']
					# xUpperLimit = Style[Data['Name'][0]]['xUpperLimit']

					yLowerLimit = 0
					yUpperLimit = 0

					xLowerLimit = 150
					xUpperLimit = 38000

					yLowerLimit = 0
					yUpperLimit = 0		
			
					plt.xlim(xLowerLimit,xUpperLimit)
				

					#set the plot limit
					# if 'xLowerLimit' in Data[DictItem] and 'xUpperLimit' in Data[DictItem]:
					# 	plt.xlim(Data[DictItem]['xLowerLimit'], Data[DictItem]['xUpperLimit'])
					
					# if 'yLowerLimit' in Data[DictItem] and 'yUpperLimit' in Data[DictItem]:
					# 	plt.ylim(Data[DictItem]['yLowerLimit'], Data[DictItem]['yUpperLimit'])

					plt.ylim(Data['y'][len(Data['y']) - 1] / 100, Data['y'][len(Data['y']) - 1] * 10)



					#===============================================================================
					#============================Start the Fitting==================================
					#===============================================================================

					dx = np.arange(xLowerLimit, xUpperLimit, 0.1)


					# if 'xLowerLimit' in Data[DictItem] and 'xUpperLimit' in Data[DictItem]:
					# 	dx = np.arange(Data[DictItem]['xLowerLimit'], Data[DictItem]['xUpperLimit'] + 500, 0.1)
					# else:
					# 	dx = np.arange(200, 10000, 0.1)

					#do the fitting and plot it		
					if bool(Data[DictItem]['yColErr']): 
						fitting_obj = kmpfit.simplefit(alpha, [ 0, 0 ], Data['x'], Data['y'], err = Data['yErr'])
						# p = np.polyfit(Data['x'], Data['y'], 2, w = Data['yErr'])
					else:
						fitting_obj = kmpfit.simplefit(alpha, [ 0, 0 ], Data['x'], Data['y'])
						# p = np.polyfit(Data['x'], Data['y'], 2)

					# add the fit line
					plt.plot(dx, alpha(fitting_obj.params, dx), c = 'black')

					# # add in a polynomial fit 
					# PolyFit = np.poly1d(p)

					# # show the poly fit
					# x_new = np.linspace(Data['x'][0], Data['x'][-1], 150)
					# y_new = PolyFit(x_new)
					# plt.plot(Data['x'],Data['y'],'o', x_new, y_new, c = 'blue')

					#===============================================================================
					#============================Plot Data Points Plot==============================
					#===============================================================================

					#Plot the data points, creating the error bars and the legend
					for xNum, x in enumerate(Data['x']):
						plt.errorbar(x, Data['y'][xNum], yerr=Data['yErr'][xNum], ls = 'none', c = 'black')
						#l = ax.plot(x, Data['y'][xNum], Data['Shape'][xNum], c = Data['Colour'][xNum], label = Data['Name'][xNum])
						l = ax.plot(x, Data['y'][xNum], Data['Shape'][xNum], c = Data['Colour'][xNum])
					
					#create the labels
					plt.xlabel('Frequency (MHz)', fontsize=fontsize)
					plt.ylabel('Integrated Flux (Jy)', fontsize=fontsize)
					plt.title('Source: ' + str(Row[Headings.index('Name')]), fontsize=fontsize)

					plt.tick_params(axis='both', which='major', labelsize=fontsize)
					plt.tick_params(axis='both', which='minor', labelsize=fontsize)

					ax.text(0.5, 0.01, r'$\alpha =' + str(format(float(fitting_obj.params[0]),'.2f')) + ' \pm' + str(format(float(fitting_obj.stderr[0]),'.2f')) + '$', verticalalignment='bottom', horizontalalignment='right', transform=ax.transAxes, fontsize=fontsize)

					#set legend position
					plt.legend(loc=3, fontsize = 'x-small')
					#plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

					#remove the old plot and save the new one
					plt.savefig(DestinationPath + DictItem + "/" + str(Row[Headings.index('Name')]) + '.eps')
					plt.close()

					if LogStats == True or str(LogStats) == "true":
						obj_printToFile(fitting_obj, str(Row[Headings.index('Name')]), Data['x'], Row, DestinationPath, DictItem)
						obj_printForTex(fitting_obj, str(Row[Headings.index('Name')]), Data['x'], Row, DestinationPath, DictItem, testcol, firstTime)


					if testcol < 4:
						testcol = testcol + 1
					else: 
						testcol = 1
						firstTime = False

					#obj_printToScreen(fitting_obj, str(Row[Headings.index('Name')]))


	# \end{multicols}
	# \clearpage

	LastTexLine = "Figure 4.1: \\textit{Continued.} "


	with open(DestinationPath + "/" + DictItem + "SpecIndex.tex", "a") as myfile:
		myfile.write(LastTexLine + "\n")
				

pbar.finish()


# if testcol < 16:
# 	TexLine = "\\end{multicols}"
# 	TexLine = TexLine + "\\clearpage"

# \end{multicols}
# \clearpage



# 	with open(DestinationPath + "/" + LogFileName + "SpecIndex.tex", "a") as myfile:
# 		myfile.write(TexLine + "\n")


# Calling the function
notify(title    = 'Spectral Index Plots',
       subtitle = 'Plot-SpecIndex.py',
       message  = 'Started: ' + str(startTime)[11:19] + ' - Time Taken: ' + str(datetime.now()-startTime)[:7] )

