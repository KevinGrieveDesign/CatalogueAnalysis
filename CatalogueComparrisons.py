#!/usr/bin/python
import os
import time
import csv
import argparse
import subprocess
import scipy.odr as odr
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from progressbar import ProgressBar    
from datetime import datetime
import numpy as np
from kapteyn import kmpfit

startTime = datetime.now()

def line(p, x):
	m, a = p
	return m*x+a

def lineErr(B, x):
    m, b = B
    return m*x+b

def alpha(p, x):
	a, b = p
	return pow( 10, a*np.log10(x) + b )

def obj_print(fitobj, DictItem):
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

def makeStatsFile(DictItem = False, Data = False, FittingObj = False):
	if DictItem == False:
		CSV_Line = ""
		CSV_Line = CSV_Line + "FileName, Name, xCol, yCol,"
		CSV_Line = CSV_Line + "xMin, xMax, xMean, xMedian, xStdDev,"
		CSV_Line = CSV_Line + "yMin, yMax, yMean, yMedian, yStdDev"

		# CSV_Line = CSV_Line + "xerror 0,xerror 1,"
		# CSV_Line = CSV_Line + "stderr 0,stderr 1,"
		# CSV_Line = CSV_Line + "chi2_min,"
		# CSV_Line = CSV_Line + "rchi2_min,"
		# CSV_Line = CSV_Line + "niter,"
		# CSV_Line = CSV_Line + "nfree,"
		# CSV_Line = CSV_Line + "dof"

		with open("Plots/LMC42/Stats.csv", "w") as myfile:
			myfile.write(CSV_Line + "\n")

	else:

		CSV_Line = ""
		CSV_Line = CSV_Line + str(DictItem) + "," +str(Data[DictItem]['Title']) + "," + Data[DictItem]['xLabel'] + "," + Data[DictItem]['yLabel'] + ","

		CSV_Line = CSV_Line + str(np.min(Data['x'])) + "," + str(np.max(Data['x'])) + "," 
		CSV_Line = CSV_Line + str(np.mean(Data['x'])) + "," + str(np.median(Data['x'])) + "," + str(np.std(Data['x'])) + "," 

		CSV_Line = CSV_Line + str(np.min(Data['y'])) + "," + str(np.max(Data['y'])) + "," 
		CSV_Line = CSV_Line + str(np.mean(Data['y'])) + "," + str(np.median(Data['y'])) + "," + str(np.std(Data['y'])) 
	
		with open("Plots/LMC42/Stats.csv", "a") as myfile:
			myfile.write(CSV_Line + "\n")

		# print Data['y']
		# print Data['xErr']
		# print Data['yErr']

		# if bool(FittingObj):
		# 	print str(FittingObj.beta)
		# 	print str(FittingObj.pprint())


# The notifier function
def notify(title, subtitle, message):
    t = '-title {!r}'.format(title)
    s = '-subtitle {!r}'.format(subtitle)
    m = '-message {!r}'.format(message)
    a = '-sound \'Glass\''
    os.system('terminal-notifier {}'.format(' '.join([m, t, s, a])))

parser = argparse.ArgumentParser(description='Make many region files from the catalogue')
parser.add_argument('-F' , '-f', '--file', help='This is the name of the Catalogue you wish to read from. Default: LMC.point.csv')
args = parser.parse_args()

FileName = ""
ExcludedCount = -1

if args.file is not None:
	FileName  = args.file
else:
	FileName = "LMC.point.csv"

ConditionalColumns = {}
ConditionalColumns['Good Source'] = "true"
ConditionalColumns['Extended Optical Source (MCELS)'] = "false"
#ConditionalColumns['Compact Optical Source at Centre'] = "true"
ConditionalColumns['Blended'] = "false"
ConditionalColumns['KnownObject'] = "false"

Data = {} 


#===========================================================================================================================================
#========================================================  Integ Flux   ====================================================================
#===========================================================================================================================================

Data['Integ-line-Err-6cm-PMN-03'] = {
'Title' : 'Integrated Flux - ATCA 4.8 GHz vs PMN',

'xCol' : '6cm Integ',
'yCol' : 'PMN Flux',
'xColErr' : '6cm Integ_Err',
'yColErr' : 'PMN Flux_Err',

'xLabel' : 'ATCA 4.8 GHz (Jy)',
'yLabel' : 'PMN (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : 0,
'UpperLimit' : 0.3,

'1-1_Line' : True,
'0_Line' : True,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr'}

Data['Integ-line-Err-6cm-PMN-08'] = {
'Title' : 'Integrated Flux - ATCA 4.8 GHz vs PMN',

'xCol' : '6cm Integ',
'yCol' : 'PMN Flux',
'xColErr' : '6cm Integ_Err',
'yColErr' : 'PMN Flux_Err',

'xLabel' : 'ATCA 4.8 GHz (Jy)',
'yLabel' : 'PMN (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : 0,
'UpperLimit' : 0.8,

'1-1_Line' : True,
'0_Line' : True,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr'}

Data['Integ-line-Err-20cm-Marx-1.4GHz'] = {
'Title' : 'Integrated Flux - 20cm vs 1.4GHz Marx',

'xCol' : '20cm Integ',
'yCol' : 'Marx 1.4GHz Integ',

'xColErr' : '20cm Integ_Err',

'xLabel' : '20cm (Jy)',
'yLabel' : 'Marx 1.4GHz (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : 0,
'UpperLimit' : 0.3,

'1-1_Line' : True,
'0_Line' : False,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr'}

Data['Integ-line-Err-20cm-Marx-2.4GHz'] = {
'Title' : 'Integrated Flux - 20cm vs 2.4GHz Marx',

'xCol' : '20cm Integ',
'yCol' : 'Marx 2.4GHz Integ',

'xColErr' : '20cm Integ_Err',

'xLabel' : '20cm (Jy)',
'yLabel' : 'Marx 2.4GHz (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : 0,
'UpperLimit' : 0.3,

'1-1_Line' : True,
'0_Line' : False,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr'}

Data['Integ-line-Err-36cm-SUMSS'] = {
'Title' : 'Integrated Flux - MOST 0.843 GHz vs SUMSS',

'xCol' : '36cm Integ',
'yCol' : 'SUMSS Integ',
'xColErr' : '36cm Integ_Err',
'yColErr' : 'SUMSS Integ_Err',

'xLabel' : 'MOST 0.843 GHz (Jy)',
'yLabel' : 'SUMSS (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : True,
'LowerLimit' : 0,
'UpperLimit' : 1.5,

'1-1_Line' : True,
'0_Line' : True,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr',
'Exclude' : ['043856-672153','053543-660204','051537-672128','052502-693840']}

# #===========================================================================================================================================
# #========================================================  Peak Flux   =====================================================================
# #===========================================================================================================================================

Data['Peak-line-Err-6cm-PMN-03'] = {
'Title' : 'Peak Flux - ATCA 4.8 GHz vs PMN',

'xCol' : '6cm Peak',
'yCol' : 'PMN Flux',
'xColErr' : '6cm Peak_Err',
'yColErr' : 'PMN Flux_Err',

'xLabel' : 'ATCA 4.8 GHz (Jy)',
'yLabel' : 'PMN (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : 0,
'UpperLimit' : 0.3,

'1-1_Line' : True,
'0_Line' : True,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr'}

Data['Peak-line-Err-6cm-PMN-08'] = {
'Title' : 'Peak Flux - ATCA 4.8 GHz vs PMN',

'xCol' : '6cm Peak',
'yCol' : 'PMN Flux',
'xColErr' : '6cm Peak_Err',
'yColErr' : 'PMN Flux_Err',

'xLabel' : 'ATCA 4.8 GHz (Jy)',
'yLabel' : 'PMN (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : 0,
'UpperLimit' : 0.8,

'1-1_Line' : True,
'0_Line' : False,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr'}

Data['Peak-line-Err-20cm-Marx-1.4GHz'] = {
'Title' : 'Peak Flux - 20cm vs 1.4GHz Marx',

'xCol' : '20cm Peak',
'yCol' : 'Marx 1.4GHz Peak',
'xColErr' : '20cm Peak_Err',

'xLabel' : '20cm (Jy)',
'yLabel' : 'Marx 1.4GHz (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : 0,
'UpperLimit' : 0.3,

'1-1_Line' : True,
'0_Line' : False,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr'}


Data['Peak-line-Err-20cm-Marx-2.4GHz'] = {
'Title' : 'Peak Flux - 20cm vs 2.4GHz Marx',

'xCol' : '20cm Peak',
'yCol' : 'Marx 2.4GHz Peak',
'xColErr' : '20cm Peak_Err',

'xLabel' : '20cm (Jy)',
'yLabel' : 'Marx 2.4GHz (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : 0,
'UpperLimit' : 0.3,

'1-1_Line' : True,
'0_Line' : True,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr'}

Data['Peak-line-Err-36cm-SUMSS'] = {
'Title' : 'Peak Flux - MOST 0.843 GHz vs SUMSS',

'xCol' : '36cm Peak',
'yCol' : 'SUMSS Peak',
'xColErr' : '36cm Peak_Err',
'yColErr' : 'SUMSS Peak_Err',

'xLabel' : 'MOST 0.843 GHz (Jy)',
'yLabel' : 'SUMSS (Jy)',

'UseErrors' : True,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : True,
'LowerLimit' : 0,
'UpperLimit' : 1.5,

'1-1_Line' : True,
'0_Line' : True,

'PointColour' : 'blue',
'PointType' : 'o',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit'    : 'lineErr',
'Exclude' : ['043856-672153','053543-660204','051537-672128','052502-693840']
}

# #===========================================================================================================================================
# #========================================================  Separation  =====================================================================
# #===========================================================================================================================================

Data['Separation-6cm-PMN'] = {
'Title' : 'Separation - ATCA 4.8 GHz vs PMN',

'xCol' : 'Sep True Delta RA 6cm PMN',
'yCol' : 'Sep Delta DEC 6cm PMN',
# 'xColErr' : '36cm Integ_Err',
# 'yColErr' : 'SUMMS Integ_Err',

'xLabel' : 'RA (arcsec)',
'yLabel' : 'DEC (arcsec)',

'UseErrors' : False,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : -35,
'UpperLimit' : 35,

'1-1_Line' : False,
'0_Line' : True,

'PointColour' : 'blue',
'PointType' : '+',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit' : ''}

# Data['Separation-20cm-Marx-1.4GHz'] = {
# 'Title' : 'Separation - 20cm vs 1.4GHz Marx',

# 'xCol' : 'Sep True Delta RA 20cm Marx 1.4',
# 'yCol' : 'Sep Delta DEC 20cm Marx 1.4',

# 'xLabel' : 'RA (arcsec)',
# 'yLabel' : 'DEC (arcsec)',

# 'UseErrors' : False,
# 'xColErrTol' : 0.1,
# 'yColErrTol' : 0.1,

# 'xScale' : '',
# 'yScale' : '',

# 'UseLimitsForAxis' : True,
# 'UseLimitsForFit' : False,
# 'LowerLimit' : -15,
# 'UpperLimit' : 15,

# '1-1_Line' : False,
# '0_Line' : True,

# 'PointColour' : 'blue',
# 'PointType' : '+',
# 'FitColour' : 'red',
# 'Type'   : 'scatter',
# 'Fit' : ''}

Data['Separation-36cm-SUMSS'] = {
'Title' : 'Separation - MOST 0.843 GHz vs SUMSS',

'xCol' : 'Sep True Delta RA 36cm SUMSS',
'yCol' : 'Sep Delta DEC 36cm SUMSS',

'xLabel' : 'RA (arcsec)',
'yLabel' : 'DEC (arcsec)',

'UseErrors' : False,
'xColErrTol' : 0.1,
'yColErrTol' : 0.1,

'xScale' : '',
'yScale' : '',

'UseLimitsForAxis' : True,
'UseLimitsForFit' : False,
'LowerLimit' : -35 ,
'UpperLimit' : 35,

'1-1_Line' : False,
'0_Line' : True,

'PointColour' : 'blue',
'PointType' : '+',
'FitColour' : 'red',
'Type'   : 'scatter',
'Fit' : ''}

proc = subprocess.Popen("wc -l < " + FileName, shell=True, stdout=subprocess.PIPE) # gather the amount of lines in the csv for the progress bar
LinesInFile = proc.communicate()[0]

makeStatsFile()

DictItems = list(Data.keys())
DictItems.sort()

for DictItem in DictItems: # loop through the different plots
	Data['x'] = []
	Data['y'] = []
	Data['xErr'] = []
	Data['yErr'] = []

	ExcludedCount = 0 

	print "============================ Plotting " + DictItem + " ============================"
	pbar = ProgressBar(maxval=int(LinesInFile)).start() # create the progress bar for each plot

	with open(FileName,'rb') as Catalogue: #open the csv
		reader = csv.reader(Catalogue)
		Headings = next(reader) #gather the headings

		pbar.update(1) #move the bar 1 place because of the headings. 

		#get the two important columns 
		GoodSource = Headings.index('Good Source')
		Blended = Headings.index('Blended')
	 
	 	for RowNum, Row in enumerate(reader):  #loop through the csv
	 		pbar.update(RowNum + 1)

	 		Exclude = False
	 		xGoodError = False
	 		yGoodError = False

	 		if 'Exclude' in Data[DictItem]: #check to see if the current line is to be ignored. 
	 			for ExcludeSource in Data[DictItem]['Exclude']:
	 				if ExcludeSource == Row[Headings.index('Name')]:
	 					Exclude = True

	 		for Conditions in list(ConditionalColumns):
	 			if Row[Headings.index(Conditions)] != ConditionalColumns[Conditions]:
	 				Exclude = True
	 				ExcludedCount += 1

	 		if Exclude == False: #check to see if this source is usable. 		
	 			#============== TO DO: add in a way to have a fixed column of data or possible have a way to bin it on conditions defined above? ========
	 			if bool(Row[Headings.index(Data[DictItem]['xCol'])]) and bool(Row[Headings.index(Data[DictItem]['yCol'])]):
	 				#if 1 == 1:
	 				if (Data[DictItem]['UseLimitsForFit'] == False) or (float(Row[Headings.index(Data[DictItem]['xCol'])]) < Data[DictItem]['UpperLimit'] and float(Row[Headings.index(Data[DictItem]['yCol'])]) < Data[DictItem]['UpperLimit']):
		 				#if the data exisits, then use it
		 				Data['x'].append(float(Row[Headings.index(Data[DictItem]['xCol'])]))
		 				Data['y'].append(float(Row[Headings.index(Data[DictItem]['yCol'])]))
		 				
		 				if 'xColErr' in Data[DictItem]: #if an error column is specified then use it. 
		 					if bool(Row[Headings.index(Data[DictItem]['xColErr'])]):
		 						#if we have a small error (i.e. smaller than a preset value in the plot def), then use the error val
		 						if (float(Row[Headings.index(Data[DictItem]['xCol'])]) * float(Data[DictItem]['xColErrTol'])) > float(Row[Headings.index(Data[DictItem]['xColErr'])]) and float(Row[Headings.index(Data[DictItem]['xColErr'])]) > 0:
		 							xGoodError = True

		 				if 'yColErr' in Data[DictItem]: #if an error column is specified then use it. 
		 					if bool(Row[Headings.index(Data[DictItem]['yColErr'])]):
		 						#if we have a small error (i.e. smaller than a preset value in the plot def), then use the error val
		 						if ( float(Row[Headings.index(Data[DictItem]['yCol'])]) * float(Data[DictItem]['yColErrTol'])) > float(Row[Headings.index(Data[DictItem]['yColErr'])]) and float(Row[Headings.index(Data[DictItem]['yColErr'])]) > 0:
		 							yGoodError = True

		 				if xGoodError == True: #get the error val
		 					Data['xErr'].append(float(Row[Headings.index(Data[DictItem]['xColErr'])]))
		 				elif 'xColErrTol' in Data[DictItem]: # assume an error value that is set in the plot definition
		 					Data['xErr'].append(float(Row[Headings.index(Data[DictItem]['xCol'])]) * float(Data[DictItem]['xColErrTol']))

	 					if yGoodError == True: #get the error val
		 					Data['yErr'].append(float(Row[Headings.index(Data[DictItem]['yColErr'])]))
		 				elif 'yColErrTol' in Data[DictItem]: # assume an error value that is set in the plot definition
		 					Data['yErr'].append(float(Row[Headings.index(Data[DictItem]['yCol'])]) * float(Data[DictItem]['yColErrTol']))

	if bool(Data['x']) and bool(Data['y']): #generate the plot if there is some data
		fig = plt.figure()
		fig.set_size_inches(5.7 , 5.5)
		
		ax = fig.add_subplot(111)
		ax.set_aspect('equal')

		if bool(Data[DictItem]['xScale']): #set x the plot scale
			ax.set_xscale(Data[DictItem]['xScale'])
		
		if bool(Data[DictItem]['yScale']): #set y the plot scale
			ax.set_yscale(Data[DictItem]['yScale'])

		if Data[DictItem]['UseLimitsForAxis'] == True: #set the plot limit and the limits for drawing the fit
			plt.xlim(Data[DictItem]['LowerLimit'], Data[DictItem]['UpperLimit'])
			plt.ylim(Data[DictItem]['LowerLimit'], Data[DictItem]['UpperLimit'])

			dx = np.arange(float(Data[DictItem]['LowerLimit']), float(Data[DictItem]['UpperLimit'] + 0.1), 0.1)
		else:
			dx = np.arange(0,5,0.1)

		#create the labels
		if 'xLabel' in Data[DictItem]:
			if Data[DictItem]['xLabel'] != '':
				plt.xlabel(Data[DictItem]['xLabel'])
		else:
			plt.xlabel(Data[DictItem]['xCol'])

		if 'yLabel' in Data[DictItem]:
			if Data[DictItem]['yLabel'] != '':
				plt.ylabel(Data[DictItem]['yLabel'])
		else:
			plt.ylabel(Data[DictItem]['yCol'])

		if 'Title' in Data[DictItem]:
			if Data[DictItem]['Title'] != '':
				plt.title(Data[DictItem]['Title'])

 
		#Plot the data points 
		#========= TO DO: add in an if statement to handle other types of plots ===================
		l = ax.scatter(Data['x'], Data['y'],  s=10, c = Data[DictItem]['PointColour'], marker=str(Data[DictItem]['PointType']), zorder=2)

		if bool(Data['xErr']) and bool(Data['yErr']) and Data[DictItem]['UseErrors'] == True:
			plt.errorbar(Data['x'], Data['y'], xerr = Data['xErr'], yerr = Data['yErr'], ls = 'none')

		#do the fitting and plot it
		if Data[DictItem]['Fit'] == "alpha":			
			if bool(Data[DictItem]['yColErr']): 
				fitting_obj = kmpfit.simplefit(alpha, [ 1, 1 ], Data['x'], Data['y'], err = Data[DictItem]['yColErr']) 
			else:
				fitting_obj = kmpfit.simplefit(alpha, [ 1, 1 ], Data['x'], Data['y'])

			plt.plot(dx, alpha(fitting_obj.params, dx), c = Data[DictItem]['FitColour'])

		elif Data[DictItem]['Fit'] == "line":
			fitting_obj = kmpfit.simplefit(line, [ 1, 1 ], Data['x'], Data['y'])
			plt.plot(dx, line(fitting_obj.params, dx), c = Data[DictItem]['FitColour'])		

		elif Data[DictItem]['Fit'] == "lineErr":
			linear = odr.Model( lineErr )
			mydata = odr.Data( np.asarray(Data['x']), np.asarray(Data['y']), wd = np.asarray(Data['xErr']), we = np.asarray(Data['yErr']) )
			myodr = odr.ODR( mydata, linear, beta0 = [ 1, 1 ] )
			myoutput = myodr.run()

			plt.plot(dx, lineErr(myoutput.beta, dx), c = Data[DictItem]['FitColour'])

			print myoutput.beta
			print myoutput.pprint()


		if '1-1_Line' in Data[DictItem]: #create a 1-1 dotted black line. 
			if Data[DictItem]['1-1_Line'] == True:
				plt.plot(dx, line([1,0], dx), '--', c = 'black', linewidth=1)

		if '0_Line' in Data[DictItem]: #create a 1-1 dotted black line. 
			if Data[DictItem]['0_Line'] == True:
				plt.plot((float(Data[DictItem]['LowerLimit']), float(Data[DictItem]['UpperLimit'])), (0,0), '--', c = 'black', linewidth=1, zorder=0)
				plt.plot((0,0), (float(Data[DictItem]['LowerLimit']), float(Data[DictItem]['UpperLimit'])), '--', c = 'black', linewidth=1, zorder=0)

		
		pbar.finish()
		#remove the old plot and save the new one
		os.system('rm Plots/LMC42/' + DictItem + '.eps')
		plt.savefig('Plots/LMC42/' + DictItem + '.eps')
		plt.close()
		
		#get the relevant info and send it to the stats file. 
		if Data[DictItem]['Fit'] == "alpha" or Data[DictItem]['Fit'] == "line":	
			makeStatsFile(DictItem,Data,fitting_obj)
		elif Data[DictItem]['Fit'] == "lineErr":
			makeStatsFile(DictItem,Data,myoutput)
		else:
			makeStatsFile(DictItem,Data)
	else:
		pbar.finish()
		print "ERROR: NO DATA TO SEND TO THE PLOT"

	# os.system('rm Plots/test.txt')
	# #print out the fitting info
	# with open("Plots/test.txt", "a") as myfile:
	# 	myfile.write("============================ Logging " + DictItem + " ============================\n")	

	# 	if Data[DictItem]['Fit'] == "alpha" or Data[DictItem]['Fit'] == "line":	
	# 		Text = obj_print(fitting_obj, DictItem)
	# 		myfile.write(Text)
	# 	elif Data[DictItem]['Fit'] == "lineErr":
	# 		myfile.write(str(myoutput.pprint()))
	

	

	print "Sources Included: " + str(len(Data['x']))

print "Sources Excluded: " + str(ExcludedCount)	


# Calling the function
notify(title    = 'Plot Comparrions Completed',
       subtitle = 'Plot-CatalogueComparrions.py',

       message  = 'Started: ' + str(startTime)[11:19] + ' - Time Taken: ' + str(datetime.now()-startTime)[:7] )























