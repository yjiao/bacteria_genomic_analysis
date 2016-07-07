# just output max score and min score for now, too confusing otherwise!!
# also time crunch = no time to think = crappy algos

from os import listdir, makedirs
import heapq as h
import fnmatch
from os.path import isfile, join, exists
from collections import defaultdict
from sys import argv

def splitfile(fh, path, filename):
    name = filename.split('.')[0]
    lines = fh.readlines()
    
    totmaxline = []
    totminline = []
    topmaxline = []
    topminline = []
    botmaxline = []
    botminline = []

    for line in lines:
	line = line.replace('\n', '') # note can't use strip here, it will kill all the trailing tabs
	line = line.split('\t')
	
	# top and bottom strands
	totmax = []
	totmin = []
	topmax = []
	topmin = []
	botmax = []
	botmin = []

	for field in line:
	    if field=='':
		topmax.append('NA')
                topmin.append('NA')
                botmax.append('NA')
		botmin.append('NA')
		totmax.append('NA')
                totmin.append('NA')
		continue

	    field = field.split(',')
	    tops = []
	    bots = []
	    tots = []
	    for subfield in field:
		if subfield == '':
		    top = 'NA'
		    bot = 'NA'
		    flag=True
		else:
		    [top, bot] = subfield.split('/')
		    flag=False
		    tops.append(top)
		    bots.append(bot)
		    tots.append(int(top)+int(bot))
	    
	    if not flag:
		tops = map(int, tops)
		bots = map(int, bots)

		totmax.append(max(tots))
		totmin.append(min(tots))
		topmax.append(max(tops))
		topmin.append(min(tops))
		botmax.append(max(bots))
		botmin.append(min(bots))
	    else:
		totmax.append('NA')
		totmin.append('NA')
		topmax.append('NA')
		topmin.append('NA')
		botmax.append('NA')
		botmin.append('NA')
    
	
	totmaxline.append('\t'.join(map(str,totmax)))
	totminline.append('\t'.join(map(str,totmin)))
	topmaxline.append('\t'.join(map(str,topmax)))
	topminline.append('\t'.join(map(str,topmin)))
	botmaxline.append('\t'.join(map(str,botmax)))
	botminline.append('\t'.join(map(str,botmin)))
    
    with open(path+name+'_total_max.txt', 'w') as f:
	f.write('\n'.join(totmaxline))
    with open(path+name+'_total_min.txt', 'w') as f:
	f.write('\n'.join(totminline))
    
    with open(path+name+'_top_max.txt', 'w') as f:
	f.write('\n'.join(topmaxline))
    with open(path+name+'_top_min.txt', 'w') as f:
	f.write('\n'.join(topminline))
    
    with open(path+name+'_bottom_max.txt', 'w') as f:
	f.write('\n'.join(botmaxline))
    with open(path+name+'_bottom_min.txt', 'w') as f:
	f.write('\n'.join(botminline))


# find the right files to use
#basepath = 'fullRun20151221/postProcess/04merged/'
basepath = argv[1]
files = [f for f in listdir(basepath) if isfile(join(basepath, f))]
filesOfInterest = set(['total_cov.txt', 'major_cov.txt', 'minor_cov.txt'])
# note total_cov.txt, minor_cov.txt, major_cov.txt is too difficult to parse for now, also time crunch = crappy code
files = filter(lambda x: x in filesOfInterest, files)
files = filter(lambda x: x.find('.txt') != -1, files)

if not exists(basepath + 'minmax_params'):
    makedirs(basepath + 'minmax_params')

for file in files:
    with open(basepath + file) as fh:
	print basepath+file
	splitfile(fh, basepath+'minmax_params/', file)


