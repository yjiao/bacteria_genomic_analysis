# just output max score and min score for now, too confusing otherwise!!
# also time crunch = no time to think = crappy algos

from os import listdir, makedirs
import heapq as h
import fnmatch
from os.path import isfile, join, exists
from collections import defaultdict
import sys

def splitfile(fh, path, filename):
    name = filename.split('.')[0]
    lines = fh.readlines()
    
    maxfile = []
    minfile = []
    for line in lines:
	line = line.replace('\n', '') # note can't use strip here, it will kill all the trailing tabs
	line = line.split('\t')
	
	maxline = []
	minline = []

	for field in line:
	    field = field.split(',')
	    field = [en if len(en)>0 else 'NA' for en in field]
	    #try:
	    if len(field[0]) > 0:
		#print 'processed:', field, len(field)
		fieldmax = [float('-Inf') if f=='NA' else f for f in field]
		fieldmax = [float('-Inf') if f=='NT' else f for f in fieldmax]
		maxfield = max([float(f) for f in fieldmax])
		fieldmin = [float('Inf') if f=='NA' else f for f in field]
		fieldmin = [float('Inf') if f=='NT' else f for f in fieldmin]
		minfield = min([float(f) for f in fieldmin])

		if maxfield == float('-Inf'):
		    maxfield = 'NA'
		if minfield == float('Inf'):
		    minfield = 'NA'

	    else:
		#print '-----skipped:', field, len(field)
		maxfield = 'NA'
		minfield = 'NA'
	    #except: print 'ERROR:', field

	    maxline.append(maxfield)
	    minline.append(minfield)
	
	maxline = '\t'.join(map(str, maxline))
	minline = '\t'.join(map(str, minline))

	maxfile.append(maxline)
	minfile.append(minline)
    with open(path+name+'_max.txt', 'w') as f:
	f.write('\n'.join(maxfile))
    with open(path+name+'_min.txt', 'w') as f:
	f.write('\n'.join(minfile))


# find the right files to use
#basepath = 'fullRun20151221/postProcess/04merged/'
basepath = sys.argv[1]# +'/postProcess/04merged/'
files = [f for f in listdir(basepath) if isfile(join(basepath, f))]
notparams = set(['mutations.txt', 'matrix.txt', 'strains.txt', 'indexhtml.txt', 'side_2_annotate_key.txt', 'total_cov.txt', 'minor_cov.txt', 'major_cov.txt', 'side_1_annotate_key.txt', 'prediction.txt', 'evidencehtml.txt', 'indexhtml.txt', 'key.txt', 'ref_base.txt'])
# note total_cov.txt, minor_cov.txt, major_cov.txt is too difficult to parse for now, also time crunch = crappy code
files = filter(lambda x: x not in notparams, files)
files = filter(lambda x: x.find('.txt') != -1, files)

if not exists(basepath + 'minmax_params'):
    makedirs(basepath + 'minmax_params')

for file in files:
    with open(basepath + file) as fh:
	print basepath+file
	splitfile(fh, basepath+'minmax_params/', file)


