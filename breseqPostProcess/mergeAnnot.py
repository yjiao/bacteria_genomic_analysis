# Merge annotations with their evidence, and output to a new folder
# this will then go into generateMatrices.py

from os import listdir, makedirs
import heapq as h
import fnmatch
from os.path import isfile, join, exists
from collections import defaultdict
from sys import argv


# functions
def skipComments(fh):
# skip all the lines in the beginning
# "rewind" the file pointer
    comments = []
    line = fh.readline()
    firstchar = line[0]
    dist=len(line)
    while firstchar == '#':
	comments.append(line)
	line = fh.readline()
	dist = len(line)
	if line=='':
	    return [fh, True, []] #file handle, is it done
	else:
	    firstchar = line[0]
    fh.seek(-dist,1)
    return [fh, False, comments]

def getAnnotations(fh):
    annotations = []
    line = fh.readline()
    dist = len(line)
    line = line.split('\t')
    while line[0] not in evidence_types:
	annotations.append('\t'.join(line))
	line = fh.readline()
	if line == '':
	    print '------------------------------------------------------!!!!!!!!!!!annotations empty, manually check file'
	    return fh, annotations
	dist = len(line)
	line = line.split('\t')
    fh.seek(-dist,1)
    return fh, annotations

def createEvDict(evidence):
    evdict = defaultdict(list)
    for line in evidence:
	evID = line[1]
	evdict[evID] = line
    return evdict

def merge(annot, evDict):
    for line in annot:
	evlinks = line[2].split(',')
	evlines = []
	for link in evlinks:
	    entry = evDict[link]

	    atype = line[0]
	    etype = entry[0]

	    entry = filter(lambda x: x.find('html_gene_name')==-1, entry)
	    entry = ';'.join(entry)
	    evlines.append(entry)
	    
#	    if atype == 'SUB' and etype == 'JC':
#		print line
#		print entry
	evfield = 'evidence=' + '@@@'.join(evlines)
	line.append(evfield)
    # note there's no need to put "line" back into annot, since line itself is a list and thus is a soft copy
    return annot

def writeFile(newname, data):
    with open(newname, 'w') as fh:
	fh.write('\n'.join(['\t'.join(line) for line in data]))

################################################################################################################################################
# get annotated files--snps, dels, etc only
base = argv[1]
basepath = base + '/postProcess/'
path = basepath+'01annotated/'
annotfiles = [f for f in listdir(path) if isfile(join(path, f))]
strains = [f.replace('.gd', '') for f in annotfiles]
evidence_types = set(['RA', 'MC', 'JC', 'UN'])

for strain, file in zip(strains, annotfiles):
    print strain
    with open(path + file) as fh:
	[fh, isclosed, comments ]= skipComments(fh)
	[fh, annotations]= getAnnotations(fh)

	if len(annotations) == 0:
	    # if there are no annotations, then there were only marginal calls. There's nothing to merge, so skip.
	    # need to write out an empty file?
	    writeFile(newpath+file, comments)
	    continue

	# format annotations
	annotations = [line.strip().split('\t') for line in annotations]

	evidence = fh.readlines()
	# format evidence
	evidence = [line.strip().split('\t') for line in evidence]
	
	# double check line parsing
	check1 = int(annotations[len(annotations)-1][1])
	check2 = int(evidence[0][1])
	
	if check2 - check1 !=1:
	    print '------------------------------------------------------!!!!!!!!!!!LINE READ ERROR, manually check file in 01annotated'
	    print 'This may require copying evidence from 00 to 01 and 02'
	    print 'Remember to also change the second field to increasing indices (e.g. 1, 1 to 1,2)'
	
	# create dictionary of evidence lines based on evidence ID (for recall later by annotation lines)
	evDict = createEvDict(evidence)
	
	# merge evidence into the annotation line as an extra field
	# format: @@@ separates evidences, ; replaces tabs in original
	annotations = merge(annotations, evDict)

	# output file
	if not exists(basepath + '05annotEvMerged'):
	    makedirs(basepath + '05annotEvMerged')
	newpath = basepath + '05annotEvMerged/'
	writeFile(newpath+file, annotations)

	
	