# i/o operations will be a lot easier if read evidence and mutation calls were split into different files
# 2016 01 25: also split into annotated SNP/DEL/SUB etc only and evidence only

from os import listdir, makedirs
from os.path import isfile, join, exists
from sys import argv

def skiptoevidence(fh):
    annotations = []
    line = fh.readline()
    dist = len(line)
    line = line.split('\t')
    while line[0] not in evidence_types:
	annotations.append('\t'.join(line))
	line = fh.readline()
	if line == '':
	    print 'manually check this case!! Defect in breseq gd file processing may have caused evidence to be missing.'
	    return fh, annotations
	dist = len(line)
	line = line.split('\t')
    fh.seek(-dist,1)
    return fh, annotations

def skipComments(fh):
# skip all the lines in the beginning
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

def writefile(name, content):
    with open(name,'w') as f:
	f.write(''.join(content)) # newline char are still there

base = argv[1]
basepath = base + '/postProcess/'
path = basepath+'01annotated/'
files = [f for f in listdir(path) if isfile(join(path, f))]
strains = [f.replace('.gd', '') for f in files]
evidence_types = set(['RA', 'MC', 'JC', 'UN'])

fharray = [open(join(path,f)) for f in files]

if not exists(basepath+'02evidence'):
    makedirs(basepath + '02evidence')

if not exists(basepath+'03annotOnly'):
    makedirs(basepath + '03annotOnly')

boolarray = []
for i, fh in enumerate(fharray):
    evcomments = ['#Evidence file split from annotations in 01annotated\n']
    ancomments = ['#Annotations file split from annotations in 01annotated\n']

    strain = strains[i]
    print strain
    evName = basepath + '02evidence/' + strain + '_evidence.gd'
    annotName = basepath + '03annotOnly/' + strain + '_annotation.gd'
    fh, fh_done, comments = skipComments(fh)

    if fh_done:
	print '...file empty'
	writefile(evName, comments)
	writefile(annotName, comments)
    else:
	fh, annotations = skiptoevidence(fh)
	evidence_lines = fh.readlines()

	evcomments.extend(evidence_lines)
	ancomments.extend(annotations)
	writefile(evName, evcomments)
	writefile(annotName, ancomments)

    boolarray.append(files[i].replace('.gd', '')!=strain)
    fh.close()

print 'Any mismatches between strain names and gd files:', any(boolarray)
