# NOTES
# DONE: make sure evidence is actually appended in each case
# DONE: implement INS
# DONE:implement large DEL
# DONE:implement large INS
# DONE figure out how to handle multi-base stuff such as substitutions
# run this AFTER generateMatrices.py, requires mutations.txt as well as matrix.txt
# plan: go back to 08_mutation_identificationfolder and grab the marginal evidence to see if they match any of our calls

from os import listdir, makedirs
from sys import argv
import heapq as h
import fnmatch
from os.path import isfile, join, exists
from collections import defaultdict

verbose = False
# define a print functon only when verbose = True
if verbose:
    def verboseprint(*args):
	print ' '.join(map(str, args))
else:
    verboseprint = lambda *a: None

# classes
class mutation():
    def __init__(self, allkeys, line):
	line = line.replace('\n', '')
	line = line.split('\t')
	for i, key in enumerate(allkeys):
	    self.__dict__[key] = line[i]
	self.position = int(self.position)
	self.size = int(self.size)
	if len(self.repeat_ref_copies) > 0:
	    self.repeat_ref_copies = int(self.repeat_ref_copies)

	self.evidence = []

    def __cmp__(self, other):
	return cmp(self.position, other.position)

class entry:
    def __init__(self, evidence_url, index_url, name, evidence, allkeys):
	for key in allkeys:
	    self.__dict__[key] = []

	self.evidence = evidence_url # this is a list
	self.index = index_url
	self.name = name
	self.evidenceType = []
	self.evidenceCount = len(evidence)

	# add evidence annotations
	for line in evidence:
	    for i, ev in enumerate(line):
		ev = ev.split('=')
		if i == 0:
		    self.evidenceType.append(ev[0])
		if len(ev) < 2:
		    continue
		key = ev[0]
		val = ev[1]
		if key not in allkeys:
		    continue
		self.__dict__[key].append(val)

class evidence:
    def __init__(self, line, allkeys, strain):
	self.strain = strain

	for key in allkeys:
	    self.__dict__[key] = ''
	
	self.evidenceType = line[0]
	if self.evidenceType == 'RA':
	    self.position = int(line[4])
	    self.alt = line[6]
	    self.alt_minor = line[7]
	
	elif self.evidenceType == 'MC':
	    self.start = int(line[4])
	    self.end = int(line[5])
	    self.startRange = int(line[6])
	    self.endRange = int(line[7])
	    self.position = self.start
	
	elif self.evidenceType == 'UN':
	    self.start = int(line[4])
	    self.end = int(line[5])
	    self.position = self.start

	for i, ev in enumerate(line):
	    ev = ev.split('=')
	    if len(ev) < 2:
		continue
	    key = ev[0]
	    val = ev[1]
	    if key not in allkeys:
		continue
	    self.__dict__[key] = val
	
	if self.evidenceType == 'JC':
	    side1pos = int(self.key.split('__')[1])
	    side2pos = int(self.key.split('__')[4])
	    self.side1pos = min(side1pos, side2pos)
	    self.side2pos = max(side1pos, side2pos)
	    self.sidepos = [side1pos, side2pos]
	    self.position = self.side1pos

    def __cmp__(self, other):
	return cmp(self.position, other.position)

# functions
# file I/O
def skiptoevidence(fh):
    annotations = []
    line = fh.readline()
    dist = len(line)
    line = line.split('\t')
    while line[0] not in evidence_types:
	annotations.append('\t'.join(line))
	line = fh.readline()
	if line == '':
	    print 'check this case!!'
	    return [fh, True]
	dist = len(line)
	line = line.split('\t')
    fh.seek(-dist,1)
    return [fh, False]

def parseMatrix(file):
    matrix = []
    for line in file.readlines():
	line = line.strip().split('\t')
	line = map(int, line)
	matrix.append(line)
    return matrix

def skipComments(fh):
# skip all the lines in the beginning
# "rewind" the file pointer
    line = fh.readline()
    firstchar = line[0]
    dist=len(line)
    while firstchar == '#':
	line = fh.readline()
	dist = len(line)
	if line=='':
	    return [fh, True, []] #file handle, is it done
	else:
	    firstchar = line[0]
    fh.seek(-dist,1)
    return [fh, False]

# output 
def write_matrix(outname, mut, strains):
    with open(outname, 'a') as f:
	newline = [0 for s in strains]
	for evidence in mut.evidence:
	    newline[strains.index(evidence.strain)] = 1
	
	f.write('\t'.join(map(str, newline))+'\n')

def write_mutations(outname, mut, allkeys):
    with open(outname, 'a') as f:
	newline = []
	for key in allkeys:
	    newline.append(str(mut.__dict__[key]))
	newline = '\t'.join( newline ) + '\n'
	f.write(newline)

def write_url(outname, mut, strains):
    with open(outname, 'a') as f:
	newline = ['' for s in strains]
	for strain in mut.strains:
	    evs = ','.join(strain.evidence)
	    newline[strains.index(strain.name)] = evs
	f.write('\t'.join(newline)+'\n')

def write_evidence_params(mut, strains, allkeys, path):

    for param in allkeys:
	outname = path + param + '.txt'
	with open(outname, 'a') as f:
	    newline = [[] for s in strains]
	    for evidence in mut.evidence:
		field = evidence.__dict__[param] # we actually only record one field?
		#field = ','.join(evidence.__dict__[param])
		newline[strains.index(evidence.strain)].append(field)
	    
	    newline2 = ['' for s in strains]
	    for strain, field in zip(strains, newline):
		newline2[strains.index(strain)] = ','.join(map(str, field))
	    f.write('\t'.join(newline2)+'\n')

def write_index(outname, mut, strains):
    with open(outname, 'a') as f:
	newline = ['' for s in strains]
	for strain in mut.strains:
	    newline[strains.index(strain.name)] = strain.index
	
	f.write('\t'.join(map(str, newline))+'\n')

# force calling functions
def forcecall(mutation, fharray, allevkeys):
    verboseprint('---------------------------------------- mutation position:', mutation.position)
    verboseprint('mutation pos:', mutation.position, 'mutation size:',  mutation.size, mutation.ref, mutation.alt)
    # NOTE: can set this lineto be equal to 1 as a positive check--the calls should all be found
    idx = [i for i in range(len(fharray)) if row[i] == 0]
    total = len(idx)
    nfound = 0
    for i in idx:
	strain = strains[i]
	fh = fharray[i]
	if fh.closed:
	    continue
	found = False

	##################### SNP #####################
	if mutation.type == 'SNP':
	    method = 'SNP: findRA'
	    found = findRA(mutation, fh, strain)

	##################### INS #####################
	elif mutation.type == 'INS':
	    if isinstance(mutation.repeat_ref_copies, basestring):
		# ostensibly not a repeat (except when breseq messes up)
		method='INS: findIndelInterval'
		found = findIndelInterval(mutation, fh, strain)
		if not found:
		    method='INS: findSUB2'
		    found = findSUB2(mutation, fh, strain)

	    elif isinstance(mutation.repeat_ref_copies, int):
		# this is a repeat
		method='INS: findRepeat'
		found = findRepeat(mutation, fh, strain)

	    if not found:
		# look for junctions
		method = 'INS: findJunction'
		found = findJunction(mutation, fh, strain)
	    else:
		method = 'INS: none'

	##################### SUB #####################
	elif mutation.type == 'SUB':
	    if len(mutation.alt) > len(mutation.ref):
		method='SUB: findSUB2'
		found = findSUB2(mutation, fh, strain)
	    else:
		method='SUB: findSUB'
		found = findSUB(mutation, fh, strain)

	##################### DEL #####################
	elif mutation.type=='DEL':
	    found = False
	    if isinstance(mutation.repeat_ref_copies, basestring):
		# ostensibly not a repeat (except when breseq messes up)
		method='DEL: findIndelInterval'
		found = findIndelInterval(mutation, fh, strain)
		if not found:
		    method='DEL: findSUB2'
		    found = findSUB2(mutation, fh, strain)

	    elif isinstance(mutation.repeat_ref_copies, int):
		# this is a repeat
		method='DEL: findRepeat'
		found = findRepeat(mutation, fh, strain)

	    if not found:
		method='DEL: findJunction'
		found = findJunction(mutation, fh, strain)

	    if not found:
		method = 'DEL: findMC'
		found = findMC(mutation, fh, strain)

	##################### OTHER?? #####################
	else:
	    method='??: none'
	    #print mutation.type, mutation.size, mutation.position
	    continue

	if found:
	    nfound += 1
	else:
	    pass
	    #print method, strain, 'not found'

    print nfound, ' / ', total, mutation.type
#    if len(mutation.evidence)!= total:
#	print mutation.__dict__
    verboseprint(nfound, 'found of', total)
    write_matrix(outpath+'/matrix.txt', mutation, strains)
    write_evidence_params(mutation, strains, allevkeys, outpath)

def findMC(mutation, fh, strain):
    # start positon should be the same
    # ev.end = start.size?
    fhstart = fh.tell() # get starting fh position if we need to rewind later
    line = fh.readline()
    ev = evidence(line.strip().split('\t'), allevkeys, strain)
    if isinstance(mutation.repeat_ref_copies, basestring):
	# this is not a repeat
	min_position = mutation.position - mutation.size
	copies = mutation.size

    else: 
	# this is a repeat
	copies = max(int(mutation.repeat_ref_copies)*len(mutation.repeat_seq), int(mutation.repeat_new_copies)*len(mutation.repeat_seq))
	if mutation.type=='DEL':
	    min_position = mutation.position - copies
	elif mutation.type=='INS':
	    min_position = mutation.position - len(mutation.repeat_seq)*int(mutation.repeat_ref_copies) - len(mutation.repeat_seq)*int(mutation.repeat_new_copies) + 1
    
    if ev.position > mutation.position + copies:
	fh.seek(fhstart)
	return False

    # get rid of useless evidence lines: RA, MC
    while ev.evidenceType == 'RA':
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    return False
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
   
    evs = [ev]
    while ev.position <= mutation.position + copies:
    #for i in range(copies):
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    break
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
	evs.append(ev)
    MCs = [e for e in evs if e.evidenceType=='MC']

    test_positions = [e.position for e in MCs]
    test_positions_pass = [1 if min_position <= t <= min_position + copies else 0 for t in test_positions]
    
    # assume that junction evidence in permitted range is used to support the insertion...
    if sum(test_positions_pass) > 0:
	fh.seek(fhstart)
	mutation.evidence.append(ev)
	return True
    else:
	fh.seek(fhstart)
	return False

def findJunction(mutation, fh, strain):
    # when finding RAs fail, then there might be junction evidence. Unfortuantely these have weird annotations...
    fhstart = fh.tell() # get starting fh position if we need to rewind later
    line = fh.readline()
    ev = evidence(line.strip().split('\t'), allevkeys, strain)
    if isinstance(mutation.repeat_ref_copies, basestring):
	# this is not a repeat
	min_position = mutation.position - mutation.size
	copies = mutation.size

    else: 
	# this is a repeat
	copies = max(int(mutation.repeat_ref_copies)*len(mutation.repeat_seq), int(mutation.repeat_new_copies)*len(mutation.repeat_seq))
	if mutation.type=='DEL':
	    min_position = mutation.position - copies
	elif mutation.type=='INS':
	    min_position = mutation.position - len(mutation.repeat_seq)*int(mutation.repeat_ref_copies) - len(mutation.repeat_seq)*int(mutation.repeat_new_copies) + 1
    
    # get rid of useless evidence lines: RA, MC
    while ev.evidenceType == 'RA':
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    return False
	ev = evidence(line.strip().split('\t'), allevkeys, strain)

    while ev.evidenceType == 'MC':
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    return False
	ev = evidence(line.strip().split('\t'), allevkeys, strain)

    if ev.evidenceType != 'JC':
	fh.seek(fhstart)
	return False
    
    if ev.position > mutation.position + copies:
	fh.seek(fhstart)
	return False

    while ev.side2pos < min_position:
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    return False
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
	if ev.evidenceType != 'JC':
	    break
    
    if ev.evidenceType != 'JC':
	fh.seek(fhstart)
	return False
    
    evs = [ev]
    while ev.position <= mutation.position + copies:
    #for i in range(copies):
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    break
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
	if ev.evidenceType != 'JC':
	    fh.seek(fhstart)
	    break
	evs.append(ev)
    JCs = [e for e in evs if e.evidenceType=='JC']
    
    test_positions_pass = []
    for j in JCs:
	if j.side1pos - int(j.coverage_minus) <= mutation.position <= (j.side2pos + int(j.coverage_plus)): #abs(int(j.alignment_overlap))):
	    test_positions_pass.append(1)
	else:
	    test_positions_pass.append(0)
#    test_positions = [e.position for e in JCs]
#    test_positions.extend([t  + int(e.alignment_overlap) for t in test_positions])
#    test_positions_pass = [1 if min_position <= t <= min_position + copies else 0 for t in test_positions]
    
    # assume that junction evidence in permitted range is used to support the insertion...
    if sum(test_positions_pass) > 0:
	fh.seek(fhstart)
	mutation.evidence.append(ev)
	return True
    else:
	fh.seek(fhstart)
	return False

def findSUB2(mutation, fh, strain):
    # when alt is bigger than reference, sometimes it's shifted to the wrong direction
    #	    0123456
    # ref: |...A...| len(ref) = 1, len(alt) = 4, size=3
    # alt: |TTTT...|
    #      |.TTTT..|
    #      |..TTTT.|
    #      |...TTTT|
    # criteria:
    # 1. continuous
    # 2. match alt sequence exactly
    # 3. fall into this range
    # don't require the ref base to be = '.' for now, since we don't know how it will align

    min_pos = mutation.position - len(mutation.alt) + 1
    max_pos = mutation.position + len(mutation.alt) - len(mutation.ref)

    fhstart = fh.tell() # get starting fh position if we need to rewind later
    line = fh.readline()
    dist = len(line)
    ev = evidence(line.strip().split('\t'), allevkeys, strain)
    
    if ev.position > max_pos:
	fh.seek(fhstart)
	return False

    while ev.position < min_pos:
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    return False
	dist = len(line)
	ev = evidence(line.strip().split('\t'), allevkeys, strain)

    evs = [ev]
    while ev.position <= max_pos:
    #for i in range(max_pos - min_pos + 1):
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    break
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
	evs.append(ev)
    
    RAs = [e for e in evs if e.evidenceType=='RA']

    test_positions = [e.position for e in RAs]
    test_positions_pass = [1 if min_pos <= t <= max_pos else 0 for t in test_positions]

    evs_pass = [e for i,e in enumerate(RAs) if test_positions_pass[i]==1]
    test_alt = ''.join([e.alt for e in evs_pass])
    test_alt_minor = ''.join([e.alt_minor for e in evs_pass])
    
    test_positions = [t for t in test_positions if min_pos <= t <= max_pos]
    # occasionally there is only 1 position for multiple insertions, and the calls may not be in order. In this case, we need to see if all the alts at one location are permuations of mutation.alt
    if len(set(test_positions)) == 1: # the RAs might be mixed up here, check for permutation
	if sorted(test_alt) == sorted(mutation.alt) or sorted(test_alt_minor) == sorted(mutation.alt):
	    fh.seek(fhstart)
	    mutation.evidence.append(ev)
	    return True

    # since this case is miscalled, we're not guaranteed exactly 1 deletion in this stretch of repeats. look for at least 1 '.' instead.
    if mutation.alt in test_alt or mutation.alt in test_alt_minor:
	fh.seek(fhstart)
	mutation.evidence.append(ev)
	return True
    else:
	fh.seek(fhstart)
	return False

def findIndelInterval(mutation, fh, strain):
    # when repeat is not explicitly stated, but breseq might have missed this case, so we want to make sure we also look for repeats
    # assume mutation.size = 1!!!
    refbase = genome[mutation.position-1] # this is where the difference is, supposedly
    start = mutation.position-1
    
    left_shift = 0
    while genome[start - left_shift] == refbase:
	left_shift += 1
    left = start - left_shift + 1

    right_shift = 0
    while genome[start + right_shift] == refbase:
	right_shift += 1
    right = start + right_shift
    
    # genome is 1-indexed
    min_pos = left
    max_pos = right + 1

    fhstart = fh.tell() # get starting fh position if we need to rewind later
    line = fh.readline()
    dist = len(line)
    ev = evidence(line.strip().split('\t'), allevkeys, strain)
    
    if ev.position > max_pos:
	fh.seek(fhstart)
	return False

    while ev.position < min_pos:
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    return False
	dist = len(line)
	ev = evidence(line.strip().split('\t'), allevkeys, strain)

    evs = [ev]
    while ev.position <= right:
    #for i in range(right - left - 1):
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    break
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
	evs.append(ev)
    
    RAs = [e for e in evs if e.evidenceType=='RA']

    test_positions = [e.position for e in RAs]
    test_positions_pass = [1 if min_pos <= t <= max_pos else 0 for t in test_positions]

    evs_pass = [e for i,e in enumerate(RAs) if test_positions_pass[i]==1]
    test_alt = [e.alt for e in evs_pass]
    test_alt_minor = [e.alt_minor for e in evs_pass]

    #print genome[left:right], min_pos, max_pos, test_positions, test_positions_pass, test_alt, test_alt_minor
    # since this case is miscalled, we're not guaranteed exactly 1 deletion in this stretch of repeats. look for at least 1 '.' instead.
    if mutation.type == 'DEL':
	if '.' in test_alt or '.' in test_alt_minor:
	    fh.seek(fhstart)
	    mutation.evidence.append(ev)
	    return True
	else:
	    fh.seek(fhstart)
	    return False
    elif mutation.type == 'INS':
	if mutation.alt in test_alt or mutation.alt in test_alt_minor:
	    fh.seek(fhstart)
	    mutation.evidence.append(ev)
	    return True
	else:
	    fh.seek(fhstart)
	    return False

def findRepeat(mutation, fh, strain):
    # many of these are 1 bp deletions
    # relevant fields:
    # repeat_length
    # repeat_new_copies
    # repeat_ref_copies
    # repeat_seq
    # we can have shifted deletions:
    # ref: |ABABAB|
    # alt: |AB....|
    #      |..AB..|
    #      |....AB|
    # ref: |AAAAA| 5 copies
    # alt: |AAA..| 3 copies = 1 + (5-3) extra positions = 3 total positions
    #      |.AAA.|
    #      |..AAA|
    #       012345
    # ref: |..AA..| 2 copies
    # alt: |AAA...| 3 copies = 4 total positions
    #      |.AAA..|
    #      |..AAA.|
    #      |...AAA|
    # In theory we should require RA at one of these positions with major or minor allele = '.', and nowhere else
    
    fhstart = fh.tell() # get starting fh position if we need to rewind later
    line = fh.readline()
    ev = evidence(line.strip().split('\t'), allevkeys, strain)
    
    copies = max(int(mutation.repeat_ref_copies)*len(mutation.repeat_seq), int(mutation.repeat_new_copies)*len(mutation.repeat_seq))
    if mutation.type=='DEL':
	min_position = mutation.position - copies
    elif mutation.type=='INS':
	min_position = mutation.position - len(mutation.repeat_seq)*int(mutation.repeat_ref_copies) - len(mutation.repeat_seq)*int(mutation.repeat_new_copies) + 1
    
    if ev.position > mutation.position + copies:
	fh.seek(fhstart)
	return False

    # get rid of useless evidence lines
    while ev.position < min_position:
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    return False
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
   
    evs = [ev]
    while ev.position <= mutation.position + copies:
    #for i in range(copies):
	line = fh.readline()
	if line== '':
	    fh.seek(fhstart)
	    break
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
	evs.append(ev)
    RAs = [e for e in evs if e.evidenceType=='RA']

    test_positions = [e.position for e in RAs]
    test_positions_pass = [1 if min_position <= t <= min_position + copies else 0 for t in test_positions]
    
    change_len = abs(int(mutation.repeat_length) * (int(mutation.repeat_ref_copies) - int(mutation.repeat_new_copies)))
    
    if sum(test_positions_pass) != change_len: # need the exact number of calls in this region to match change_len
	fh.seek(fhstart)
	return False
    
    evs_pass = [e for i,e in enumerate(RAs) if test_positions_pass[i]==1]
    test_alt = ''.join([e.alt for e in evs_pass])
    test_alt_minor = ''.join([e.alt_minor for e in evs_pass])

    # because we expect no other calls in this region, the sequence in the region that passed should exactly match the number of '.' based on change_len
    if mutation.type == 'DEL':
	if test_alt == '.'*change_len or test_alt_minor == '.'*change_len:
	    fh.seek(fhstart)
	    mutation.evidence.append(ev)
	    return True
	else:
	    fh.seek(fhstart)
	    return False
    elif mutation.type == 'INS':
	if test_alt == mutation.repeat_seq*change_len or test_alt_minor == mutation.repeat_seq*change_len:
	    fh.seek(fhstart)
	    mutation.evidence.append(ev)
	    return True
	else:
	    fh.seek(fhstart)
	    return False

def findSUB(mutation, fh, strain):
    # this will try to match the exact number of basepairs that the mutation was called for. That means if the mutation is a sub of length 3, we will try to find 3 RAs that start and end at the same position.
    # after this, the function will REWIND the evidence file back to where it started, in case there are SNPs that overlap the SUB later!
    # ref: AAAAAAAAAA
    # alt: GGG
    # ref: |AAAAAAAAAA|
    # alt: |GGG.......|
    # alt: |.GGG......|
    # alt: |..GGG.....|
    # alt: |...GGG....|
    # alt: |....GGG...|
    # alt: |.....GGG..|
    # alt: |......GGG.|
    # alt: |.......GGG|
    fhstart = fh.tell()
    line = fh.readline()
    dist = len(line)
    ev = evidence(line.strip().split('\t'), allevkeys, strain)

    if ev.position > mutation.position:
	newdist = fh.tell()-dist
	fh.seek(newdist)
	return False
    
    # get right up to the sub
    while ev.position < mutation.position:
	line = fh.readline()
	dist = len(line)
	if line== '':
	    fh.close()
	    return False
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
    

    # now the next line either has to equal the start of the sub, or return false
    if not ev.position == mutation.position:
	newdist = fh.tell()-dist
	fh.seek(newdist)
	return False

    evs=[]
    # this first line matches in position.
    # sub must pass 2 tests: all pos must be consecutive, and alt must be in the set of ok alts
    # 1. dealing with deletions: construct set of ok alts
    diff = max([len(mutation.ref) - len(mutation.alt), 0])
    alt_match = set()
    for shift in range(diff+1):
	newmatch = ['.' for i in range(mutation.size)]
	newmatch[shift:shift+len(mutation.alt)] = list(mutation.alt)
	alt_match.add(''.join(newmatch))

    # 2. match for consecutive RAs
    consecutive_match = [i for i in range(mutation.position, mutation.position + mutation.size)]

    fhstart =  fh.tell()
    verboseprint('--------')
    verboseprint('file location before stepping:', fh.tell())

    # 3. grab the next mutation.size evs
    evs.append(ev)
    for i in range(1, mutation.size):
	line = fh.readline()
	if line== '':
	    fh.close()
	    fh.seek(fhstart)
	    verboseprint('file location after stepping:', fh.tell())
	    return False
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
	evs.append(ev)

    # 4. test for consecutive
    test_consec = [e.position for e in evs]
    test1 = test_consec == consecutive_match

    # 5. test for correct alt: need to look at both major and minor alleles--one of them must match
    test_alt = ''.join([e.alt for e in evs])
    test_alt_minor = ''.join([e.alt_minor for e in evs])
    test2 = (test_alt in alt_match) or (test_alt_minor in alt_match)
    
    fh.seek(fhstart)
    verboseprint(test_consec, consecutive_match, test_alt, test_alt_minor, alt_match)
    verboseprint('file location after stepping:', fh.tell())
    if test1 and test2:
	verboseprint('match!')
	mutation.evidence.extend(evs)
	return True
    else:
	verboseprint('nope')
	return False

def findRA(mutation, fh, strain):
    fhstart = fh.tell() # get starting fh position if we need to rewind later
    line = fh.readline()
    dist = len(line)
    ev = evidence(line.strip().split('\t'), allevkeys, strain)

    if ev.evidenceType != 'RA':
	fh.seek(fhstart)
	return False

    if ev.position > mutation.position:
	fh.seek(fhstart)
	return False

    while ev.position < mutation.position:
	line = fh.readline()
	if line== '':
	    fh.close()
	    return False
	dist = len(line)
	ev = evidence(line.strip().split('\t'), allevkeys, strain)
	if ev.evidenceType != 'RA':
	    fh.seek(-dist, 1)
	    return False

    if ev.position == mutation.position and (ev.alt == mutation.alt or ev.alt_minor == mutation.alt):
	mutation.evidence.append(ev)
	return True
    if ev.position > mutation.position:
	fh.seek(fhstart)
    return False


################################################################################################################################################
# relevant paths
runTag = argv[1]
mutations_path = runTag + '/postProcess/04merged/mutations.txt'
strains_path = runTag + '/postProcess/04merged/strains.txt'
matrix_path = runTag + '/postProcess/04merged/matrix.txt'
evpath_prefix = runTag + '/'
evpath_suffix = '/output/output.gd'

outpath = runTag + '/postProcess/06forceCalled/'
if not exists(outpath):
    makedirs(outpath)


evidence_types = set(['RA', 'MC', 'JC', 'UN'])
allevkeys = set(['side_1_read_count', 'left_inside_cov', 'left_outside_cov', 'coverage_minus', 'new_junction_coverage', 'prediction', 'ks_quality_p_value','frequency', 'variant_frequency', 'alignment_overlap', 'max_left_plus', 'bias_e_value', 'max_min_left', 'side_2_overlap', 'flanking_left', 'max_left', 'neg_log10_pos_hash_p_value', 'max_min_right_plus', 'total_non_overlap_reads', 'coverage_plus', 'major_frequency', 'fisher_strand_p_value', 'max_left_minus', 'side_2_redundant', 'side_2_coverage', 'side_1_continuation', 'side_1_possible_overlap_registers', 'minor_cov', 'side_1_redundant', 'side_2_possible_overlap_registers', 'right_outside_cov', 'new_junction_read_count', 'max_pos_hash_score', 'max_min_left_minus', 'side_2_continuation', 'pos_hash_score', 'polymorphism_score', 'junction_possible_overlap_registers', 'side_1_overlap', 'bias_p_value', 'max_min_right', 'flanking_right', 'max_right', 'side_1_coverage', 'right_inside_cov', 'max_right_plus', 'major_cov', 'side_2_read_count', 'max_min_right_minus', 'max_right_minus', 'consensus_score', 'total_cov', 'max_min_left_plus', 'ref_base', 'key'])

print '1) Loading input files'
# read mutations, remember position, ref, alt
mutationArr = []
with open(mutations_path) as f:
    header = f.readline()
    header = header.replace('\n', '') # need to be careful about using strip here, since later lines might be empty. Not a big deal for this line
    header = header.split('\t')

    for line in f.readlines():
	mutationArr.append(mutation(header, line))

    positions = [m.position for m in mutationArr]
    if positions != sorted(positions):
	print '----------------------------ERROR mutations are not sorted. This indicates a serious problem with generateMatrices.py and/or breseq gd files.'
	exit()

# read strains
with open(strains_path) as f:
    strains = f.readline().strip().split('\t')

# read matrix
with open(matrix_path) as f:
    strains_check = f.readline().strip().split('\t')
    for a,b in zip(strains, strains_check):
	if a != b:
	    print '----------------------------ERROR STRAIN MISMATCH CHECK GENERATEMATRICES.PY'
	    exit()
    matrix = parseMatrix(f)

# open file handle array for evidence.gd files
fharray = []
for st in strains:
    path = evpath_prefix + st + evpath_suffix
    if not exists(path):
	print '----------------------------ERROR evidence path cannot be found, possibly because of breseq error? Double check:'
	print path
	exit()
    fharray.append(open(path))

# skip comment lines for all files
print '2) Skipping comment and annotation lines in all gd files'
for fh in fharray:
    # skip the comment lines
    [fh, fh_done] = skipComments(fh)
    if fh_done:
	fh.close()
    [fh, fh_done] = skiptoevidence(fh)
    if fh_done:
	fh.close()

# read in genome.fasta
filename = '/groups/kishony/Reference_Genomes/SaureusNCTC8325/genome.fasta'
genome = []
with open(filename) as f:
    f.readline()
    for line in f.readlines():
	genome.append(line.strip())
genome = ''.join(genome)

print '3) Start force calling'
# start force calling
count = 0
print len(mutationArr)
print len(matrix)
for mutation, row in zip(mutationArr, matrix):
    forcecall(mutation, fharray, allevkeys)
    count +=1
print 'total processed:', count