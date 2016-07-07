# This script parses breseq outputs for each strain in the dataset and collects them in matrices where rows are mutations and columns are strains.

# NOTES
# 2016 01 26 fixed mutations.txt output bug where strains were printed as objects intead of names

from os import listdir, makedirs
import heapq as h
import fnmatch
from os.path import isfile, join, exists
from collections import defaultdict
from sys import argv

# classes
class mutation():
    def __init__(self, annotList, genome, allkeys):
	for key in allkeys:
	    self.__dict__[key] = ''

	self.type = annotList[0]
	self.position = int(annotList[4])

	if self.type == 'SNP':
	    self.alt = annotList[5]
	    self.size = 0
	elif self.type == 'SUB':
	    self.alt = annotList[6]
	    self.size = int(annotList[5])
	elif self.type == 'DEL':
	    self.alt = '.'
	    self.size = int(annotList[5])
	elif self.type == 'INS':
	    self.alt = annotList[5]
	    self.size = len(self.alt)
	
	start = self.position - 1
	if self.size == 0:
	    stop = start + 1
	else:
	    stop = start + self.size
	self.ref = genome[start:stop]

	self.ID = '_'.join([str(self.position), self.ref, self.alt])
	self.strains = []
	
	# add all other annotations
	if self.type == 'SNP':
	    start = 6
	elif self.type == 'SUB':
	    start = 7
	elif self.type == 'DEL':
	    start = 6
	elif self.type == 'INS':
	    start = 6

	for i in range(start, len(annotList)):
	    pair = annotList[i].split('=')
	    if len(pair) != 2:
		continue
	    key = pair[0]
	    val = pair[1]
	    if key not in allkeys:
		continue
	    self.__dict__[key] = val

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
    def __repr__(self):
	return self.name


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

def write_matrix(outname, mut, strains):
    with open(outname, 'a') as f:
	newline = [0 for s in strains]
	for strain in mut.strains:
	    newline[strains.index(strain.name)] = 1
	
	f.write('\t'.join(map(str, newline))+'\n')

def write_mutations(outname, mut, allkeys):
    with open(outname, 'a') as f:
	newline = []
	for key in allkeys:
	    if key == 'strains':
		strains = []
		for strain in mut.strains:
		    strains.append(strain.name)
		newline.append(','.join(strains))
	    else:
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
	    newline = ['' for s in strains]
	    for strain in mut.strains:
		field = ','.join(strain.__dict__[param])
		newline[strains.index(strain.name)] = field
	    f.write('\t'.join(newline)+'\n')

def write_index(outname, mut, strains):
    with open(outname, 'a') as f:
	newline = ['' for s in strains]
	for strain in mut.strains:
	    newline[strains.index(strain.name)] = strain.index
	
	f.write('\t'.join(map(str, newline))+'\n')

################################################################################################################################################
# get annotated files--snps, dels, etc only
base = argv[1]
basepath = base + '/postProcess/'
path = basepath+'05annotEvMerged/'
#annotfiles = [f for f in listdir(path) if isfile(join(path, f))][0:10]
annotfiles = [f for f in listdir(path) if isfile(join(path, f))]
strains = [f.replace('.gd', '') for f in annotfiles]
evpath_prefix = base + '/'#'fullRun20151221/'
evpath_suffix = '/output/evidence/'

# outfile names
outMutationsPath = basepath + '04merged/mutations.txt'
strainsPath = basepath + '04merged/strains.txt'
outMatrixPath = basepath + '04merged/matrix.txt'
outIndexPath = basepath + '04merged/indexhtml.txt'
outHTMLPath = basepath + '04merged/evidencehtml.txt'
if not exists(basepath+'04merged'):
    makedirs(basepath + '04merged')

# read in genome.fasta
filename = '/groups/kishony/Reference_Genomes/SaureusNCTC8325/genome.fasta'
genome = []
with open(filename) as f:
    f.readline()
    for line in f.readlines():
	genome.append(line.strip())

genome = ''.join(genome)

print '-----------------'
print len(annotfiles), 'annotfiles will be opened for merging'

# create tuples of annot and norm files and store in fharray
print '1) Creating tuples of normalized, annotated, and evidence files'
fharray = []
for s, annotfile in zip(strains, annotfiles):
    path = basepath + '05annotEvMerged/' + s + '.gd'
    fharray.append((open(path)))

# skip comment lines for all files
print '2) Skipping comment lines in all files'
for fh in fharray:
    # skip the comment lines
    [fh, fh_done, comments] = skipComments(fh)
    if fh_done:
	fh.close()

print '3) main loop'
Q = []
# all keys used
allkeys = set(['snp_type', 'aa_ref_seq', 'alt', 'size', 'locus_tag', 'repeat_length', 'gene_strand', 'gene_product', 'codon_position', 'ref', 'gene_list', 'aa_position', 'codon_ref_seq', 'repeat_ref_copies', 'ID', 'gene_position', 'type', 'codon_new_seq', 'strains', 'aa_new_seq', 'position', 'repeat_new_copies', 'gene_name', 'repeat_seq'])
allevkeys = set(['side_1_read_count', 'left_inside_cov', 'left_outside_cov', 'coverage_minus', 'new_junction_coverage', 'prediction', 'ks_quality_p_value','frequency', 'variant_frequency', 'side_2_annotate_key', 'alignment_overlap', 'max_left_plus', 'bias_e_value', 'max_min_left', 'side_2_overlap', 'flanking_left', 'max_left', 'neg_log10_pos_hash_p_value', 'max_min_right_plus', 'total_non_overlap_reads', 'coverage_plus', 'major_frequency', 'fisher_strand_p_value', 'side_1_annotate_key', 'max_left_minus', 'side_2_redundant', 'side_2_coverage', 'side_1_continuation', 'side_1_possible_overlap_registers', 'minor_cov', 'side_1_redundant', 'side_2_possible_overlap_registers', 'right_outside_cov', 'new_junction_read_count', 'max_pos_hash_score', 'max_min_left_minus', 'side_2_continuation', 'pos_hash_score', 'polymorphism_score', 'junction_possible_overlap_registers', 'side_1_overlap', 'bias_p_value', 'max_min_right', 'flanking_right', 'max_right', 'side_1_coverage', 'right_inside_cov', 'max_right_plus', 'major_cov', 'side_2_read_count', 'max_min_right_minus', 'max_right_minus', 'consensus_score', 'total_cov', 'max_min_left_plus'])


# write headers to the mutation out file
with open(strainsPath, 'w') as f:
    header = '\t'.join(strains)
    f.write(header)
with open(outMutationsPath, 'w') as f:
    header = '\t'.join(allkeys)+ '\n'
    f.write(header)
# write headers to the matrix out file
with open(outMatrixPath, 'w') as f:
    header = '\t'.join(strains)+ '\n'
    f.write(header)
# write headers to the evidence url out file
with open(outHTMLPath, 'w') as f:
    header = '\t'.join(strains)+ '\n'
    f.write(header)
# write headers to the index url out file
with open(outIndexPath, 'w') as f:
    header = '\t'.join(strains)+ '\n'
    f.write(header)

h.heapify(Q)
allmutations = set() # paranoid checking
mutdict = {}
while True:
    # get one line from each file
    for strain, fh in zip(strains,fharray):
	if fh.closed:
	    continue

	line = fh.readline()
	if line == '':
	    fh.close()
	    continue
	
	# create mutation object and insert it into dictionary using mut ID
	line = line.strip().split('\t')
	mut = mutation(line, genome, allkeys)
	ID = mut.ID
	if mut.ID not in mutdict:
	    # this mutation is new, so it must not be in the queue
	    if mut.ID in allmutations:
		print 'ERROR!!!! INVARIANT VIOLATED'
	    mutdict[ID] = mut
	    h.heappush(Q, mutdict[ID])
	else:
	    # this mutation must already be in the queue, ignore
	    del mut
	allmutations.add(ID)

	##### find evidence url
	url = evpath_prefix + strain + evpath_suffix + mutdict[ID].type + '_' + line[1] + '.html' # this is not valid for missing coverage
	evlinks = line[2].split(',')
	if not isfile(url):
	    url = []
	    path = evpath_prefix + strain + evpath_suffix
	    for link in evlinks:
		matches = fnmatch.filter(listdir(path), '*'+link+'*.html')
		keep = [i for i in range(len(matches)) if matches[i].find('SIDE') == -1]
		url.extend([path + matches[k] for k in keep])
	else:
	    url = [url]
#	if len(url) == 0:
#	    print strain, line
	# double check urls are correct
	for u in url:
	    if not isfile(u):
		print u, line

	### find index.html url
	indexurl = evpath_prefix + strain + '/output/index.html'
	if not isfile(indexurl):
	    indexurl = 'none'
	
	# get evidence field
	evidence = filter(lambda x: x.find('evidence') != -1, line)
	if len(evidence) != 1:
	    print '----------------------------weird evidence entry???'
	evidence = evidence[0]
	evidence = evidence.replace('evidence=','')
	evidence = evidence.split('@@@')
	evidence = [evline.split(';') for evline in evidence]

	# add strain to mutation object in mutdict
	col = entry(url, indexurl, strain, evidence, allevkeys)
	mutdict[ID].strains.append(col)
    
    earliest_mutation = h.heappop(Q)
    write_matrix(outMatrixPath, earliest_mutation, strains)
    write_mutations(outMutationsPath, earliest_mutation, allkeys)
    write_url(outHTMLPath, earliest_mutation, strains)
    write_index(outIndexPath, earliest_mutation, strains)
    write_evidence_params(earliest_mutation, strains, allevkeys, basepath + '04merged/')

    mutdict.pop(earliest_mutation.ID)


    if len(mutdict) == 0:
	break

print 'All files closed:', all([fh.closed for fh in fharray])