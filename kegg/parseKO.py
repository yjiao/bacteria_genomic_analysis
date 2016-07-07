import re

A = ''
B = ''
C = ''

fid = open('sao_parsed.tsv', 'w')
with open('sao00001.keg','r') as f:
    for line in f.readlines():
	if line[0] != '#' and line[0] !='+' and line[0] !='%' and line[0] !='!':
	    line = re.sub(r'<[^>]*>', '\t', line)
	    line = line.strip()

	    if len(line) == 0:
		continue

	    # this file has great formatting
	    level = line[0]
	    
	    if level == 'A':
		A = line[1:].strip()
	    elif level == 'B':
		B = line[1:].strip()
	    elif level == 'C':
		C = line[1:].strip()
	    else:
		if level != 'D':
		    print 'ERROR: WRONG LEVEL??'
		    continue

		lineorig = line
		line = line[1:]
		line = line.strip().split(' ')
		locus = line[0]
		
		line = ' '.join(line[1:])
		line = line.split('\t')
		product = line[0]
		
		line = ' '.join(line[1:])
		line = line.split(' ')
		KO = line[0]
		
		if len(KO) == 0:
		    gene = 'na'
		    KO = 'na'
		    descr = 'na'
		else:
		    line = ' '.join(line[1:])
		    line = line.split(';')
		    gene = line[0]
		    try:
			descr = line[1].strip()
		    except:
			descr = 'na'


		fid.write('\t'.join([locus, gene, product, descr, KO, C, B, A]))
		fid.write('\n')




		