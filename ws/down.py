#!/usr/bin/env python
import logging, os, sys, json, zipfile, uuid, urllib

logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
logging.debug('START_SCRIPT:down.py')

# init output
cwd      = os.getcwd()
datapath = cwd + '/../data'
tmppath  = cwd + '/../tmp'
outfname = str(uuid.uuid4()) + '.zip'
outfile  = tmppath + '/' + outfname

# get input parameters
# Warning: we don't check if is defined
inparam = json.load(sys.stdin)

# zip input files
with zipfile.ZipFile(outfile, 'w') as myzip:
	for indat in inparam:
		# check data attributes
		if 'raw' not in indat:
			raise ValueError("'raw' param does not exist")
		if 'vseq' not in indat:
			raise ValueError("'vseq' param does not exist")
		if 'path' not in indat['vseq']:
			raise ValueError("'vseq.path' param does not exist")
		if 'imgfile' not in indat['vseq']:
			raise ValueError("'vseq.imgfile' param does not exist")
		
		# check dirs and files
		os.chdir(datapath)
		indir = indat['vseq']['path']
		if not os.path.exists(indir):
			raise ValueError("'vseq.path' directory does not exist")
		infile = indir + '/' + indat['vseq']['imgfile']
		if not os.path.isfile(infile):
			raise ValueError("'imgfile' file does not exist")
		
		# zip files
		myzip.write(infile)

myzip.close()

# print response
print "Content-type: application/json\n"
result = {
	'success':'true',
	'message':'The Command Completed Successfully',
	'data': 'tmp/'+outfname
}
print json.dumps(result)
