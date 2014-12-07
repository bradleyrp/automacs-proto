#!/usr/bin/python

import os,re

bests = []
cwd = os.path.abspath(os.getcwd())
for root,modeldirs,filenames in os.walk('./',followlinks=True): break
for md in modeldirs:
	sign = re.findall('^model\-v[0-9]+\-(.+)',md)[0]
	with open(md+'/script-single.log','r') as fp: lines = fp.readlines()
	#---assume summary is found after "Summary of..." line and before "Total CPU time"
	lstart = [i for i,l in enumerate(lines) if re.search('Summary of successfully produced models',l)][0]
	lend = [i for i,l in enumerate(lines) if re.search('Total CPU time',l)][0]
	#---assume that DOPE score is in the third column
	rows = [l.strip('\n').split() for l in lines[lstart:lend] if re.match('.+\.pdb\s',l)]
	scores = [float(j[2]) for j in rows]
	bests.append([cwd+'/'+md+'/'+r[0] for i,r in enumerate(rows) if scores[i]==min(scores)][0])
with open('../repo/batch_file_list.txt','w') as fp:
	for b in bests: fp.write(b+'\n')
