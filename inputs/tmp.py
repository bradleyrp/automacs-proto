#!/usr/bin/python -i

import re
import json

blank = {}
regex = '^([^\s]+)\s?=\s?([^\s]+)\s?$'
with open('../sources/cgmd-bilayer-construct/input-em-steep-in.mdp','r') as fp:
	for line in fp:
		 if re.match(regex,line):
		 	a = re.findall(regex,line)[0]
		 	blank[a[0]] = a[1]
print re.sub(r'"','\'',json.dumps(blank,indent=4,separators=(',',':')))

