# -*- coding: iso-8859-1 -*-
import sys, os

filename=sys.argv[1]
f = open(filename, 'r')
out = 'Kernel'
g = open(out, 'w')
buff=''
Nline=0

for line in f:
	buff=line[0:len(line)-1]+str('\\n\\')+'\n'
	g.write(buff)
	Nline+=1
g.close()
print 'number of lines'
print Nline
print 'number of characters'
g.close()
g = open(out, 'r')
buff=g.read()
Nchar=len(buff)-3*Nline
print Nchar
g.close()