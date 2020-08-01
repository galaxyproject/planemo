# reverse order of columns in a tabular file
import argparse
parser = argparse.ArgumentParser()
a = parser.add_argument
a('--infile',default='')
a('--outfile',default=None)
a('--prefix',default=None)
args = parser.parse_args()
inp = args.infile
outp = args.outfile
i = open(inp,'r').readlines()
o = open(outp,'w')
for row in i:
	rs = row.rstrip()
	rs = list(rs)
	rs.reverse()
	o.write('%s:%s' % (args.prefix,''.join(rs)))
	o.write('\n')
o.close()
