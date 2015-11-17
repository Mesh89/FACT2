import sys
import re;

n, k = int(sys.argv[1]), int(sys.argv[2])
print "xread", 10, n
print ""

for i in range(n):
	print "random_taxa_%d"%i, "0"*10
print ";"
print ""

print "tread ''"

trees = "*\n".join([line.split('(',1)[1].rsplit(')',1)[0] for line in open(sys.argv[3]) if line.startswith('tree')])

trees = trees.replace(":2.2222", "")
trees = trees.replace(",", " ")

trees = re.sub('\d+', 	lambda x:`int(x.group())-1`, trees)
print trees

print ";"
print ""

print "collapse-;"
print ""

print "freqdifs*;"

print "tread ''"
factres = sys.argv[4]
factres = factres.replace(":2.2222", "")
factres = factres.replace(",", " ")

factres = re.sub('\d+', lambda x:`int(x.group())-1`, factres)
print factres

print ";"
print ""

print "export- results.tre;"

print "quit;"


