import os,numpy

files = [
"0.000 Strained InGaAs.txt",
"0.100 Strained InGaAs.txt",
"0.200 Strained InGaAs.txt",
"0.240 Strained InGaAs.txt",]

p = os.path.split(__file__)[0]
for f in files:
    path = os.path.join(p,f)
    
    out_n = path.replace(".txt"," n.txt")
    out_p = path.replace(".txt"," k.txt.txt")
    nm,n,k = numpy.loadtxt(path, unpack = True)
    numpy.savetxt(out_n, numpy.array((nm,n)).T)
    numpy.savetxt(out_p, numpy.array((nm,k)).T)