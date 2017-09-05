import os,numpy

files = [f for f in os.listdir(".") if "InGaSb" in f]

for fn in files:
    n,data = numpy.loadtxt(fn, unpack=True)
    numpy.savetxt(fn, numpy.array((n*1e-9,data)).transpose())
