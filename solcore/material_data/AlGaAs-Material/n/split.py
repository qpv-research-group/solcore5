# import os,numpy
#
# files = ["0.000 AlGaAs.txt",
# "0.099 AlGaAs.txt",
# "0.198 AlGaAs.txt",
# "0.315 AlGaAs.txt",
# "0.419 AlGaAs.txt",
# "0.491 AlGaAs.txt",
# "0.590 AlGaAs.txt",
# "0.700 AlGaAs.txt",
# "0.804 AlGaAs.txt",
# "1.000 AlGaAs.txt",]
#
# p = os.path.split(__file__)[0]
# for f in files:
#     path = os.path.join(p,f)
#
#     out_n = path.replace(".txt"," n.txt")
#     out_p = path.replace(".txt"," k.txt.txt")
#     nm,n,k.txt = numpy.loadtxt(path, unpack = True)
#     numpy.savetxt(out_n, numpy.array((nm,n)).T)
#     numpy.savetxt(out_p, numpy.array((nm,k.txt)).T)


import numpy
from solcore3 import *
from solcore3.graphing import *


x,y = numpy.loadtxt("../k.txt/0.000 AlGaAs k.txt.txt", unpack=True)

Graph(GraphData(x,y), edit=lambda x,y: (x*1e9,y), yscale="log").draw()