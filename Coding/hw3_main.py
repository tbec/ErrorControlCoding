from utils import *
import decoders

print
# ----------------------------------------#
# ------ BCJR Test case from notes ------ #
# --------------------------------------- #
# nkvTuple = (2,1,1)
# gens = (3,2)
# r = [0.8, 0.1, 1, -0.5, -1.8, 1.1, 1.6, -1.6]
# Es_N0 = 0.25
#
# bcjr = decoders.BCJR(nkvTuple,gens,recursive=True)
# msg = bcjr.decode(r,Es_N0=Es_N0)
# LLR = bcjr.last_decoded_msg_energy
# print 'LRR:',LLR
# print 'msg:',msg


# --------------------------------------- #
# ------ BCJR Email test ---------------- #
# --------------------------------------- #
nkvTuple = (2,1,2)
gens = (7,5)
r = [0.8, -0.6, 1.2, -0.5, 2, -1, -1.5, 2.0, 1.3, -0.7]
Es_N0 = 0.2
bcjr = decoders.BCJR(nkvTuple,gens,recursive=True)
msg = bcjr.decode(r,Es_N0=Es_N0)
LLR = bcjr.last_decoded_msg_energy
print 'LRR:',LLR
print 'msg:',msg


# --------------------------------------- #
# ------ BCJR Email test (2) ------------ #
# --------------------------------------- #
nkvTuple = (2,1,2)
gens = (7,5)
r = [0.8, -0.6, 1.2, -0.5, 2, -1, -1.5, 2.0, 1.3, -0.7]
La = [2, -1, 0, 0.5, -0.7]
Es_N0 = 0.2
bcjr = decoders.BCJR(nkvTuple,gens,recursive=True)
msg = bcjr.decode(r,La=La,Es_N0=Es_N0)
LLR = bcjr.last_decoded_msg_energy
print 'LRR:',LLR
print 'msg:',msg


# --------------------------------------- #
# ------ 12.30(a) log-MAP --------------- #
# --------------------------------------- #
# LRR: [  1.72898863   1.5267149    1.51964064   0.22103251   7.77038962
#         7.16099204   3.48839034   2.91098357 -42.54640126]
# msg: [1 1 1 1 1 1 1 1 0]
#
# nkvTuple = (2,1,1)
# gens = (2,3)
# r = [1.5339,0.639,-0.6737,-3.0183,1.5096,
#      0.7664,-0.4019,0.3185,2.7121,-0.7304,
#      1.4169,-2.0341,0.8971,-0.3951,1.6254,
#      -1.1768,2.6954,-1.0575]
# Es_N0 = 0.5
# bcjr = decoders.BCJR(nkvTuple,gens,recursive = True)
# LRR,msg = bcjr.decode(r,Es_N0=Es_N0)
# print 'LRR:',LRR
# print 'msg:',msg


# --------------------------------------- #
# ------ 12.30(b) Max-log-MAP ----------- #
# --------------------------------------- #
# LRR: [  1.46560000e+00   1.46560000e+00   8.49400000e-01   2.00000000e-02
#         7.60740000e+00   7.60740000e+00   3.72020000e+00   3.48940000e+00
#        -4.24942000e+01]
# msg: [1 1 1 1 1 1 1 1 0]


# --------------------------------------- #
# ------ (3) ---------------------------- #
# --------------------------------------- #
# LRR: [  2.74557952e+00   2.30681009e+00  -2.83605009e-02   6.04598533e-01
#         7.90506018e+00   6.25086130e+00   4.03494233e+00   1.58764800e+00
#        -4.26772162e+01]
# msg: [1 1 0 1 1 1 1 1 0]
#
# nkvTuple = (2,1,1)
# gens = (2,3)
# r = [1.5339,0.639,-0.6737,-3.0183,1.5096,
#      0.7664,-0.4019,0.3185,2.7121,-0.7304,
#      1.4169,-2.0341,0.8971,-0.3951,1.6254,
#      -1.1768,2.6954,-1.0575]
# La = [1,-0.5,-1.5,0,0.8,-1.2,2,-1.8]
# Es_N0 = 0.5
# bcjr = decoders.BCJR(nkvTuple,gens,recursive = True)
# LRR,msg = bcjr.decode(r,La=La,Es_N0=Es_N0)
# print 'LRR:',LRR
# print 'msg:',msg


# --------------------------------------- #
# ------ (4,5) simulations.py ----------- #
# --------------------------------------- #

# --------------------------------------- #
# ------ Other test stuff --------------- #
# --------------------------------------- #
print
