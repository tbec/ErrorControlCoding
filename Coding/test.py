import numpy as np
from utils import *
from encoders import *
from decoders import *
from simulations import msgTxRx

file = open('./test_file.txt','w+')
file.write("This is a file test.\n")
file.close()

# u = np.random.randint(0,2,2000)
#
# nkvTuple = (2,1,2)
# gens = (5,7)
#
# for i in xrange(10000):
#     print msgTxRx(2000,nkvTuple,gens,"NSFFCC",lin(2))

# encoder = SRCC(nkvTuple,gens)
# v = encoder.encode(u)
#
# v[np.where(v==0)] = -1
#
# v_bpsk = BPSKencode(v,lin(4))
#
# decoder = BCJR(nkvTuple,gens,recursive=True)
# u_hat = decoder.decode(v_bpsk)
#
# _,PbE = BER(u,u_hat)

#print PbE