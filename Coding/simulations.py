import numpy as np
import sys
from matplotlib import pyplot as plt
from utils import *
from encoders import *
from decoders import *
import Queue
import threading

# Simulate (5,7) and (27,31) codes for
#   nonsystematic nonrecursive and for
#   systematic recursive. Information
#   blocks are 2000 bits, BPSK modulated,
#   then AWGN is added according to SNR.
#   r is then demodulated. Bit error is
#   calculated for various SNR values.
def main():


    # Eb/N0 to loop over
    Eb_N0 = np.arange(0,4.25,0.25)
    Eb_N0_lin = lin(Eb_N0)

    # Number of iterations per Eb/N0
    ITER_PER_SNR = 500
    NUM_INFO_BITS = 2000

    # nkv = [(2,1,2),(2,1,2),(2,1,4),(2,1,4)]
    # gens = [(5,7),(5,7),(27,31),(27,31)]
    # encoders = ['NSFFCC','SRCC','NSFFCC','SRCC']
    # file_names = []
    # for i in xrange(len(nkv)):
    #     file_names.append('./dat/%s_%s'%(encoders[i],"_".join(map(str,gens[i]))))
    #
    # q = Queue.Queue()
    #
    # for i in xrange(len(nkv)):
    #     args = (ITER_PER_SNR,
    #             NUM_INFO_BITS,
    #             encoders[i],
    #             file_names[i],
    #             nkv[i],
    #             gens[i],
    #             i)
    #     t = threading.Thread(target=estimateCodeBER,args=args)
    #     t.daemon = False
    #     t.start()


    PbE_vs_SNR_NSFFCC_5_7 = []
    PbE_vs_SNR_SRCC_5_7 = []
    PbE_vs_SNR_NSFFCC_27_31 = []
    PbE_vs_SNR_SRCC_27_31 = []

    nkvTuple = (2,1,2)
    gens = (5,7)

    # ---- NSFFCC (5,7) ---- #

    print '\nNSFFCC (5,7)'

    for EbN0linear in Eb_N0_lin:
        print '\tEb/N0(dB):',dB(EbN0linear)
        PbE_avg = []
        for i in xrange(ITER_PER_SNR):
            PbE = msgTxRx(NUM_INFO_BITS,nkvTuple,gens,"NSFFCC",EbN0linear)
            PbE_avg.append(PbE)
        PbE_vs_SNR_NSFFCC_5_7.append(sum(PbE_avg)/len(PbE_avg))

    np.save('./PbE_vs_SNR_NSFFCC_5_7',np.array(PbE_vs_SNR_NSFFCC_5_7))

    # ---- SRCC (5,7) ---- #

    print '\nSRCC (5,7)'

    for EbN0linear in Eb_N0_lin:
        print '\tEb/N0(dB):',dB(EbN0linear)
        PbE_avg = []
        for i in xrange(ITER_PER_SNR):
            PbE = msgTxRx(NUM_INFO_BITS,nkvTuple,gens,"SRCC",EbN0linear)
            PbE_avg.append(PbE)
        PbE_vs_SNR_SRCC_5_7.append(sum(PbE_avg)/len(PbE_avg))

    np.save('./PbE_vs_SNR_SRCC_5_7',np.array(PbE_vs_SNR_SRCC_5_7))

    # ---- NSFFCC (27,31) ---- #

    nkvTuple = (2,1,4)
    gens = (27,31)

    print '\nNSFFCC (27,31)'

    for EbN0linear in Eb_N0_lin:
        print '\tEb/N0(dB):',dB(EbN0linear)
        PbE_avg = []
        for i in xrange(ITER_PER_SNR):
            PbE = msgTxRx(NUM_INFO_BITS,nkvTuple,gens,"NSFFCC",EbN0linear)
            PbE_avg.append(PbE)
        PbE_vs_SNR_NSFFCC_27_31.append(sum(PbE_avg)/len(PbE_avg))

    np.save('./PbE_vs_SNR_NSFFCC_27_31',np.array(PbE_vs_SNR_NSFFCC_27_31))

    # ---- SRCC (27,31) ---- #

    print '\nSRCC (27,31):'

    for EbN0linear in Eb_N0_lin:
        print '\tEb/N0(dB):',dB(EbN0linear)
        PbE_avg = []
        for i in xrange(ITER_PER_SNR):
            PbE = msgTxRx(NUM_INFO_BITS,nkvTuple,gens,"SRCC",EbN0linear)
            PbE_avg.append(PbE)
        PbE_vs_SNR_SRCC_27_31.append(sum(PbE_avg)/len(PbE_avg))

    np.save('./PbE_vs_SNR_SRCC_27_31',np.array(PbE_vs_SNR_SRCC_27_31))


# ------------------------------------------------- #
# ----------------- Methods ----------------------- #
# ------------------------------------------------- #
ENCODER_ERROR_MSG = \
"""
Encoder parameter not recognized.
Please choose a supported encoder type:
    - NSFFCC
    - SRCC
"""
DECODER_ERROR_MSG = \
"""
Decoder parameter not recognized.
Please choose a supported decoder type:
    - BCJR
    - Viterbi
"""

def estimateCodeBER(ITER_PER_SNR,NUM_INFO_BITS,encoder,file_name,nkvTuple,generator_polynomials,threadNum):

    # Eb/N0 to loop over
    Eb_N0 = np.arange(0,4.25,0.25)
    Eb_N0_lin = lin(Eb_N0)

    arr = []

    for EbN0linear in Eb_N0_lin:
        PbE_avg = []
        for i in xrange(ITER_PER_SNR):
            if i%10==0:
                print ('-'*50)+'\nThread: %d\tEncoder: %s\tGens: (%d,%d)\n\tEb/N0: %d\tIteration: %d'% \
                    (threadNum,encoder,generator_polynomials[0],generator_polynomials[1],dB(EbN0linear),i)
            PbE = msgTxRx(NUM_INFO_BITS,nkvTuple,generator_polynomials,encoder,EbN0linear)
            PbE_avg.append(PbE)
        avg = sum(PbE_avg)/len(PbE_avg)
        arr.append(avg)

    np.save(file_name,np.array(arr))

def msgTxRx(msg_size,nkvTuple,gens,encoder,Eb_N0):

    if encoder=="NSFFCC":
        rec = False
        encoder = NSFFCC(nkvTuple,gens)
    elif encoder=="SRCC":
        rec = True
        encoder = SRCC(nkvTuple,gens)
    else: sys.exit()

    decoder = BCJR(nkvTuple,gens,recursive=rec)

    u = np.random.randint(0,2,msg_size)

    v = encoder.encode(u)

    v_bpsk_AWGN = BPSK_AWGN(v,Eb_N0)

    decoder = BCJR(nkvTuple,gens,recursive=rec)
    u_hat = decoder.decode(v_bpsk_AWGN)

    _,PbE = BER(u,u_hat)

    return PbE


if __name__=='__main__':
    main()
