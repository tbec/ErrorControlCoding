"""
Convolutional code helper functions
"""

# Imports
import numpy as np
import time
from math import sqrt
from numpy.random import randn
from numpy.random import seed
import sys
from matplotlib import pyplot as plt
#---------------------------------------------------

def bits(binVal):
    """
    Get number of bits in an integer
    :param binVal: integer input
    :return: bits in integer
    """
    return int(np.ceil(np.log2(binVal)))

def int2BinList(intVal):
    """
    Does octal conversion only
    :param intVal: int < 8 (oct)
    :return: binary list
    """

    if intVal>7: sys.exit()

    return [b for b in bin(intVal)[2:].zfill(3)]

def bin2Int(bitList):
    """
    Convert binary list to integer
    :param bitList: list of binary digits
    :return: integer equivelent
    """
    intVal = 0
    for bit in bitList:
        intVal = (intVal << 1) | bit
    return intVal

def stripLeadingZeros(binList):
    """
    Find number with fewest leading zeros and ret
    :param binList:
    :return: binList with leading zeros stripped
    """
    ind = 0
    gba = np.asarray(binList).astype(int)
    while not gba[:,ind].sum(): ind+=1
    return gba[:,ind:]

def intStr2BinArray(intBlock):
    """
    Converts octal integers into binary array
    :param intStr:
    :return: binary list of octals
    """

    # Ints to int list
    ints = [int(i) for i in str(intBlock)]

    binList = []

    # List of Strings
    for i in ints:
        binList += int2BinList(i)

    #Convert to ints and return
    return map(int,binList)

def mod2Mult(A,B):
    """
    Matrix multiply then mod 2, corresponding
    to binary addition with no carry.
    :param A: first matrix
    :param B: second matrix
    :return: mod 2 matrix multiply of A,B
    """
    return np.mod(np.dot(A,B),2)

def mod2Convolve(A,B):
    return np.mod(np.convolve(A,B),2)

def parseCodeword(v,n):
    """
    Takes output of v = uG and converts
    into individual codewords
    :param v: codeword input
    :param n: number of outputs
    :return: list of codewords
    """
    return [v[i::n] for i in xrange(n)]

def int2bitArray(intVal,bits):
    """
    Convert integer to binary number
    :param intVal: integer number
    :param bits: number of bits (cannot truncate)
    :return: numpy binary array (LSB at arr[0])
    """
    binary = bin(intVal)[2:].zfill(bits)[::-1]
    return np.asarray([int(i) for i in binary])

def int2EsArray(intVal,bits):
    temp = int2bitArray(intVal,bits)
    temp[np.where(temp == 0)] = -1
    return temp

def bitArray2int(bitArray):
    """
    [1 1 0] is a 3
    [0 1 1] is a 6
    :param bitArray:
    :return:
    """
    intVal = 0
    for i in xrange(max(bitArray.shape)):
        intVal += bitArray[i]<<i
    return intVal

def hamming_distance(A,B):
    return np.bitwise_xor(A,B).sum()

def hamming_int(Aint,Barr):
    a = int2bitArray(Aint,len(Barr))
    return hamming_distance(a,Barr)

def poly2trellis2(nkvTuple,GenPolys,recursive=False):
    """
    Create the trellis structure for decoding.
    :return:
    """

    n = nkvTuple[0]
    k = nkvTuple[1]
    v = nkvTuple[2]

    g = np.asarray([intStr2BinArray(i) for i in GenPolys],int)
    g = stripLeadingZeros(g)

    num_input_symbols = 2**k
    num_states = 2**v

    # State/Output Table Initialization
    state_table = np.zeros((num_input_symbols,num_states),int)
    output_table = state_table.copy()

    # Loop through states
    for state in xrange(num_states):

        bit_state = int2bitArray(state,v)

        # Loop through inputs (0,1) for k=1
        for ipt in xrange(num_input_symbols):

            if recursive:

                # compute feedback
                fb = (ipt+(bit_state&g[0,1:]).sum())%2

                # Compute next state
                state_table[ipt,state] = (state*2+fb) % 2**v

                # rest of v's
                bit_output = [((fb&g[1,0])+(bit_state&g[1,1:]).sum())%2]
                # bit_output = [((fb&g[i,0])+(bit_state&g[i,1:]).sum())%2 for i in xrange(1,n)]

                # insert v0 - systematic
                bit_output.insert(0,ipt)

            else:
                # Compute next state
                state_table[ipt,state] = (state*2+ipt) % 2**v

                # Compute output
                bit_output = [(((ipt&g[i,0])+(bit_state&g[i,1:]).sum())%2) for i in xrange(n)]

            output_table[ipt,state] = bitArray2int(np.asarray(bit_output,int))

    return state_table, output_table

def poly2trellis(nkvTuple,GenPolys,recursive=False):
    """
    Create the trellis structure for decoding.
    :return:
    """

    n = nkvTuple[0]
    k = nkvTuple[1]
    v = nkvTuple[2]

    g = np.asarray([intStr2BinArray(i) for i in GenPolys],int)
    g = stripLeadingZeros(g)

    num_input_symbols = 2**k
    num_states = 2**v

    # State/Output Table Initialization
    state_table = np.zeros((num_input_symbols,num_states),int)
    output_table = state_table.copy()

    # Loop through states
    for state in xrange(num_states):

        bit_state = int2bitArray(state,v)

        # Loop through inputs (0,1) for k=1
        for ipt in xrange(num_input_symbols):

            if recursive:

                # compute feedback
                fb = (ipt+(bit_state&g[0,1:]).sum())%2

                # rest of v's
                bit_output = [((fb&g[i,0])+(bit_state&g[i,1:]).sum())%2 for i in xrange(1,n)]

                # insert v0 - systematic
                bit_output.insert(0,ipt)

                # Compute next state
                state_table[ipt,state] = (state*2+fb) % 2**v

            else:
                # Compute next state
                state_table[ipt,state] = (state*2+ipt) % 2**v

                # Compute output
                bit_output = [(((ipt&g[i,0])+(bit_state&g[i,1:]).sum())%2) for i in xrange(n)]

            output_table[ipt,state] = bitArray2int(np.asarray(bit_output,int))

    return state_table, output_table

def maxs(logList):
    # Check for single element
    if len(logList)==1:
        return logList[0]

    # Check for 1 (and NaN) or 2 elements
    if len(logList)==2:
        if   np.isnan(logList[0]): return logList[1]
        elif np.isnan(logList[1]): return logList[0]
        else: return maxs2(logList[0],logList[1])

    # Otherwise start multi-element max* algorithm
    else:
        a = logList.pop()
        return max_rec(a,logList)

def max_rec(a,logList):
    if len(logList) > 1:
        return maxs2(a,max_rec(logList.pop(),logList))
    else:
        return maxs2(a,logList[0])

def maxs2(a,b):
    return max(a,b) + np.log(1+np.exp(-abs(a-b)))

def parentStates(stateTable):
    """
    return a matrix of a states parent states and inputs
    For k=1, there can only be two parent states for each state
    Format of matrix is:
           ----------------
          /    state n    /|
         /---------------/ |
        /      ...      /| |
       /---------------/ | |
      /   state 1     /| | /
     / --------------/ | |/
    /    state 0    /| | /
   |----------------|| |/
   | parent | input || /
   |--------|-------||/
   | parent | input |/
    ----------------
    """
    _, states = stateTable.shape
    parentTable = np.zeros((2,2,states),int)

    for state in xrange(states):

        # get indices of parents: indices are (input,parent state)
        parents = np.asarray(np.where(stateTable==state))
        parentTable[:,:,state] = parents

    return parentTable

def toBits(arr):
    arr = np.sign(arr).astype(int)
    arr[np.where(arr == -1)] = 0
    return arr

def BPSKencode(msg,Eb_N0):
    msg[np.where(msg==0)]=-1
    Es = sqrt(2*Eb_N0)
    return msg*Es

def BPSK_AWGN(x,Eb_N0_linear):
    """
    add AWGN to a channel
    :param msg: message to add AWGN to
    :param SNR: SNR in dB
    :return:    AWGN signal
    """

    x[np.where(x==0)]=-1
    #sigma = sqrt(1/(2*n*Eb_N0_linear))
    sigma = sqrt(1/(Eb_N0_linear))
    noise = sigma*randn(max(x.shape))
    y = x + noise
    #
    # red = x[np.where(x==1)]
    # red = red + randn(0,sigma,max(red.shape))
    # blue = x[np.where(x==-1)]
    # blue = blue + randn(0,sigma,max(blue.shape))
    # plt.scatter(np.linspace(0,1,max(blue.shape)),blue)
    # plt.scatter(np.linspace(0,1,max(red.shape)),red,color="red")
    # plt.show()

    return y



def BER(u,u_hat):
    """
    Calculate bit error rate
    :param u: encoded message
    :param v: decoded message ( can include terminating bits )
    :return: (total weight (d), BER (weight/length))
    """
    L = max(u.shape)
    diff = np.subtract(u,u_hat[:L])        # Find the difference
    d = np.sum(np.abs(diff))               # Calculate total weight
    ber = d/float(L)                       # BER is ratio of weight/length
    return d,ber

def lin(dB):
    return 10**(dB/10.0)

def dB(lin):
    return 10*np.log10(lin)