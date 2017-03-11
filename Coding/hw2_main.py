import encoders,decoders
import numpy as np

class_header = \
"""
Thomas Becnel
u0668404
ECE 6521 - Error Control Coding
Homework 02
"""

p1_header = \
"""
(1) Non-systematic feed-foward encoder on the
    NASA standard (2,1,6) code with generator
    polynomials (117,155).
    Run on the input sequence:
    u = [1 0 1 1 1 0 1 0 0 1 1 0]
"""

p2_header = \
'\n'+'='*50+"""
(2) Systematic recursive encoder (SRCC)."""

p2a_header = \
"""
    (a) SRCC (2,1,1) code with generator
        polynomials (3,1).
        Run on the input sequence:
        u = [1 0 1 1 1 0 1 0 0 1 1 0]"""

p2b_header = \
"""
    (b) SRCC (2,1,4) code with generator
        polynomials (37,21).
        Run on the input sequence:
        u = [1 0 1 1 1 0 1 0 0 1 1 0]
"""

p3_header = \
'\n'+'='*50+"""
(3) Assume a codeword is transmitted over
    the DMC of Problem 12.4. Use the Viterbi
    Algorithm to decode the received sequence:
    r = [[6,7],[6,0],[2,0],[0,5],[6,1],[2,7],[2,1]]
"""

p4_header = \
'\n'+'='*50+"""
(4) A codeword from the code of 12.5 is transmitted
    over a continuous AWGN channel. Use the Viterbi
    Algorithm to decode the received sequence:
    r = [1.72,0.93,2.34,-3.42,-0.14,-2.84,-1.92,0.23,0.78,-0.63,-0.05,2.95,-0.11,-0.55]
"""

def p1():
    print p1_header

    code = (2,1,6)
    gens = (117,155)

    nsff = encoders.SRCC(code,gens)
    msg = np.array((1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0),int)
    v = nsff.encode(msg)

    print 'Codeword:'
    print v.T.flatten()

def p2():
    print p2_header
    msg = [1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0]

    ## (a)
    print p2a_header

    code = (2,1,1)
    gens = (3,1)

    coder = encoders.SRCC(code,gens)
    v = coder.encode(msg)
    print 'v:'
    print v.T.flatten()

    ## (b)
    print p2b_header

    code = (2,1,4)
    gens = (37,21)

    coder = encoders.SRCC(code,gens)
    v = coder.encode(msg)
    print 'v:'
    print v
    print v.T.flatten()

def p3():
    print p3_header

    # Decoding
    code = (2,1,3)
    gens = (13,17)  # ( 1011 1111 )
    #
    # # Q-ary
    c1 = 2.6999
    c2 = 4.28
    parr = np.array([0.434,0.197,0.167,0.111,0.058,0.023,0.008,0.002])
    r = [[6,7],[6,0],[2,0],[0,5],[6,1],[2,7],[2,1]]

    parr2 = parr[::-1]
    p_table = np.vstack((parr,parr2))
    vit = decoders.VirterbiDecoder(code,gens,channel='Qary',metric_params=(p_table,c1,c2))
    v_h = vit.decode(r)
    print 'message:'
    print v_h

def p4():
    print p4_header

    # Decoding
    nkvTuple = (2,1,3)
    gens = (13,17)  # ( 1011 1111 )
    r = [1.72,0.93,2.34,-3.42,-0.14,-2.84,-1.92,0.23,0.78,-0.63,-0.05,2.95,-0.11,-0.55]
    #r = [-0.72, 0.93, -0.34, -3.42, -1.14, -2.84, 2.92, 0.23, -1.78, -0.63, -1.05, 2.95, -2.11, -0.55]

    vit = decoders.VirterbiDecoder(nkvTuple,gens,channel='AWGN',recursive=True)
    v_h = vit.decode(r)
    print 'message:'
    print v_h
    print
    print '='*50

def main():
    #print class_header
    p1()
    p2()
    #p3()
    #p4()

if __name__=="__main__":
    main()