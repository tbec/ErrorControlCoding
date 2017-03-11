from utils import *
"""
Classes of nonsystematic codes
"""


class NSFFCC(object):
    """
    Nonsystematic, feed-forward encoder
    """

    def __init__(self,nkvTuple,generatorPolynomials):
        self.n = nkvTuple[0]
        self.k = nkvTuple[1]
        self.v = nkvTuple[2]
        self.genList = list(generatorPolynomials)

        # Get generators in binary
        self.generators = stripLeadingZeros(map(intStr2BinArray,self.genList))
        self.m = self.generators.shape[1]-1

    def encode(self,msg):
        """
        Messgae can be binary list or numpy array
        :param u: message to encode
        :return: encoded message (v)
        """
        u = np.array(msg,int)

        codes = []
        for g in self.generators:
            codes.append(mod2Convolve(u,g))

        return np.array(codes,int).T.flatten()

class SRCC(object):
    """
    Systematic recursive convolutional code
    """
    def __init__(self, nkvTuple, generatorPolynomials):
        self.n = nkvTuple[0]
        self.k = nkvTuple[1]
        self.v = nkvTuple[2]
        self.g = generatorPolynomials
        self.gBits = np.asarray([intStr2BinArray(i) for i in generatorPolynomials],int)
        self.gBits = stripLeadingZeros(self.gBits)

        # Compute trellis structure
        self.state_table,self.output_table = poly2trellis(nkvTuple,self.g,recursive=True)

    def encode(self,msg):

        h = len(msg)
        u = np.asarray(msg,int)

        # Create array for encoded message
        encoded_msg = np.zeros((self.n,self.v+h),int)

        # Trellis starts at zero state
        state = 0
        ipt = 0
        # Loop through message bits to create encoded message
        for t in xrange(h):

            ipt = u[t]

            # compute v1-vn at time t
            encoded_msg[:,t] = int2bitArray(self.output_table[ipt,state],self.n)

            # Set current state equal to next state (increment time)
            state = self.state_table[ipt,state]

        # Force trellis to zero state
        for t in xrange(h,h+self.v):

            # find input to force a zero (feedback + input = 0) into LSB

            # Force a zero into next state LSB
            next_state = (state*2) % (2**self.v)

            # Find input corresponding to state -> next state transition
            ipt = np.argwhere(self.state_table[:,state] == next_state)[0,0]

            encoded_msg[:,t] = int2bitArray(self.output_table[ipt,state],self.n)

            state = next_state
        #
        #     # Find input (equal to feedback):
        #     ipt = (self.gBits[0,1:]&bit_state).sum() % 2
        #
        #     v1 = (bit_state&self.gBits[1,1:]).sum() % 2
        #
        #     encoded_msg[0,i] = ipt
        #     encoded_msg[1,i] = v1
        #
        #     msg_in.append(ipt)
        #
        #     b = (bitArray2int(bit_state)<<1)%2**self.v
        #
        #     bit_state = int2bitArray(b,self.v)
        #
        #     print 'next: ',bitArray2int(bit_state)
        #     print


        return encoded_msg.T.flatten()


