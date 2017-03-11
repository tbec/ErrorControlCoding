from utils import *

class BCJR(object):
    def __init__(self, nkvTuple, GeneratorPolynomials, channel='AWGN', metric_params=None,recursive=False):

        self.terminated = True

        # Internal metrics
        self.n = nkvTuple[0]
        self.k = nkvTuple[1]
        self.v = nkvTuple[2]

        self.recursive = recursive

        # Turn octal generators into binary arrays
        self.g = np.asarray([intStr2BinArray(i) for i in GeneratorPolynomials],int)
        self.g = stripLeadingZeros(self.g)

        # Compute state/output tables
        self.m = self.g.shape[1]-1
        self.v = self.m

        # Create state table and output table
        self.stateTable,self.outputTable = poly2trellis(nkvTuple,GeneratorPolynomials,recursive=recursive)

        # Helper table - gives parents and inputs for states (2 x 2 x states)
        #   where each row is [parent,input]
        self.parentTable = parentStates(self.stateTable)

    def createGammaMat(self,La,Es_N0=1):
        """

        :param La: a priori information
        :param Es_N0:
        :return:
        """

        # --- Some initial Setup ---- #

        # Lc scalar
        Lc_d2 = 0.5 * 4*Es_N0

        # If La wasn't given, assume it's all zeros ( no a priori info )
        if La is None:
            La = np.zeros((self.K,))

        # make the matrix (0 input in front, 1 input in back)
        self.gamma = np.zeros((2**self.v,self.K,2))

        rows,cols,_ = self.gamma.shape

        # ---- gamma for message bits ---- #
        for l in xrange(self.K):
            # Extract the current r sequence
            r_seq = self.r[l*self.n:l*self.n+(self.n)]

            for state in xrange(rows):

                # --- Input 0 --- #
                out = self.outputTable[0,state]

                # Convert output to binary array with 0 --> -1
                outBits = int2EsArray(out,self.n)

                # Calculate gamma according to (12.128a) ( pg. 568 )
                self.gamma[state,l,0] = (-1*La[l]*0.5) + Lc_d2*np.dot(r_seq,outBits)

                # --- Input 1 --- #
                out = self.outputTable[1,state]

                # Convert output to binary array with 0 --> -1
                outBits = int2EsArray(out,self.n)

                # Calculate gamma according to (12.128a) ( pg. 568 )
                self.gamma[state,l,1] = (La[l]*0.5) + Lc_d2*np.dot(r_seq,outBits)


    def createAlphaMat(self):
        """
        Alpha at state s and trellis section l is:
            alpha_l^s = max*(gammas+alphas)
        :return:
        """

        # make the matrix
        self.alpha = np.full((2**self.v,self.K),-50,dtype=np.float64)

        # seed alpha_0(0) = 0
        self.alpha[0,0] = 0

        # and filler up
        rows,cols = self.alpha.shape
        for l in xrange(1,cols):    # For alpha start at l=1, not l=0
            for state in xrange(rows):

                # 'metrics' holds the two branches to put into max*
                metrics = []
                # Get states that go into current state
                parents = self.parentTable[:,:,state]

                for i in xrange(2):  # Two parent states for k = 1

                    parent,ipt = parents[1,i],parents[0,i]

                    # handles every state in trellis
                    #   either puts partial path (gamma_l-1 + alpha_l-1) or NaN in 'metrics'
                    metrics.append(self.gamma[parent,l-1,ipt]+self.alpha[parent,l-1])

                # Perform ln max on values to take weighted max
                self.alpha[state,l] = maxs(metrics)

        return self.alpha

    def createBetaMat(self):
        """
        Alpha at state s and trellis section l is:
            alpha_l^s = max*(gammas+alphas)
        :return:
        """

        # make the matrix
        if self.terminated:
            beta_initial = -50
        else: beta_initial = 0
        self.beta = np.full((2**self.v,self.K),beta_initial,dtype=np.float64)

        # seed beta_h(0) = 0
        self.beta[0,-1] = 0

        # and filler up
        rows,cols = self.beta.shape
        for l in xrange(cols-2,-1,-1):    # For beta start at l=h-1, not h
            for state in xrange(rows):

                beta_max = []
                for ipt in xrange(2):  # Loop through the 2 inputs

                    child_state = self.stateTable[ipt,state]
                    gamma = self.gamma[state,l+1,ipt]
                    beta_ll = self.beta[child_state,l+1]

                    beta_max.append(gamma+beta_ll)

                # Perform ln max on values to take weighted max
                self.beta[state,l] = maxs(beta_max)

        self.beta

    def decode(self, received_sequence, La=None, Es_N0=1):
        """

        :param received_sequence: AWGN received sequence
        :param La: a priori information on received sequence
        :param Es_N0: Symbol to noise ratio (linear)
        :return: LLR and msg
        """

        # Total bits
        self.K = len(received_sequence)/self.n

        # Message bits
        self.h = self.K - self.m

        msg = np.zeros((self.K,))

        self.r = np.asarray(received_sequence)

        # ---- Create gamma matrix ---- #
        self.createGammaMat(La,Es_N0)

        # ---- Create alpha matrix ---- #
        self.createAlphaMat()

        # ---- Create beta matrix ---- #
        self.createBetaMat()

        # ---- Decode received message ---- #
        #   max*( +1 transitions ) - max*( -1 transitions )
        for l in xrange(0,self.K):    # Loop through trellis sections

            # Positive and negative transition lists
            #   alpha_l + gamma_l + beta_l+1
            pos_trans,neg_trans = [],[]

            # Loop through states
            for state in xrange(2**self.v):

                alpha = self.alpha[state,l]

                # Calculate positive transitions
                ipt = 1
                gamma = self.gamma[state,l,ipt]
                next_state = self.stateTable[ipt,state]
                beta = self.beta[next_state,l]

                pos_trans.append(alpha+beta+gamma)

                # Calculate negative transitions
                ipt = 0
                gamma = self.gamma[state,l,ipt]
                next_state = self.stateTable[ipt,state]
                beta = self.beta[next_state,l]

                neg_trans.append(alpha+beta+gamma)

            # Compute L(u_l)
            #   max*( +1 transitions ) - max*( -1 transitions )
            msg_l = maxs(pos_trans)-maxs(neg_trans)
            msg[l] = msg_l

            self.last_decoded_msg_energy = msg
            self.last_decoded_msg = toBits(msg)

        # print 'alpha:'
        # print self.alpha
        # print 'beta:'
        # print self.beta

        return self.last_decoded_msg

class VirterbiDecoder(object):

    def __init__(self, nkvTuple, GeneratorPolynomials, channel, metric_params=None,recursive=False):

        # Internal metrics
        self.n = nkvTuple[0]
        self.k = nkvTuple[1]
        self.v = nkvTuple[2]

        self.recursive = recursive

        # Turn octal generators into binary arrays
        self.g = np.asarray([intStr2BinArray(i) for i in GeneratorPolynomials],int)
        self.g = stripLeadingZeros(self.g)

        # Compute state/output tables
        self.m = self.g.shape[1]-1
        # self.stateTable,self.outputTable = self.createTrellis()
        self.stateTable,self.outputTable = poly2trellis(nkvTuple,GeneratorPolynomials,recursive=recursive)

        # Correct output table for 4a recursive
        # self.outputTable = np.array([[0,2,0,2,2,0,2,0],[1,3,1,3,3,1,3,1]],int)

        # Set channel type ('Binary','Qary','AWGN')
        self.int_table = None
        self.change_channel(channel,metric_params)

    def _prob2int_table(self,metric_params):
        metric_table = metric_params[0] # probability table
        c1 = metric_params[1]
        c2 = metric_params[2]
        return self._Qary(metric_table,c1,c2)

    def _Qary(self,p,c1,c2):
        return c2*(np.log10(p)+c1)

    def _newList(self):
        """
        :return: None type list for states
        """
        return [None]*2**self.v

    # def createTrellis(self):
    #     """
    #     Create the trellis structure for decoding.
    #     :return:
    #     """
    #
    #     numInputSymbols = 2**self.k
    #     numStates = 2**self.v
    #
    #     # State/Output Table Initialization
    #     nextStates = np.zeros((numInputSymbols,numStates),int)
    #     outputs = np.zeros((numInputSymbols,numStates),int)
    #
    #     # Loop through states
    #     for state in xrange(numStates):
    #
    #         bitState = int2bitArray(state,self.v)
    #
    #         # Loop through inputs (0,1) for k=1
    #         for ipt in xrange(numInputSymbols):
    #
    #             # Compute next state
    #             nextStates[ipt,state] = (state*2+ipt) % 2**self.v
    #
    #             # Compute output
    #             bitOutput = [(((ipt&self.g[i,0])+(bitState&self.g[i,1:]).sum())%2) for i in xrange(self.n)]
    #             outputs[ipt,state] = bitArray2int(np.asarray(bitOutput,int))
    #
    #     return nextStates,outputs

    def _branch_metric_Hamming(self,r_seq):

        bm_t = self.outputTable.copy()

        for yx in np.ndindex(self.outputTable.shape):
            bm_t[yx] = hamming_int(bm_t[yx],r_seq)

        return bm_t

    def _branch_metric_Qary(self,r_seq):

        bm_t = self.outputTable.copy()

        for yx in np.ndindex(bm_t.shape):
            bm = 0
            out = int2bitArray(self.outputTable[yx],max(r_seq.shape))
            for i in xrange(max(r_seq.shape)):
                bm += self.int_table[int(out[i]),int(r_seq[i])]
            bm_t[yx] = bm

        return bm_t

    def _branch_metric_AWGN(self,r_seq):

        bm_t = self.outputTable.copy().astype(np.float)
        for yx in np.ndindex(bm_t.shape):
            bm = 0
            out = int2bitArray(self.outputTable[yx],max(r_seq.shape))
            for i in xrange(max(r_seq.shape)):
                if out[i]==0:
                    bm -= r_seq[i]
                else: bm += r_seq[i]
            bm_t[yx] = bm

        return bm_t

    def change_channel(self,channel,metric_params=None):
        self.channel = channel
        if channel=='Qary':
            self.int_table = self._prob2int_table(metric_params).astype(int)

    def decode(self,coded_msg):
        """
        Decode input message
        :param coded_msg: received message
        :return: decoded msg
        """

        # Get length of the desired message
        if self.channel=='Qary':
            len_r = len(coded_msg)      # code is sent as list of lists
        elif self.channel=='AWGN':
            len_r = len(coded_msg)/self.n   # code is sent as single list, must be divided into v's

        coded_msg = np.asarray(coded_msg,np.float)

        h = len_r - self.v

        # Make the path metric table (states x time units x (previous state,path metric))
        pmTable = np.full((2**self.v,len_r+1,2),np.nan,dtype=np.float64)

        # Seed with zero-state  0 weight / no parent
        pmTable[0,0] = np.array([0,np.nan])

        # Some checking
        if self.channel!='Qary' and self.channel!='AWGN':
            print"""Channel type not supported. Please select a valid channel type:
                    'Qary' | 'AWGN'"""
            sys.exit(0)

# ------ Begin Loop -------------------------

        # Loop through trellis sections
        # for t in xrange(pmTable.shape[1]-1):
        for t in xrange(h):

            # Set trellis step r, and create branch metric table

            if self.channel == 'Qary':
                r_seq = coded_msg[t]
                bm_t = self._branch_metric_Qary(r_seq)

            elif self.channel == 'AWGN':
                r_seq = coded_msg[2*t:2*t+2]
                bm_t = self._branch_metric_AWGN(r_seq)

            inputs,states = bm_t.shape
            for i in xrange(inputs):
                for s in xrange(states):

                    # partial path metric at next state depends on
                    #   path metric at current state + branch metric
                    next_state = self.stateTable[i,s]

                    if not np.isnan(pmTable[s,t,0]):

                        pm = pmTable[s,t,0]+bm_t[i,s]

                        # Haven't added first partial path metric
                        if np.isnan(pmTable[next_state,(t+1),0]):
                            pmTable[next_state,(t+1),0] = pm
                            pmTable[next_state,(t+1),1] = s
                        else:
                            if self.channel=='Qary' or self.channel=='AWGN':
                                if pm>pmTable[next_state,(t+1),0]:
                                    pmTable[next_state,(t+1),0] = pm
                                    pmTable[next_state,(t+1),1] = s

                            else:   # Case for hamming distance
                                pmTable[next_state,(t+1),0] = min(pm,pmTable[next_state,(t+1),0])


# ----- Force trellis to zero state -----------------------------
        for t in xrange(h,pmTable.shape[1]-1):

            if self.channel == 'Qary':
                r_seq = coded_msg[t]
                bm_t = self._branch_metric_Qary(r_seq)

            elif self.channel == 'AWGN':
                r_seq = coded_msg[2*t:2*t+2]
                bm_t = self._branch_metric_AWGN(r_seq)

            # Loop through each state at time t
            for s in xrange(pmTable.shape[0]):

                # Check for NaN at current state
                if not np.isnan(pmTable[s,t,0]):

                    # Next state is a zero shifted into LSB
                    next_state = s*2 % (2**self.v)

                    # find input that forces current state to next state
                    ipt = np.argwhere(self.stateTable[:,s]==next_state)[0,0]

                    # Calculate the path metric for the next state (t+1)
                    pm = pmTable[s,t,0]+bm_t[ipt,s]

                    # Haven't added first partial path metric
                    if np.isnan(pmTable[next_state,(t+1),0]):
                        pmTable[next_state,(t+1),0] = pm
                        pmTable[next_state,(t+1),1] = s
                    else:
                        if self.channel=='Qary' or self.channel=='AWGN':
                            if pm>pmTable[next_state,(t+1),0]:
                                pmTable[next_state,(t+1),0] = pm
                                pmTable[next_state,(t+1),1] = s

                        else:   # Case for hamming distance
                            pmTable[next_state,(t+1),0] = min(pm,pmTable[next_state,(t+1),0])

        # Traceback
        #   Qary and AWGN: trace highest branch metric

        msg = []
        state = 0
        for t in xrange(pmTable.shape[1]-1,0,-1):
        #for t in xrange(4,0,-1):

            # Get the parent state
            prev_state = pmTable[state,t,1]  # each state in table also holds parent state at index 1

            if np.isnan(prev_state): break
            prev_state = int(prev_state)

            # append input that goes from parent state to current state
            ipt = np.argwhere(self.stateTable[:,prev_state]==state)[0,0]

            msg.append(ipt)

            state = prev_state

        msg = msg[::-1]

        print pmTable[:,:,0]
        print pmTable[:,:,1]

        return msg[:h]


