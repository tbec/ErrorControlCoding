# from helpers import *
#
# """
# Classes of nonsystematic codes
# """
#
#
# class NSFFCC:
#     """
#     Nonsystematic, feed-forward encoder
#     """
#     # TODO: k > 1 functionality
#
#     def __init__(self,nkvTup):
#         self.n = nkvTup[0]
#         self.k = nkvTup[1]
#         self.v = nkvTup[2]
#         self.m = None           # memory order
#         self.h = None           # Length of message
#         self.u = None           # Message (User input)
#         self.codeword = None    # Codeword from message (u)
#         self.genList = None     # List of generators
#         self.generators = None  # Generators in binary
#         self.codes = None
#
#     def setGenerators(self,genTup):
#         """
#         Generators must be octal
#         :param genTup: tuple of generators
#         :return: N/A
#         """
#         self.genList = list(genTup)
#
#         # Get generators in binary
#         self.generators = stripLeadingZeros(map(intStr2BinArray,self.genList))
#         self.m = self.generators.shape[1]-1
#
#     def generateCodeword(self,u):
#         """
#         Messgae can be binary list or numpy array
#         :param u: message to encode
#         :return: encoded message (v)
#         """
#         self.u = np.array(u,int)
#
#         self.codes = []
#         for g in self.generators:
#             self.codes.append(mod2Convolve(self.u,g))
#
#         self.codeword = np.array(self.codes).T.flatten()
#
#
#     def separateCodewords(self):
#         self.codes = parseCodeword(self.codeword,self.n)
#         return self.codes
#
# class SRCC:
#     """
#     Systematic recursive convolutional code
#     """
#     def __init__(self, nkvTup):
#         self.n = nkvTup[0]
#         self.k = nkvTup[1]
#         self.v = nkvTup[2]
#
# def poly2trellis(ConstraintLength,CodeGenerator):
#     return
#
# def srcc(msg,genPolys):
#
#     msg = np.asarray(msg,int)
#
#     # Serial input initially
#     k = 1
#     n = len(genPolys)
#
#     # Initialize feedback
#     fb = 0
#
#     # Convert octals to binary arrays
#     genArr = np.asarray([intStr2BinArray(i) for i in genPolys],int)
#     genArr = stripLeadingZeros(genArr)
#     genBool = genArr.astype(bool)
#
#     # Get generator memory order
#     m = genArr.shape[1]-1
#
#     # Initialize codeword matrix
#     v = np.zeros((n,len(msg)+m+1),int)
#
#     # Initialize shift register (include message bit)
#     reg = np.zeros((1,m+1),int)
#
#     # Loop over the message
#     for imsg in xrange(len(msg)):
#
#         # Add in message bit
#         reg[0,0] = msg[imsg]
#
#         # Systematic, so v[0,i] is same as msg
#         v[0,imsg] = reg[0,0]
#
#         print imsg,'\t',reg
#
#         # Loop over the rest of the codes
#         for code in xrange(1,n):
#             g = genBool[code,:].reshape(reg.shape)
#             v[code,imsg] = reg[g].sum()%2
#
#         # First generator g(0)
#         g = genBool[0,:].reshape(reg.shape)
#         fb = reg[g].sum()%2
#
#         # Perform shift
#         reg = np.roll(reg,1)
#
#         reg[0,1] = fb
#
#     # # Send trellis to zero state
#     print '-'*40
#     for ind in xrange(len(msg),len(msg)+m+1):
#
#         # Next input needs to equal feedback
#         # g = genBool[0,:].reshape(reg.shape)
#         # fb = reg[g].sum() % 2
#         #g = genBool[0,1:]
#         reg[0,0] = reg[0,1:].sum() % 2
#
#         v[0,ind] = reg[0,0]
#
#         # Loop over the rest of the codes
#         for code in xrange(1,n):
#             g = genBool[code,:].reshape(reg.shape)
#             v[code,ind] = (reg[0,-1]+reg[0,0])%2
#
#         # First generator g(0)
#         g = genBool[0,:].reshape(reg.shape)
#         fb = reg[g].sum() % 2
#
#         print ind,'\t',reg
#
#         # Perform shift
#         reg = np.roll(reg,1)
#
#         reg[0,1] = fb
#
#     print '-'*40
#     print v
#     print v.T.flatten()
#
#
# # --------------------------------------------------------------------
# # --------------------------------------------------------------------
#
# msg = [1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0]
# srcc(msg,(37,21))
#
#
