
########################################################################

# Distributions class encapsulates the different distributions that we
# may need to sample structures:
#   ssDomainLengthDist - distribution of ss domain lengths in nm (given length in nt)
#   dsDomainLengthDist - distribution of ds domain lengths in nm (given length in nt)
#   tetherAngleDist - distribution of angles for domains immediately attached to a tether
#   ssDomainAngleDist - distribution of angles between an ss domain and any other domain
#   dsdsDomainAngleDist - distribution of angles between two ds domains (i.e., a nick)
# Each of these distributions will be represented by a Python object with an appropriate method to call
#   sampleAngle() for angle distributions
#   sampleLengthNm(n) for length distributions, where n is length of domain in nt

class Distributions():

    def __init__(self, ssDomainLengthDist, dsDomainLengthDist, tetherAngleDist,
                 ssDomainAngleDist, dsdsDomainAngleDist):
        self.ssDomainLengthDist = ssDomainLengthDist
        self.dsDomainLengthDist = dsDomainLengthDist
        self.tetherAngleDist = tetherAngleDist
        self.ssDomainAngleDist = ssDomainAngleDist
        self.dsdsDomainAngleDist = dsdsDomainAngleDist

########################################################################
