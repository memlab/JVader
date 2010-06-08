from numpy        import any, argmax, array, cov, diag, dot, exp, log, mean, max, pi, sqrt, sum, vstack, zeros
from numpy.linalg import det, inv, LinAlgError
from numpy.random import uniform, normal
from scipy.linalg.matfuncs import sqrtm

def mvn_pdf(x, mu, sigma):
    """
    Computes the value of a multi-variate normal with
    mean mu and covariance matrix sigma at point x.
    """
    return 1./( (2*pi)**(len(x)/2.) * sqrt(det(sigma))) * \
            exp(-1./2 * dot(dot((x-mu).T, inv(sigma)), x-mu))

class GausMixHMM:
    """
    A hidden markov model, observations at each state
    are modeled by a mixture of gaussians.
    """
    def __init__(self):
        """
        Creates a GausMixHMM object.
        """
        # The number of states.
        self.numStates = 0

        # The number of mixtures per state.
        # We keep the number the same for each state
        # not only because this is easier to do, but
        # empirical results have shown this is better
        # (as opposed to different mixture numbers per
        #  state) for isolated word recognition.
        self.numMixtures = 0

        # Vector of probabilities for being in each 
        # state at time 0.
        self.pi = None

        # Matrix of state transition probabilities.
        # i.e. a[i][j] is the probability of being
        # in state i and going to state j.
        self.a = None

        # Matrix of mixture coefficients. Each row
        # is a mixture and each column is a state.
        self.c = None

        # A 3D array containing mixture means. The
        # first dimension is mixtures and the second
        # dimension is states.
        # i.e. mu[m][i] is the mean vector for the
        # mth mixture of state i.
        self.mu = None

        # A 4D array containing mixture covariance 
        # matrices. The first dimension is mixtures 
        # and the second dimension is states.
        # i.e. sigma[m][i] is the covariance matrix
        # for the mth mixture of state i.
        self.sigma = None

        # If set to True, outputs various status 
        # updates, especially during training.
        self.verbose = False

        # Is this a garbage model?
        self.garbage = False

        # The classification label for observations 
        # that this model is fit to. 
        self.label = None

        # Whether this is a left-right HMM.
        self.leftRight = False

    def _forward(self, observation, b):
        """
        Compute the forward variable for all state
        and time combinations.
        
        alpha[t][i] = P(o1,o2,o3,..,ot,qt = i | model)
        The probability of observing o1..ot and being
        in state i in the end, given the current model
        parameters. 

        A scale factor is used to prevent underflow. 
        The function returns a 2-tuple, the first component 
        is alpha, and the second component is the scale
        factor for each time point.

        observation should be a matrix, with each row 
        constituting a different time point.

        b should be a matrix of precomputed output 
        probabilities for each time point and state. 
        i.e. b = sum(self._outputProb(observation), 2)
        We pass this in from the outside for efficiency 
        sake.. it's usually already computed anyway, so 
        there's no sense in recomputing it.
        """
        scale     = zeros((observation.shape[0]))
        alpha     = zeros((observation.shape[0], self.numStates))
        alpha[0]  = self.pi * b[0]
        scale[0]  = sum(alpha[0])
        alpha[0] /= scale[0]
        for t in range(1, observation.shape[0]):
            alpha[t]  = sum(self.a.T * alpha[t-1], 1) * b[t]
            scale[t]  = sum(alpha[t])
            alpha[t] /= scale[t]

        return (alpha, scale)

    def _outputProb(self, observation):
        """
        Computes for each time step of the given 
        observation the probability of being in each 
        possible state under each mixture using the 
        current model values.

        Returns a 3D array where 
        b[t][s][m] = P(observing Ot in state s under mixture m | observation)
        """
        b = zeros((observation.shape[0], self.numStates, self.numMixtures))
        for t in range(observation.shape[0]):
            for s in range(self.numStates):
                for m in range(self.numMixtures):
                    b[t][s][m] = self.c[m][s] * mvn_pdf(observation[t], self.mu[m][s], self.sigma[m][s])
        return b 

    def _makeLeftRight(self):
        """
        Turns the current transition matrix into a
        left-righttransition matrix. That is,
        from each state, you may only jump to a higher
        state; you cannot transition backwards.
        """
        for i in range(self.a.shape[0]):
            self.a[i][0:i] = 0

    def baumWelch(self, observations, likTolerance=0.0001, maxCycles=200):
        """
        Fits the model parameters to maximize the 
        likelihood of the given data.

        observations should be a list of matrices. 
        The rows in each observation matrix are fixed 
        time steps of the observation.
        """
        # Minimum value for various parameters, including mixture 
        # coefficients and cov matrix components.
        minValue = 0.0001

        # All observations stacked one after the other in the
        # time dimension, for easy access later.
        observationsAll = zeros((0, observations[0].shape[1]))
        for observation in observations:
            observationsAll = vstack((observationsAll, observation))

        # Log likelihood of the data using the given parameters,
        # unnormalized.
        logLikelihood     = 0
        oldLogLikelihood  = 0
        logLikelihoodBase = 0

        for cycle in range(maxCycles):
            logLikelihood = 0

            # List of gamma values for each observation, where
            # gamma[t][i] = P(in state i at time t | observation)
            # Each observation is stacked one after the other in 
            # the time dimension.
            gammaAll  = zeros((0, self.numStates))

            # Sum of gamma[0] over all observations.
            gamma0Sum = zeros((1, self.numStates))

            # List of xi values for each observation, where
            # xi[t][i][j] = P(in state i at time t, in state j at t+1 | observation)
            # Each observation is stacked one after the other in the time
            # dimension.
            xiAll     = zeros((0, self.numStates, self.numStates))

            for observation in observations:
                bPerMixture = self._outputProb(observation)
                b           = sum(bPerMixture, 2)
                
                (alpha, scale) = self._forward(observation, b)  
                logLikelihood += sum(log(scale))
                #if self.verbose:
                #    print "Baum-Welch, LL = ", str(logLikelihood)

                # Computes the backward variable for all
                # state and time combinations.
                # 
                # beta[t][i] = P(ot+1, ot+2,..,oT | qt = i,model)
                # The probability of observing ot+1..the end,
                # given that we're currently in state i.
                beta = zeros((observation.shape[0], self.numStates))

                beta[observation.shape[0]-1][:] = 1
                for t in range(observation.shape[0]-2, -1, -1):
                    beta[t]  = sum(self.a * b[t+1] * beta[t+1], 1)
                    beta[t] /= scale[t]

                # Calculate gamma for this observation.
                gamma      = alpha*beta
                gamma      = (gamma.T / sum(gamma,1)).T
                gamma0Sum += gamma[0]
                gammaAll   = vstack((gammaAll, gamma))

                # Calculate xi for this observation.
                xi = zeros((observation.shape[0], self.numStates, self.numStates))
                for t in range(observation.shape[0]-1):
                    xi[t]  = (self.a.T * alpha[t]).T * b[t+1] * beta[t+1]
                    xi[t] /= sum(xi[t])
                xiAll = vstack((xiAll, xi))
           
            # If we've reached a small neighberhood of some minimum, stop.
            if cycle <= 2:
                logLikelihoodBase = logLikelihood
            else:
                likeDiff    = logLikelihood-logLikelihoodBase
                oldLikeDiff = (1+likTolerance)*(oldLogLikelihood-logLikelihoodBase)
                if self.verbose:
                    print "Baum-Welch, difference = ", str(likeDiff-oldLikeDiff)
                if likeDiff < oldLikeDiff:
                    break
            oldLogLikelihood = logLikelihood

            if self.verbose:
                print "Baum-Welch, cycle = ", str(cycle)

            # These two may need to be scaled somehow, as they sometimes
            # result in underflow. Since we're only using one mixture
            # per state now, it doesn't matter.
            bPerMixture = self._outputProb(observationsAll)
            b           = sum(bPerMixture, 2)

            # Re-estimate priors using updated forward-backward 
            # values. These are the sum of the expected 
            # frequencies of being in state i at time 0 (gamma[0]'s), 
            # normalized over all observations.
            #
            # We now use a left-right HMM instead and we always
            # start in the first state. If someone wants a separate,
            # general HMM implementation later, you can use the below 
            # line to re-estimate pi.
            #self.pi = gamma0Sum/float(len(observations))

            # Re-estimate state transition matrix. Each transition
            # probability is the ratio between the expected 
            # number of transitions from i to j (xi[:][i][j]) and the 
            # expected number of times at state i (gamma[:][i]).
            # Each term is taken over every time point except the 
            # last in each observation (there's no transition from the
            # last time point to anywhere). We recalculate gamma based 
            # on xi since our running gamma is over ALL time points and 
            # is used elsewhere below.
            self.a = (sum(xiAll, 0).T / sum(sum(xiAll,2),0)).T

            # Recompute mixtures.
            for m in range(self.numMixtures):
                # Gamma for this mixture only.
                gammaMix    = gammaAll * (bPerMixture[:,:,m] / b)
                gammaMixSum = sum(gammaMix, 0)

                # Ratio of the number of times in a state and using this
                # mixture to the total number of times in the state.
                self.c[m]   = gammaMixSum / sum(gammaAll, 0)
                self.c[m][self.c[m]<minValue] = minValue

                # The observations weighed by the mixture ratio.
                self.mu[m]  = (dot(gammaMix.T, observationsAll).T / gammaMixSum).T
          
                # Variance weighed by the mixture ratio.
                for s in range(self.numStates):
                    tmp = observationsAll - self.mu[m][s]
                    self.sigma[m][s] = dot(tmp.T * gammaMix[:,s], tmp) / gammaMixSum[s]
                    self.sigma[m][s][self.sigma[m][s]<minValue] = minValue

            # If one of the mixtures has a singular covariance
            # matrix, restart.
            for m in range(self.numMixtures):
                for s in range(self.numStates):
                    if det(self.sigma[m][s]) == 0:
                        if self.verbose:
                            print "Baum-Welch, degenerate mixture detected, restarting."
                        self.initParams(self.numStates, self.numMixtures, observations, self.isLeftRight())
                        self.baumWelch(observations, likTolerance, maxCycles)
                        return

    def distance(self, otherHmm, observation):
        """
        Computes the distance between this model and the given
        other model using the observation. The metric used is
        1/t * [log(Ot | model1) - log(Ot | model2)]
        """
        return 1.0/observation.shape[0] * \
                (self.logLikelihood(observation) - otherHmm.logLikelihood(observation))
   
    def initParams(self, numStates, numMixtures, observations, leftRight=True):
        """
        Initialize parameters with default values,
        usually done for a new model before training.
        
        observations should be the same structure later given
        to baumWelch(..). numStates it the number of states, and
        numMixtures is the number of mixtures per state (we fix
        this to be the same number for all states).
        """
        # All observations stacked one after the other in the
        # time dimension.
        observationsAll = zeros((0, observations[0].shape[1]))
        for observation in observations:
            observationsAll = vstack((observationsAll, observation))

        self.leftRight    = leftRight
        self.numStates    = numStates
        self.numMixtures  = numMixtures
        self.pi           = zeros((self.numStates))
        self.pi[0]        = 1
        self.a            = uniform(0, 1, (self.numStates, self.numStates))
        if leftRight:
            self._makeLeftRight()
        self.a            = (self.a.T/sum(self.a,1)).T # Each row adds to 1.
        self.c            = uniform(0, 1, (self.numMixtures, self.numStates))
        self.c           /= sum(self.c, 0)
        self.sigma        = diag(diag(cov(observationsAll, rowvar=0)))
        self.sigma.shape  = (1, 1, self.sigma.shape[0], self.sigma.shape[1])
        self.sigma        = self.sigma.repeat(self.numMixtures, 0).repeat(self.numStates, 1)

        dev = sqrtm(self.sigma[0][0])
        self.mu = zeros((numMixtures, numStates, observationsAll.shape[1]))
        for m in range(numMixtures):
            self.mu[m] = mean(observationsAll, 0) + dot(normal(0, 1, (numStates, observationsAll.shape[1])), dev)

    def logLikelihood(self, observation):
        """
        Returns the log likelihood of the given observation, 
        'observation' is a matrix where each row is a fixed 
        time step.
        """
        bPerMixture = self._outputProb(observation)
        b           = sum(bPerMixture, 2)

        (alpha, scale) = self._forward(observation, b)

        return sum(log(scale))

    def viterbi(self, observation):
        """
        Returns a 2-tuple. The first component is a vector 
        containing the most likely state to be in at each time 
        step of the observation, and the second component is
        the log likelihood of this (optimal) path.

        'observation' is a matrix where each row is a fixed 
        time step.
        """
        pi          = log(self.pi)
        a           = log(self.a)
        bPerMixture = self._outputProb(observation)
        b           = log(sum(bPerMixture, 2))

        # delta[t][i] = max P(q1, q2, ..., qt-1, qt = i, o1, o2, ..., ot | model)
        # the maximum probability along a single path at time t, ending at state i.
        delta = zeros((observation.shape[0], self.numStates))

        # psi[t][j] = the state, i, that maximizes delta[t-1][i]*a[i][j]. 
        # i.e. the best possible state prior to j in the sequence.
        psi = zeros((observation.shape[0], self.numStates), dtype=int)

        # Compute delta and psi.
        delta[0]  = pi + b[0]
        psi[0][:] = -1
        for t in range(1, observation.shape[0]):
            delta[t] = max(delta[t-1] + a.T, 1) + b[t]
            psi[t]   = argmax(delta[t-1] + a.T, 1)

        # Backtrack to find best state sequence.
        bestStateSeq = zeros((observation.shape[0]), dtype=int)
        bestStateSeq[observation.shape[0]-1] = argmax(delta[observation.shape[0]-1])
        for t in range(observation.shape[0]-2, -1, -1):
            bestStateSeq[t] = psi[t+1][bestStateSeq[t+1]]

        return (bestStateSeq, max(delta[-1]))

    def isGarbage(self):
        """
        Identifies whether this is a garbage model.
        """
        return self.garbage

    def isLeftRight(self):
        """
        Returns whether this is a left-right HMM.
        """
        return self.leftRight

    def getLabel(self):
        """
        Returns the classification label for observations 
        that this model was fit to.
        """
        return self.label

    def setLabel(self, newLabel):
        """
        Sets the classification label for observations 
        that this model is fit to.
        """
        self.label = newLabel

    def setGarbage(self, newGarbage):
        """
        Identifies whether this is a garbage model.
        """
        self.garbage = newGarbage

    def setVerbose(self, newVerbose):
        """
        If set to True, outputs various status updates,
        especially during training.
        """
        self.verbose = newVerbose

