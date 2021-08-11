import cProfile
import queue

import numpy as np
from numba import float64, jit, njit
from numba.experimental import jitclass
from numpy.polynomial.polynomial import polyfit


@njit
def distance(a, b):
    return np.sqrt(np.dot(a-b, a-b))


def radiusOfGyration(aggregate, monomerRadius):
    # requires array of size (n,3)
    if len(aggregate) == 1:
        return np.sqrt(3/5) * monomerRadius

    center = np.mean(aggregate, axis=0)
    offsets = aggregate - center

    # computes dot product along each coordinate (x,y,z) and sums them [equivalent to x^2 + y^2 + z^2]
    # normalizes by number of monomers in aggregate
    # computes square root to get final radius of gyration
    return np.sqrt(np.square(offsets).sum() / len(aggregate))


class FractalGraph():
    # don't forget to type annotate everything

    def __init__(self, filepath: str = None):
        with open(filepath, 'r') as file:
            filecontents = file.read().split('\n')[:-1]

            self.monomerRadius = np.double(filecontents[0].split()[0])

            # temp = [line.split()[1:] for line in filecontents]
            # self.coordinates = np.array([[float(i) for i in line] for line in temp])
            self.coordinates = np.array([[float(i) for i in line.split()[1:]] for line in filecontents])
            self.aggregateSize = len(self.coordinates)

            self.graphNeighbors = []
            for i in range(self.aggregateSize):
                # there is a functional/range-based way to do this, i think
                neighbors = [
                    a for a in range(self.aggregateSize)
                    if distance(self.coordinates[i], self.coordinates[a]) < 2.01*self.monomerRadius
                    and a != i
                ]
                self.graphNeighbors.append(neighbors)

        self.dimension = None
        self.prefactor = None

    def fractalBFS(self, nodeID, parentID, rootID, subaggregate, subaggregateIDs, boxDims):
        for neighborID in self.graphNeighbors[nodeID]:
            if neighborID == parentID or neighborID in subaggregateIDs:
                continue
            elif not np.any((np.abs(self.coordinates[neighborID] - self.coordinates[rootID]) - boxDims) > 0):
                subaggregate.append(self.coordinates[neighborID])
                subaggregateIDs.append(neighborID)

                self.fractalBFS(neighborID, nodeID, rootID, subaggregate, subaggregateIDs, boxDims)

    def getFractalParameters(self, boxLow, boxHigh, iterations=5000):
        if self.dimension is not None and self.prefactor is not None:
            return self.dimension, self.prefactor

        nMonomerArray = []
        radiusGyrationArray = []
        for nodeID in np.random.randint(0, self.aggregateSize, size=iterations):
            boxDims = np.random.randint(low=boxLow, high=boxHigh, size=3)

            subaggregate = [self.coordinates[nodeID]]
            subaggregateIDs = [nodeID]

            self.fractalBFS(nodeID, -1, nodeID, subaggregate, subaggregateIDs, boxDims)

            nMonomerArray.append(len(subaggregateIDs))
            radiusGyrationArray.append(radiusOfGyration(np.array(subaggregate), self.monomerRadius))

        radiusGyrationArray = np.array(radiusGyrationArray)
        nMonomerArray = np.array(nMonomerArray)
        # print(radiusGyrationArray)
        # print(nMonomerArray)

        logNs = np.log10(nMonomerArray)
        logrg = np.log10(radiusGyrationArray) - np.log10(self.monomerRadius)

        regression = polyfit(logrg, logNs, 1)

        self.dimension = regression[1]
        self.prefactor = np.power(10, regression[0])

        return self.dimension, self.prefactor
