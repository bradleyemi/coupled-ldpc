### New implementation of LDPC and Spatially-Coupled LDPC for EE 376A
### Passes codewords through BEC (instead of BSC)

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy, copy
from time import time
import cPickle as pickle

def loadLDPC(name):
    """
    :param - name: the name of the file containing LDPC matrices
  
    return values:
    G: generator matrix
    H: parity check matrix
    """
    A = sio.loadmat(name)
    G = A['G']
    H = A['H']
    return G, H

#G, H = loadLDPC("ldpc36-1600.mat")
print "Loading H"
s = time()
H = np.load("coupled_h_large.npy")
print "Loaded, took {} s".format(time()-s)

class Factor(object):
    def __init__(self, node_type, id, neighbors):
        self.node_type = node_type
        self.id = id
        self.neighbors = neighbors

# A class to hold the factor graph structure of H and initialize messages
class Graph(object):
    def __init__(self, H):
        self.H = H

        self.n_vars = H.shape[1]
        self.n_checks = H.shape[0]

        self.vars = []
        self.checks = []

        for i in range(self.n_vars):
            print "Constructing variable factor", i
            neighbors = np.nonzero(H[:,i])
            self.vars.append(Factor("var", i, neighbors[0]))

        for i in range(self.n_checks):
            print "Constructing check factor", i
            neighbors = np.nonzero(H[i,:])
            self.checks.append(Factor("check", i, neighbors[0]))
        

class ClusterGraph(object):
    def __init__(self, graph, codeword):
        self.graph = graph
        self.orig_codeword = codeword
        self.decoded_codeword = copy(codeword)

        self.var_to_check_msgs = {}
        self.check_to_var_msgs = {}

        for var in self.graph.vars:
            for neighbor in var.neighbors:
                self.var_to_check_msgs[(var.id, neighbor)] = -1 #Message("var_to_check", -1)

        for check in self.graph.checks:
            for neighbor in check.neighbors:
                self.check_to_var_msgs[(check.id, neighbor)] = -1 #Message("check_to_var", -1)

    def decode(self, iterations, yield_errors=False):
        errors = len([c for c in self.decoded_codeword if c == -1])
        print "decoding for error probability", float(errors) / self.graph.n_vars
        for iteration in xrange(iterations):
            # First iterate over var -> check messages
            for source_dest, msg in self.var_to_check_msgs.iteritems():
                source = source_dest[0]
                dest = source_dest[1]
                source_neighbors = self.graph.vars[source].neighbors
                in_msgs = [self.check_to_var_msgs[(neigh, source)] for neigh in source_neighbors if neigh != dest]
                # If any incoming message from a check is not -1, we know the variable must be equal to the message value
                # so we can pass it forward and decode the variable bit
                for in_msg in in_msgs:
                    if in_msg != -1:
                        self.var_to_check_msgs[source_dest] = in_msg
                        self.decoded_codeword[source] = msg
                        break
                if self.decoded_codeword[source] != -1:
                    self.var_to_check_msgs[source_dest] = self.decoded_codeword[source]
            
            # Now iterate over check -> var messages
            for source_dest, msg in self.check_to_var_msgs.iteritems():
                source = source_dest[0]
                dest = source_dest[1]
                source_neighbors = self.graph.checks[source].neighbors
                in_msgs = [self.var_to_check_msgs[(neigh, source)] for neigh in source_neighbors if neigh != dest]
                # A check must have all incoming messages non-erased to be certain of the value of the destination
                # so -1 cannot be in the incoming messages. 
                # If all incoming messages are known, pass value that makes parity check = 0 (mod 2) and decode bit
                if -1 not in in_msgs:
                    self.check_to_var_msgs[source_dest] = sum(in_msgs) % 2
                    self.decoded_codeword[dest] = sum(in_msgs) % 2
            new_errors = len([c for c in self.decoded_codeword if c == -1])
            if yield_errors:
                yield self.decoded_codeword
            if (new_errors == errors) or (new_errors == 0):
                break
            else:
                errors = new_errors
        #return float(errors) / self.graph.H.shape[1]

def transmit(codeword, epsilon):
    # Take the convention that -1 is an erasure bit
    new_codeword = copy(codeword)
    for i, bit in new_codeword:
        if (i % 2) == 0 and (i % 50) != 0:
        new_codeword[i] = -1
    #n_erasures = int(np.floor(len(codeword) * epsilon))
    #erasure_indices = np.random.choice(range(len(codeword)), n_erasures, replace=False)
    #for idx in erasure_indices:
    #    new_codeword[idx] = -1
    return new_codeword

def get_error():
    epsilons = np.linspace(0.3, 0.7, num=41)
    channel_capacities = [1.0-eps for eps in epsilons]
    graph = Graph(H)
    errors = [ClusterGraph(graph, transmit(np.zeros(H.shape[1]), eps)).decode(int(1e7)) for eps in epsilons]
    with open("LDPC-uncoupled.pkl", 'w') as f:
        pickle.dump([channel_capacities, errors], f)

def plot_error():
    with open("LDPC-coupled.pkl", 'r') as f:
        coupled_c, coupled_e = pickle.load(f)
    with open("LDPC-uncoupled.pkl", 'r') as f:
        uncoupled_c, uncoupled_e = pickle.load(f)

    plt.plot(coupled_c, coupled_e, 'r', label="Error for Coupled LDPC")
    plt.plot(uncoupled_c, uncoupled_e, 'b', label="Error for Uncoupled LDPC")
    plt.title("Error over BEC with a (3,6)-regular & (3,6,64,3)-coupled LDPC")
    plt.xlabel("Channel capacity")
    plt.ylabel("Fraction of bits not recovered")
    plt.ylim(0.0,1.0)
    plt.plot((0.488,0.488),(0.0,1.0), 'r--', label="Shannon Limit, Coupled")
    plt.plot((0.5,0.5),(0.0,1.0), 'b--', label="Shannon Limit, Uncoupled")
    plt.legend()
    plt.show()

def plot_decoding_error():
    epsilon = 0.47
    all_block_errors = []
    all_iterations = []
    graph = Graph(H)
    for trial in range(1):
        iteration = 0
        decoder = ClusterGraph(graph, transmit(np.zeros(H.shape[1]), epsilon)).decode(int(1e7), yield_errors=True)
        for codeword in decoder:
            for i, bit in enumerate(codeword):
                if bit == -1:
                    all_block_errors.append(i % 64)
                    all_iterations.append(iteration)
            iteration += 1
    plt.hist2d(all_block_errors, all_iterations, bins=[16,max(all_iterations)/2])
    plt.colorbar()
    plt.title("Codeword Errors by Iteration")
    plt.xlabel("Block Number")
    plt.ylabel("Number of Errors")
    plt.show()

def main():
    plot_decoding_error()

if __name__ == "__main__":
    main()
