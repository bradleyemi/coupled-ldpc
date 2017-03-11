### New implementation of LDPC and Spatially-Coupled LDPC for EE 376A
### Passes codewords through BEC (instead of BSC)

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
from copy import copy

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

G, H = loadLDPC("ldpc36-128.mat")

class Factor(object):
    def __init__(self, node_type, id, neighbors):
        self.node_type = node_type
        self.id = id
        self.neighbors = neighbors


class ClusterGraph(object):
    def __init__(self, H, codeword):
        self.H = H
        self.orig_codeword = codeword
        self.decoded_codeword = copy(codeword)
        
        self.n_vars = H.shape[1]
        self.n_checks = H.shape[0]

        self.vars = []
        self.checks = []
        self.var_to_check_msgs = {}
        self.check_to_var_msgs = {}

        for i in range(self.n_vars):
            neighbors = [k for k in range(self.n_checks) if H[k,i] == 1]
            self.vars.append(Factor("var", i, neighbors))

        for i in range(self.n_checks):
            neighbors = [k for k in range(self.n_vars) if H[i,k] == 1]
            self.checks.append(Factor("check", i, neighbors))

        for var in self.vars:
            for neighbor in var.neighbors:
                self.var_to_check_msgs[(var.id, neighbor)] = -1 #Message("var_to_check", -1)

        for check in self.checks:
            for neighbor in check.neighbors:
                self.check_to_var_msgs[(check.id, neighbor)] = -1 #Message("check_to_var", -1)

    def decode(self, iterations):
        print "decoding"
        for iteration in xrange(iterations):
            # First iterate over var -> check messages
            for source_dest, msg in self.var_to_check_msgs.iteritems():
                source = source_dest[0]
                dest = source_dest[1]
                source_neighbors = self.vars[source].neighbors
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
                source_neighbors = self.checks[source].neighbors
                in_msgs = [self.var_to_check_msgs[(neigh, source)] for neigh in source_neighbors if neigh != dest]
                # A check must have all incoming messages non-erased to be certain of the value of the destination
                # so -1 cannot be in the incoming messages. 
                # If all incoming messages are known, pass value that makes parity check = 0 (mod 2) and decode bit
                if -1 not in in_msgs:
                    self.check_to_var_msgs[source_dest] = sum(in_msgs) % 2
                    self.decoded_codeword[dest] = sum(in_msgs) % 2

        errors = len([c for c in self.decoded_codeword if c == -1])
        return float(errors) / H.shape[1]

def transmit(codeword, epsilon):
    # Take the convention that -1 is an erasure bit
    new_codeword = []
    for c in codeword:
        if np.random.uniform() > epsilon:
            new_codeword.append(c)
        else:
            new_codeword.append(-1)
    return new_codeword

def plot_error():
    epsilons = np.linspace(0.0, 1.0, num=100)
    channel_capacities = [1.0-eps for eps in epsilons]
    errors = [ClusterGraph(H, transmit(np.zeros(H.shape[1]), eps)).decode(20) for eps in epsilons]
    plt.plot(channel_capacities, errors, 'black', label="Error")
    plt.title("Error for varying capacity BECs with a (3,6)-regular LDPC coding scheme")
    plt.xlabel("Channel capacity")
    plt.ylabel("Fraction of bits not recovered")
    plt.ylim(0.0,1.0)
    plt.plot((0.5,0.5),(0.0,1.0), 'r--', label="Shannon limit")
    plt.plot((0.571,0.571),(0.0,1.0), 'g--', label="BP limit")
    plt.plot((0.512,0.512),(0.0,1.0), 'b--', label="MAP limit")
    plt.legend()
    plt.show()



def main():
    plot_error()

if __name__ == "__main__":
    main()
