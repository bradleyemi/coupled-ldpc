import numpy as np
import random
import cPickle as pickle

M = 512
VAR_OUTDEG = 3
CHECK_OUTDEG = 6
N_BLOCKS = 64
COUPLING_WIDTH = 3

def generate_coupled_h(M, var_outdeg, check_outdeg, n_blocks, coupling_width):

    n_checks_per_block = M * var_outdeg / check_outdeg
    check_blocks = [range(n_checks_per_block*block_number, n_checks_per_block*(block_number+1)) \
              for block_number in range(n_blocks + coupling_width - 1)]
    n_checks_eff = (n_blocks + coupling_width - 1)*n_checks_per_block
    
    var_blocks = [range(M*block_number, M*(block_number+1)) \
                  for block_number in range(n_blocks)]
    n_vars = M*n_blocks
    
    H = np.zeros((n_checks_eff, n_vars), dtype=int)

    #edges_filled[i] = the number of edges check i has filled
    edges_filled = [0] * n_checks_eff
    full_checks = set([])

    for var_block_n, var_block in enumerate(var_blocks):
        print var_block_n
        edges = range(M * var_outdeg)
        edge_to_var = {}
        for edge_n, edge in enumerate(edges):
            edge_to_var[edge_n] = var_block[edge_n / var_outdeg]

        available_checks = set([])
        for check_block in check_blocks[var_block_n : var_block_n + coupling_width]:
            available_checks = available_checks.union(set(check_block))
        available_checks = available_checks.difference(full_checks)

        for edge in np.random.permutation(edges):
            var = edge_to_var[edge]
            while True:
                check = np.random.choice(list(available_checks))
                if H[check, var] != 1:
                    H[check, var] = 1
                    break
            edges_filled[check] += 1
            if edges_filled[check] == check_outdeg:
                available_checks.remove(check)
                full_checks.add(check)
    return H

H = generate_coupled_h(M, VAR_OUTDEG, CHECK_OUTDEG, N_BLOCKS, COUPLING_WIDTH)
print H.shape, H.sum()
np.save("coupled_h_medium.npy", H)
            
