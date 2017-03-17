import numpy as np
import random

N_VARS = 256
N_CHECKS = 152
VAR_OUTDEG = 3
CHECK_OUTDEG = 6
N_BLOCKS = 16
COUPLING_WIDTH = 4

def generate_coupled_h(n_vars, n_checks, var_outdeg, check_outdeg, n_blocks, coupling_width):

    assert n_checks % (n_blocks + coupling_width - 1) == 0
    n_checks_per_block = n_checks / (n_blocks + coupling_width - 1)
    check_blocks = [range(n_checks_per_block*block_number, n_checks_per_block*(block_number+1)) \
              for block_number in range(n_blocks + coupling_width -1)]
    
    assert n_vars % n_blocks == 0
    n_vars_per_block = n_vars / n_blocks
    n_vars_eff = n_vars_per_block * (n_blocks + 2 * (coupling_width - 1)) 
    var_blocks = [range(n_vars_per_block*block_number, n_vars_per_block*(block_number+1)) \
                  for block_number in range(n_blocks + 2*(coupling_width - 1))]
    
    H = np.zeros((n_checks, n_vars_eff), dtype=int)

    assert len(var_blocks) - len(check_blocks) == coupling_width - 1
    assert n_vars_per_block * var_outdeg == n_checks_per_block * check_outdeg 
    
    #edges_filled[i] = the number of edges variable i has filled
    edges_filled = [0] * n_vars_eff
    full_vars = set([])

    for check_block_n, check_block in enumerate(check_blocks):
        # Spread out the edges in a given check block
        edges = range(n_checks_per_block * check_outdeg)

        # The first check_outdeg edges go to the first check in the block, etc.
        edge_to_check = {}
        for edge_n, edge in enumerate(edges):
            edge_to_check[edge_n] = check_block[edge_n / check_outdeg]
        
        # Get all possible variables (before accounting for full variables)
        available_variables = set([])
        # We can use all variables in the coupling_width blocks that we are coupling the checks to
        for var_block in var_blocks[check_block_n : check_block_n + coupling_width]:
            available_variables = available_variables.union(set(var_block))
        # Get rid of full variables
        available_variables = available_variables.difference(full_vars)     

        # Shuffle the edges and assign them one at a time
        for edge in np.random.permutation(edges):
            check = edge_to_check[edge]
            # Choose a variable from available variables randomly, assign an edge. If an edge already exists pick again
            while True:
                var = np.random.choice(list(available_variables))
                if H[check, var] != 1:
                    H[check, var] = 1
                    break
            edges_filled[var] += 1
            if edges_filled[var] == var_outdeg:
                available_variables.remove(var)
                full_vars.add(var)
    return H

H = generate_coupled_h(N_VARS, N_CHECKS, VAR_OUTDEG, CHECK_OUTDEG, N_BLOCKS, COUPLING_WIDTH)
for line in H:
    print sum(line)
np.savetxt("coupled_h.txt", H)
            
