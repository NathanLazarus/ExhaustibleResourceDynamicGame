from casadi import *
from helper_functions import *
import numpy as np
import sys, os
import itertools
import netCDF4 as nc
import argparse

# def has_valid_swap(state, lbounds, ubounds):


Nplayers = 3
max_DS_list = np.array([4, 4, 3])
min_DS_list = np.array([2, 2, 2])

print(f'{len(np.unique(max_DS_list)) =}' )

exploit_symmetry = False
enforce_strict_zeros = True


Nstates = np.prod(max_DS_list - min_DS_list + 1)

if exploit_symmetry == True:
    if len(np.unique(max_DS_list)) == 1:
        Niter = (np.math.factorial(max_DS - min_DS + Nplayers) /
            (np.math.factorial(Nplayers)*np.math.factorial(max_DS - min_DS)))
    else:
        print('_')
#         howbig = 0
#         for tot in range(max_DS * Nplayers, 1, -1):
#             for state in states_summing_to_x(Nplayers, tot, [1] * Nplayers, [max_DS] * Nplayers):
#                 state_np_array = np.array(state)
#                 if np.array_equal(state_np_array, np.sort(state_np_array)):
#                     howbig += 1

#         print(f'{howbig =}')
#         Niter = howbig
else:
    Niter = Nstates


# indices = np.zeros([max_DS + 2] * Nplayers, dtype="int32") + Niter
all_the_states_in_order = np.zeros([Niter, Nplayers], dtype="int64")
valuearray = np.zeros([Niter + 1, Nplayers], dtype="float64") + 9999

# min_DS = 1

# max_DS_list = [5, 2, 4]
# min_DS_list = [3, 1, 1]

# print(list(states_summing_to_x(3, 11, min_DS_list, max_DS_list)))
# print(list(states_summing_to_x(3, 10, min_DS_list, max_DS_list)))
# print(list(states_summing_to_x(3, 9, min_DS_list, max_DS_list)))

print(list(states_summing_to_x(2, 9, min_DS_list, max_DS_list)))

counter = 0
for tot in range(sum(max_DS_list), sum(min_DS_list) - 1, -1):
    # print(f'{tot =}')
    # print(f'{counter =}')
    for state in states_summing_to_x(Nplayers, tot, min_DS_list, max_DS_list):
        # print(state)
        state_np_array = np.array(state)
        if not exploit_symmetry or np.array_equal(state_np_array, np.sort(state_np_array)):
            all_the_states_in_order[counter, :] = state_np_array
            print(np.array_equal(state_np_array, np.sort(state_np_array)))
            # indices[state] = counter
            counter += 1

print(all_the_states_in_order)
sys.exit()
