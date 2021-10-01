from casadi import *
from helper_functions import *
import numpy as np
import sys, os
import itertools
import netCDF4 as nc
import argparse


def split_into_segments(ary, Nsections, bounds=False):

    Ntotal = len(ary)
    Neach_section, extras = divmod(Ntotal, Nsections)
    section_sizes = ([0] +
                     extras * [Neach_section+1] +
                     (Nsections-extras) * [Neach_section])
    div_points = np.array(section_sizes).cumsum()

    sub_arys = []
    for i in range(Nsections):
        st = div_points[i]
        end = div_points[i + 1]
        if bounds == True:
            sub_arys.append([ary[st], ary[end - 1]])
        else:
            sub_arys.append(ary[st:end])

    return sub_arys


parser = argparse.ArgumentParser()
parser.add_argument('Nplayers', type = int, help = 'number of players')
parser.add_argument('Nblocks', type = int, help = 'number of blocks')
parser.add_argument('entire_problem_max_DS', type = int, help = 'size of blocks')
parser.add_argument('thisStage', type = int, help = 'current stage')
parser.add_argument('--maxworkers', type = int, help = 'maximum number of worker jobs')
parser.add_argument('-i', '--workerid', type = int, help = 'worker ID')

args = parser.parse_args()

globals().update(args.__dict__)


rho = 0.05
delta = 0.1
A = 1
Pop = 1
kappas = np.arange(1, Nplayers + 1)

exploit_symmetry = False
enforce_strict_zeros = True

settings = [enforce_strict_zeros]



Nblocks_by_dim = np.array([Nblocks] * Nplayers)
grand_scheme_min_DS_list = np.array([1] * Nplayers)

grand_scheme_jobs = []

counter = 0
for tot in range(sum(Nblocks_by_dim), sum(grand_scheme_min_DS_list) - 1, -1):
    grand_scheme_jobs.append(list(states_summing_to_x(Nplayers, tot, grand_scheme_min_DS_list, Nblocks_by_dim)))


thisStageJobs = grand_scheme_jobs[thisStage - 1]
thesejobs = split_into_segments(thisStageJobs, min(maxworkers, len(thisStageJobs)))[workerid - 1]


for thisjob in thesejobs:


    entire_problem_max_DS_list = np.array([entire_problem_max_DS] * Nplayers)
    entire_problem_min_DS_list = np.array([1] * Nplayers)

    DS_bounds = [split_into_segments(np.arange(lb, ub + 1), n, bounds=True)
        for lb, ub, n in zip(entire_problem_min_DS_list, entire_problem_max_DS_list, Nblocks_by_dim)]


    min_DS_list, max_DS_list = np.vstack([np.array(DS_bounds[i][thisjob[i] - 1]) for i in range(Nplayers)]).T

    print('max_DS_list', max_DS_list)
    print('min_DS_list', min_DS_list)

    params = [rho, delta, A, Pop, kappas]

    widths = max_DS_list - min_DS_list + 1
    Nstates = np.prod(widths)

    print(f'{Nstates = }')


    if exploit_symmetry == True:
        if len(np.unique(max_DS_list)) == 1 and len(np.unique(min_DS_list)) == 1:
            Niter = int((np.math.factorial(max(max_DS_list) - min(min_DS_list) + Nplayers) /
                (np.math.factorial(Nplayers)*np.math.factorial(max(max_DS_list) - min(min_DS_list)))))
        else:
            print("implementing symmetry on a state space that isn't a cube would take a lot more thought")
            sys.exit()
    else:
        Niter = Nstates

    # **** Getting upstream values
    # print(grand_scheme_jobs[thisStage - 2])

    # print(['_'.join(map(str, job)) for job in grand_scheme_jobs[thisStage - 2]])

    # for job in grand_scheme_jobs[thisStage - 2]:
    #     id = '_'.join(str(x) for x in list(job)) for job
    #     np.genfromtxt('upstreamvalues_' + str(Nplayers) + '_' + id + '.csv', delimiter=",")


    if thisStage > 1:
        # upstream_values_deprecated = np.vstack([
        #     np.genfromtxt('valuearray_' + str(Nplayers) + '_' + '_'.join(map(str, job)) + '.csv', delimiter=",")
        #      for job in grand_scheme_jobs[thisStage - 2]])
        filenames = ['valuearray_' + str(Nplayers) + '_' + '_'.join(map(str, job)) + '.nc' for job in grand_scheme_jobs[thisStage - 2]]
        upstream_values = np.vstack([
            nc.Dataset(filename)['Values'][:]
            if os.path.exists(filename)
            else sys.exit('Upstream File ' + filename + ' not found')
            for filename in filenames])

    else:
        upstream_values = np.array([])

    face_size = 0
    for i in range(Nplayers):
        face_size += np.prod(np.delete(widths, i)) # == sum(np.prod(widths) / i, i = 1, Nplayers)

    faces = np.zeros([face_size, Nplayers], dtype="int64")

    valuearray = np.zeros([Niter + face_size + 1, Nplayers * 2], dtype="float64") + 9999

    indices = np.zeros(widths + 1, dtype="int32") + Niter + face_size
    all_the_states_in_order = np.zeros([Niter, Nplayers], dtype="int64")


    face_counter = 0
    for which_face in range(Nplayers):
        for tot in range(sum(np.delete(max_DS_list, which_face)), sum(np.delete(min_DS_list, which_face)) - 1, -1):
            for state in states_summing_to_x(Nplayers - 1, tot, np.delete(min_DS_list, which_face), np.delete(max_DS_list, which_face)):
                state_np_array = np.insert(np.array(state), which_face, max_DS_list[which_face] + 1)
                faces[face_counter, :] = state_np_array
                face_counter += 1


    for num, face_el in enumerate(faces):
        # print('face_el', face_el)
        if np.any(entire_problem_max_DS_list - face_el < 0):
            valuearray[num + Niter, :Nplayers] = face_el
            indices[tuple(face_el - min_DS_list)] = num + Niter
        else:
            upval = upstream_values[np.all(upstream_values[:, :Nplayers] == face_el, axis=1)][:, Nplayers:]
            # print('upval', upval)
            valuearray[num + Niter, Nplayers:] = upval
            valuearray[num + Niter, :Nplayers] = face_el
            indices[tuple(face_el - min_DS_list)] = num + Niter

    counter = 0
    for tot in range(sum(max_DS_list), sum(min_DS_list) - 1, -1):
        for state in states_summing_to_x(Nplayers, tot, min_DS_list, max_DS_list):
            state_np_array = np.array(state)
            if not exploit_symmetry or np.array_equal(state_np_array, np.sort(state_np_array)):
                all_the_states_in_order[counter, :] = state_np_array
                indices[tuple(state_np_array - min_DS_list)] = counter
                counter += 1

    # print(all_the_states_in_order)
    # print(indices)
    # print(valuearray)

    # if Nplayers == 3:
    #     print('indices[0,0,0]', indices[0,0,0], valuearray[indices[0,0,0]])
    #     print('indices[1,1,1]', indices[1,1,1], valuearray[indices[1,1,1]])
    #     print('indices[2,2,2]', indices[2,2,2], valuearray[indices[2,2,2]])
    #     print('indices[2,2,1]', indices[2,2,1], valuearray[indices[2,2,1]])
    #     print('indices[1,2,2]', indices[1,2,2], valuearray[indices[1,2,2]])
    #     print('indices[1,1,2]', indices[1,1,2], valuearray[indices[1,1,2]])
    #     print('indices[0,1,2]', indices[0,1,2], valuearray[indices[0,1,2]])
    #     print('indices[1,2,1]', indices[1,2,1], valuearray[indices[1,2,1]])
    #     print('indices[2,1,1]', indices[2,1,1], valuearray[indices[2,1,1]])

    # if Nplayers == 2:
    #     print('indices[0,0]', indices[0,0], valuearray[indices[0,0]])
    #     print('indices[2,2]', indices[2,2], valuearray[indices[2,2]])
    #     print('indices[1,1]', indices[1,1], valuearray[indices[1,1]])
    #     print('indices[0,1]', indices[0,1], valuearray[indices[0,1]])
    #     print('indices[1,0]', indices[1,0], valuearray[indices[1,0]])
    #     print('indices[1,2]', indices[1,2], valuearray[indices[1,2]])
    #     print('indices[2,1]', indices[2,1], valuearray[indices[2,1]])
    #     print('indices[0,2]', indices[0,2], valuearray[indices[0,2]])
    #     print('indices[2,0]', indices[2,0], valuearray[indices[2,0]])


    initialguesses = np.array([0.6 / Nplayers, 10])

    x0_0 = np.repeat(initialguesses, Nplayers)

    Nvar = Nplayers * 2
    solarray = np.zeros([Niter, Nvar], dtype="float64") + x0_0

    # donearray = np.zeros([Nstates, 1],dtype='bool')

    for i in range(Niter):

        # if exploit_symmetry == True:
        #     if donearray[i] == True:
        #         continue

        DS = all_the_states_in_order[i, :].T
        vplusses = get_vplusses(Nplayers, DS, min_DS_list, valuearray, indices, exploit_symmetry)

        # print('**')
        # print(DS)
        smallest_error = 999

        bestsolution = []
        for j in solarray[: min(i + 1, Niter), :][::-1]:
            x0 = j
        
            solution = do_optimization(Nplayers, DS, entire_problem_max_DS_list, x0, params, vplusses, settings)

            # print("total constraint violation = ", sum1(fabs(solution["g"])))

            ee_error = sum1(fabs(solution["g"]))

            if ee_error / Nplayers < 1e-7:
                solarray[i, :] = solution["x"].T
                valuearray[i, Nplayers:] = solution["x"].T[Nplayers:]
                valuearray[i, :Nplayers] = DS
                break
            else:
                if ee_error < smallest_error:
                    smallest_error = ee_error
                    bestsolution = solution["x"].T

        else:
            if smallest_error / Nplayers > 1e-5:
                print('No solution found, minimum error was', smallest_error)
                sys.exit()
            else:
                solarray[i, :] = bestsolution
                valuearray[i, Nplayers:] = bestsolution[Nplayers:]
                valuearray[i, :Nplayers] = DS
        # print(solution['x'])

        # print_endog(Nplayers, DS, entire_problem_max_DS_list, solarray[i, :], params, vplusses, settings)

    # np.savetxt('states.csv', all_the_states_in_order, delimiter=',')
    netcdfWrite('solarray_' + str(Nplayers) + '_' + '_'.join(map(str, thisjob)) + '.nc', 'Sols', np.hstack([all_the_states_in_order, solarray]), Niter, Nplayers * 3)
    # np.savetxt('valuearray_' + str(Nplayers) + '_' + '_'.join(map(str, thisjob)) + '.csv', valuearray[:Niter], delimiter=',')
    netcdfWrite('valuearray_' + str(Nplayers) + '_' + '_'.join(map(str, thisjob)) + '.nc', 'Values', valuearray[:Niter], Niter, Nplayers * 2)
    # print(solarray)
    # print(valuearray)

