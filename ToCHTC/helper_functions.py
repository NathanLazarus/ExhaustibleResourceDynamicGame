from casadi import *
import numpy as np
import sys, os
import netCDF4 as nc


def states_summing_to_x(length, total_sum, min_vals, max_vals):
    if length == 1:
        yield (total_sum,)
    else:
        range_min = max(min_vals[len(min_vals) - length], total_sum - sum(max_vals[len(max_vals) - length + 1:]))
        range_max = min(total_sum - sum(min_vals[len(min_vals) - length + 1:]), max_vals[len(max_vals) - length])
        for value in range(range_min, range_max + 1):
            for permutation in states_summing_to_x(length - 1, total_sum - value, min_vals, max_vals):
                yield (value,) + permutation

def pricing(Nplayers, control, A, Pop):

    total_quantity = sum1(control)

    return A - total_quantity / Pop

def cost(Nplayers, control, kappas):

    return kappas * control ** 2 / 2


def d_cost(Nplayers, control, kappas):

    return kappas * control


def conditional_transition_haz(Nplayers, control, delta):

    return delta * control


def jac_transition(Nplayers, control, delta):

    jac_mat = DM.zeros((Nplayers, Nplayers))

    # # Could make this neater. Figure out how to assign to the diagonal.
    # for i in range(Nplayers):
    #     if DS[i] < entire_problem_max_DS_list[i]:
    #         jac_mat[i, i] = tran_par_1 * (RDstock[i] / sum1(RDstock)) ** tran_par_2
    #     else:
    #         jac_mat[i, i] = 0

    jac_mat[0:Nplayers ** 2:Nplayers+1] = delta

    return jac_mat

def profit_func(Nplayers, control, price, kappas):

    revenue = price * control
    return revenue - cost(Nplayers, control, kappas)

def d_revenue(Nplayers, control, price, A, Pop):

    dp_dq = -1
    return price + control * dp_dq # A - control - sum1(control)

def remove_structural_zeros(vec, zero_locs, offset=0):
    all_indices = np.arange(vec.shape[0])
    return vec[np.delete(all_indices, offset + zero_locs)]

def replace_structural_zeros(vec, zero_locs, offset=0):
    return DM(np.insert(vec, offset + zero_locs - np.arange(len(zero_locs)), 0))


def get_vplusses(Nplayers, DS, min_DS_list, valuearray, indices, exploit_symmetry):

    if exploit_symmetry == True:
        future_states = DS + np.eye(Nplayers, dtype = 'int32')
        mysort = np.argsort(future_states)
        reversesort = np.argsort(mysort)
        sorted_future_states = np.take_along_axis(future_states, mysort, axis = 1)
        # vplusses = np.take_along_axis(valuearray[indices[tuple((sorted_future_states - min_DS_list).T)], Nplayers:], reversesort, axis = 1).T
        # Need to test

    else:
        future_states = DS + np.eye(Nplayers, dtype = 'int32')
        vplusses = valuearray[indices[tuple((future_states - min_DS_list).T)], Nplayers:].T

    return vplusses


def do_optimization(Nplayers, DS, entire_problem_max_DS_list, x0, params, vplusses, settings, presolved_prices = None):

    [rho, delta, A, Pop, kappas] = params
    [enforce_strict_zeros] = settings

    where_at_max_DS = np.argwhere(DS == entire_problem_max_DS_list)[:, 0]

    control = SX.sym('control', Nplayers, 1)
    value = SX.sym('value', Nplayers, 1)

    if enforce_strict_zeros:
        x_list = [remove_structural_zeros(control, where_at_max_DS), value, ]
        control[where_at_max_DS] = 0
        x0 = remove_structural_zeros(x0, where_at_max_DS, offset=0)
    else:
        x_list = [control, value, ]
        x0[Nplayers + where_at_max_DS] = 1e-9

    x = vertcat(*x_list)

    obj = -sum1(value)

    price = pricing(Nplayers, control, A, Pop)
    profits = profit_func(Nplayers, control, price, kappas)

    transition_hazard_rates = conditional_transition_haz(Nplayers, control, delta)

    jac_tran = jac_transition(Nplayers, control, delta)

    # Constraints

    FOC = d_revenue(Nplayers, control, price, A, Pop) - d_cost(Nplayers, control, kappas) + sum2((vplusses - value) * jac_tran)

    if enforce_strict_zeros:
        FOC = remove_structural_zeros(FOC, where_at_max_DS)
    else:
        FOC[where_at_max_DS] = control[where_at_max_DS]

    Bellman = rho * value - (profits + (vplusses - value) @ transition_hazard_rates)

    constraint_list = [FOC, Bellman]
    constraint = vertcat(*constraint_list)

    nlp = {
        "x": x,
        "f": obj,
        "g": constraint,
    }

    # with suppress_stdout():
    solver = nlpsol("solver", "ipopt", nlp, {"ipopt.print_level": 0, "ipopt.tol": 1e-10, 'print_time':0}) #, "ipopt.max_iter": 10000, "ipopt.acceptable_constr_viol_tol": 1e-6, "ipopt.gamma_theta": 1e-4})
    # solver = nlpsol("solver", "ipopt", nlp, {"ipopt.print_level": 6}) #, "ipopt.max_iter": 10000, "ipopt.acceptable_constr_viol_tol": 1e-6, "ipopt.gamma_theta": 1e-4})
    solution = solver(
        x0=x0,
        lbg=-1e-13,
        ubg=1e-13,
        lbx=0,
        ubx=1e12,
    )

    if enforce_strict_zeros:
        solution['x'] = replace_structural_zeros(solution['x'], where_at_max_DS, offset=0)

    # print('Bellman', Bellman)
    # print('solution', solution)

    return solution



def netcdfWrite(filename, variable_name, array, Nrows, Ncolumns):

    data_out = nc.Dataset(filename, 'w', format='NETCDF4')
    X_dim = data_out.createDimension('X_dim', Nrows)
    Y_dim = data_out.createDimension('Y_dim', Ncolumns)
    data_out.createVariable(variable_name, 'f8', ('X_dim', 'Y_dim',))

    data_out[variable_name][:] = array
    


def print_endog(Nplayers, DS, entire_problem_max_DS_list, solution, params, vplusses, settings):


#     numberedVarnames = [
#         vectorvar + str(ind) for vectorvar in vectorvars for ind in range(1, Nplayers + 1)
#     ]
#     numberedVarnamesWithTol = ["tol"] + [
#         vectorvar + str(ind) for vectorvar in vectorvars for ind in range(1, Nplayers + 1)
#     ]

#     # print(dict(zip(numberedVarnames, np.array(solution["x"]).flatten().tolist())))

    vectorvars = ['control', 'value']

    ind = 0
    for vectorvar in vectorvars:
        globals()[vectorvar] = DM(solution[ind : ind + Nplayers])
        ind += Nplayers

    [rho, delta, A, Pop, kappas] = params

    where_at_max_DS = np.argwhere(DS == entire_problem_max_DS_list)[:, 0]


    price = pricing(Nplayers, control, A, Pop)
    profits = profit_func(Nplayers, control, price, kappas)

    transition_hazard_rates = conditional_transition_haz(Nplayers, control, delta)

    jac_tran = jac_transition(Nplayers, control, delta)

    # Constraints

    FOC = d_revenue(Nplayers, control, price, A, Pop) - d_cost(Nplayers, control, kappas) + sum2((vplusses - value) * jac_tran)

    FOC[where_at_max_DS] = control[where_at_max_DS]

    Bellman = rho * value - (
        profits
        - cost(Nplayers, control, kappas)
        + (vplusses - value) @ transition_hazard_rates
    )

    constraint_list = [FOC, Bellman]
    constraint = vertcat(*constraint_list)

    print('constraint', constraint)

    [x, y, z] = [control[0], control[1], control[2]]
    [a, b, c] = [value[0], value[1], value[2]]
    Q = sum1(control)

    # mathematica_constraints = [
    #     (-1 + Q) * x + 0.5 * x ** 2 + a * (0.05 + 0.1 * x + 0.1 * y + 0.1 * z),     
    #       0.127394 * y + 0.122472 * z, 1 * a + 10 * Q + 20 * x,     10, 
    #      10 * Q * y + 10 * y ** 2 + b * (0.5 + 1 * x + 1 * y + 1 * z),     
    #       0.801677 * x + 10 * y + 0.730289 * z, 1 * b + 10 * Q + 30 * y,     10, 
    #      10 * Q * z + 15 * z ** 2 + c * (0.5 + 1 * x + 1 * y + 1 * z),     
    #       0.522321 * x + 0.495687 * y + 10 * z, 1 * c + 10 * Q + 40 * z,     10, 
    #      Q,     x + y + z]

    # print('mathematica_constraints', mathematica_constraints)

    # bellman1_mathematica = [a/20, (1/10) * (1.99988 * 10e-7 - a) * x + (1 - Q) * x - x ** 2/2 + 1/10 * (1.27394 - a) * y + 1/10 * (1.22472 - a) * z]
    # print('bellman1_mathematica', bellman1_mathematica)
    # mybellman = [(0.05*a), (1 - Q) * x-(x ** 2/2)+(1.99988 * 10e-7 - a)*(0.1*x)+(1.27394-a)*(0.1*y)+(1.22472-a)*(0.1*z)]
    # print('mybellman', mybellman)
    
#     for i in range(len(numberedVarnames)):
#         globals()[numberedVarnames[i]] = np.array(solution["x"]).flatten().tolist()[i]


#     # DS = np.array(DS)

#     RDstock = state_to_RDstock(Nplayers, DS, gridstart, gridspace)
#     productivities = psi(Nplayers, RDstock, psi_elasticity)
#     quantity = labor * productivities
#     profits = profit_func(Nplayers, price, labor, productivities, w)

#     transition_hazard_rates = transition_haz(
#         Nplayers, investment, RDstock, DS, entire_problem_max_DS_list, tran_par_1, tran_par_2
#     )

#     total_transition_hazard_rate = sum1(transition_hazard_rates)
#     jac_tran = jac_transition_investment(
#         Nplayers, investment, RDstock, DS, entire_problem_max_DS_list, tran_par_1, tran_par_2
#     )

#     # l_basic = labor_stock - sum1(labor)
#     # c_basic = l_basic - sum1(investment_cost(Nplayers, investment, RDstock, kappa))
#     # y_total = ((B + 1) / B) * c_basic
#     y_total = w * labor_stock + sum1(profits) - investment_cost(Nplayers, investment, RDstock, kappa)

#     aggregate_p = sum1(price ** (1 - chi)) ** (1 / (1 - chi))


#     theta = chi * (1 - (price / aggregate_p) ** (1 - chi))

#     c_technological = y_total / (aggregate_p * (B + 1))


























#     # with np.printoptions(
#     #     precision=4, suppress=True, formatter={"float": "{:0.4f}".format}, linewidth=100
#     # ):
#     #     print("value ", value[:])
#     #     print("investment ", investment[:])
#     #     print("labor ", labor[:])
#     #     print("price ", price[:])
#     #     print("RDstock ", RDstock[:])
#     #     print("productivities ", productivities[:])
#     #     print("quantity ", quantity[:])
#     #     print("profits ", profits[:])
#     #     print("transition_hazard_rates ", transition_hazard_rates[:])
#     #     print("total_transition_hazard_rate ", total_transition_hazard_rate[:])
#     #     print("jac_tran ", jac_tran[:])
#     #     print("l_basic ", l_basic[:])
#     #     print("c_basic ", c_basic[:])
#     #     print("y_total ", y_total[:])
#     #     print("aggregate_p ", aggregate_p[:])
#     #     print("theta ", theta[:])
#     #     print("c_technological ", c_technological[:])