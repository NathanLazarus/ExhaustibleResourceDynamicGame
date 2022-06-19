import numpy as np
from casadi import *


Nplayers = 3
DS = [1, 2, 3]
entire_problem_max_DS_list = [4, 4, 4]
x0 = np.array([1, 1, 1, 1, 1, 1])
params = [0.04, 10, 1, 1, np.arange(1, Nplayers + 1)] # rho, delta, A, Pop, kappas
vplusses = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
settings = [False, 'CasADi']



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



[rho, delta, A, Pop, kappas] = params
[enforce_strict_zeros, solver] = settings

print(np.argwhere(DS == entire_problem_max_DS_list))
where_at_max_DS = np.argwhere(DS == entire_problem_max_DS_list)[:, 0]


if True: #solver == 'CasADi':
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

    get_sol_CasADi = nlpsol("get_sol_CasADi", "ipopt", nlp, {"ipopt.print_level": 0, "ipopt.tol": 1e-10, 'print_time':0}) #, "ipopt.max_iter": 10000, "ipopt.acceptable_constr_viol_tol": 1e-6, "ipopt.gamma_theta": 1e-4})
    # get_sol_CasADi = nlpsol("get_sol_CasADi", "ipopt", nlp, {"ipopt.print_level": 6}) #, "ipopt.max_iter": 10000, "ipopt.acceptable_constr_viol_tol": 1e-6, "ipopt.gamma_theta": 1e-4})
    solution = get_sol_CasADi(
        x0=x0,
        lbg=-1e-13,
        ubg=1e-13,
        lbx=0,
        ubx=1e12,
    )

    print(solution)

    if enforce_strict_zeros:
        solution['x'] = replace_structural_zeros(solution['x'], where_at_max_DS, offset=0)
