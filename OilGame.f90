program OilGame

	use netcdf
	implicit none


    ! dimensions for state variables

    ! number of control variables for the optimization tool
    INTEGER :: MAX_NCtrlVar

    ! constants
    double precision, PARAMETER :: g_CloseZero = 1.0d-12
    double precision, PARAMETER :: g_pInfty = 1.0d+15



    ! derivative option for objective and constraint functions
    integer :: g_derOpt
    INTEGER, PARAMETER :: DERIVATIVE_OPT = 1

    integer :: g_NVar
    double precision, dimension(:), allocatable :: x0
    double precision, dimension(:), allocatable :: xsol
    double precision :: fvalue
    ! bound for states and controls
    double precision, dimension(:), allocatable :: g_LB, g_UB ! lower and upper bounds for variables for the optimization tool 

    ! NPSOL specific

	double precision, DIMENSION(:, :), allocatable :: g_linA 
	double precision, DIMENSION(:), allocatable :: g_linL, g_linU, g_nlnL, g_nlnU ! constraint for NPSOL: 
										! kmin <= F(k,L)-c+eps <= kmax

	INTEGER :: g_nclin ! number of linear constraints
	INTEGER :: g_ncnln ! number of nonlinear constraints
	INTEGER :: g_ldA ! row dimension of linear constraints 
					 ! g_linL <= g_linA*x <= g_linU
	INTEGER :: g_ldJ ! dimension of Jacobian for nonlinear constraints 
					 ! g_nlnL <= ConFun(x) <= g_nlnU
	INTEGER :: g_ldR ! row dimension of Hessian matrix
	! length of work vectors
	INTEGER :: g_leniw
	INTEGER :: g_lenrw
	INTEGER, DIMENSION(:), allocatable :: g_iw ! work vector
	double precision, DIMENSION(:), allocatable :: g_rw ! work vector
	INTEGER :: g_nctotl
	INTEGER, DIMENSION(:), allocatable :: g_istate 
	double precision, DIMENSION(:), allocatable :: g_bl, g_bu, g_clambda, g_grad, g_Con
	double precision, DIMENSION(:,:), allocatable :: g_cJac, g_Hess




	INTEGER, DIMENSION(:), allocatable :: DS
	character :: numstr*20
    integer :: max_DS, Nplayers, Nstates

	! character :: priceFilename*100
	double precision, dimension(:,:), allocatable :: valuearray, controlarray

    double precision, dimension(:,:), allocatable :: solutionarray
    double precision, dimension(:,:), allocatable :: vplusses
    integer, dimension (:,:), allocatable :: all_the_states_in_order

	double precision, dimension(:), allocatable :: LB_value, LB_investment, UB_value, UB_investment
	double precision :: rho, delta, A, Pop

	integer :: i, j, k, tot, counter
	integer :: n_past_iterations, n_iterations_to_do
	logical :: is_initialization_step
    integer :: solver_inform

    double precision, dimension(:), allocatable :: v_I_go_up, v_you_go_up
    integer, dimension(:,:), allocatable :: indices

	max_DS = 3
	rho = 0.95d0
	delta = 0.1
	A = 1
	Pop = 1

	call GETARG(1, numstr)
	read(numstr, *) Nplayers

	call GETARG(2, numstr)
	read(numstr, *) max_DS

	call GETARG(3, numstr)
	if (numstr == 'init' .or. numstr == '0') then
		is_initialization_step = .true.
		n_past_iterations = 0
	else
		is_initialization_step = .false.
		read(numstr, *) n_past_iterations
	end if

    Nstates = max_DS ** Nplayers

    allocate(DS(Nplayers))
    allocate(all_the_states_in_order(Nstates, Nplayers))

	allocate(valuearray(Nstates, Nplayers))
	allocate(controlarray(Nstates, Nplayers))

	! if (is_initialization_step .eqv. .true.) then
	! 	call DoInitialization()
	! else
	! 	call netcdfRead("Values.nc", "Values", max_DS, valuearray)
	! end if

	g_NVar = Nplayers * 2
    MAX_NCtrlVar = g_NVar
	g_nclin = 0
	g_ncnln = Nplayers * 2
    allocate(xsol(g_NVar))
    allocate(x0(g_NVar))
    allocate(solutionarray(Nstates, g_NVar))
    allocate(g_LB(g_NVar), g_UB(g_NVar))
    allocate(LB_investment(Nplayers), LB_value(Nplayers))
    allocate(UB_investment(Nplayers), UB_value(Nplayers))
    allocate(vplusses(Nplayers, Nplayers))

    LB_investment = 0
    LB_value = -g_pInfty

    UB_investment = g_pInfty
    UB_value = g_pInfty

    g_LB = (/LB_investment, LB_value/)
    g_UB = (/UB_investment, UB_value/)

    g_derOpt = 0
    call OptimInit()


	x0(1:Nplayers) = 1
	x0(Nplayers + 1:2 * Nplayers) = 10

	vplusses = 0


    allocate(v_I_go_up(Nplayers), v_you_go_up(Nplayers))
    allocate(indices(max_DS, max_DS))

    counter = 1
    do tot = max_DS * 2, 1, -1
        do i = max(tot - max_DS, 1), min(max_DS, tot - 1)
            
            all_the_states_in_order(counter,:) = (/tot - i, i/)
            indices(tot - i, i) = counter
            counter = counter + 1


        end do
    end do


    do i = 1, Nstates
            
        DS = all_the_states_in_order(i,:)

        write(*,*) 'starting i = ', i, ' DS = ', DS

        ! if (DS(1) == max_DS) then
        !     v_I_go_up = 9999
        ! else
        !     v_I_go_up = valuearray(indices(DS(1) + 1, DS(2)),:)
        ! end if
        ! if (DS(2) == max_DS) then
        !     v_you_go_up = 9999
        ! else
        !     v_you_go_up = valuearray(indices(DS(1), DS(2) + 1),:)
        ! end if
        vplusses(1,:) = 0 ! v_I_go_up
        vplusses(2,:) = 0 ! v_you_go_up
        write(*, *) 'vplusses'
        write(*, "(*(g0))") ((vplusses(j, k)," ",j=1, Nplayers), new_line("A"), k=1, Nplayers)

        write(*,*) 'x0'
        write(*,*) x0

        call OptimMethod(g_NVar, x0, xsol, fvalue, solver_inform)

        print *, "Inform is ", solver_inform
        ! print *, "Solved to within ", sum(abs(checkConstraints(g_ncnln, g_NVar, xsol)))
        ! print *, "Errors of ", abs(checkConstraints(g_ncnln, g_NVar, xsol))
        if (solver_inform > 0) then
            x0 = (/2, 2, 10, 10/)
            write(*,*) 'x0'
            write(*,*) x0
            call OptimMethod(g_NVar, x0, xsol, fvalue, solver_inform)
            print *, "Inform is ", solver_inform
            if (solver_inform > 0) then
                stop
            end if
        end if

        write(*,*) 'xsol'
        write(*,*) xsol

        write(*,*) 'Nstates'
        write(*,*) Nstates
        write(*,*) 'finished i = ', i

        valuearray(i,:) = xsol(Nplayers + 1:2*Nplayers)
        solutionarray(i,:) = xsol
    end do
    ! is_initialization_step = .false.

	call netcdfWrite("Values.nc", "Values", valuearray, Nstates, Nplayers)
	call netcdfWrite("Policies.nc", "Policies", controlarray, Nstates, Nplayers)

contains

	subroutine OptimInit()

		g_ldA = max(1, g_nclin)  ! row dimension of linear constraints 
								 ! g_linL <= g_linA*x <= g_linU
		g_ldJ = max(1, g_ncnln)  ! dimension of Jacobian for nonlinear constraints 
								 ! g_nlnL <= ConFun(x) <= g_nlnU
		g_ldR = g_NVar           ! row dimension of Hessian matrix
		allocate(g_linA(g_ldA, g_NVar))
		allocate(g_linL(g_ldA), g_linU(g_ldA))
		allocate(g_nlnL(g_ldJ), g_nlnU(g_ldJ)) 

		g_nctotl = g_NVar+g_nclin+g_ncnln 
		allocate(g_istate(g_nctotl))
		allocate(g_bl(g_nctotl), g_bu(g_nctotl), g_clambda(g_nctotl)) 
		allocate(g_grad(g_NVar), g_cJac(g_ldJ,g_NVar), g_Hess(g_ldR,g_NVar), g_Con(g_ldJ))

		! length of work vectors
		g_leniw = 3*g_NVar+g_nclin+2*g_ncnln
		g_lenrw = 2*g_NVar**2+g_NVar*g_nclin+2*g_NVar*g_ncnln+20*g_NVar+11*g_nclin+21*g_ncnln
		allocate(g_iw(g_leniw), g_rw(g_lenrw)) ! work vector

		g_nlnL(1:g_ldJ) = 0 ! 0<= F(k,L)-c-Ekplus <= 0
		g_nlnU(1:g_ldJ) = 0

		! set all bounds
		g_bl(1:g_NVar) = g_LB(1:g_NVar)
		g_bu(1:g_NVar) = g_UB(1:g_NVar)
		if (g_nclin > 0) then
			g_bl((g_NVar+1):(g_NVar+g_nclin)) = g_linL(1:g_nclin)
			g_bu((g_NVar+1):(g_NVar+g_nclin)) = g_linU(1:g_nclin)
		end if
		if (g_ncnln > 0) then
			g_bl((g_NVar+g_nclin+1):g_nctotl) = g_nlnL(1:g_ncnln)
			g_bu((g_NVar+g_nclin+1):g_nctotl) = g_nlnU(1:g_ncnln)
		end if

		! g_istate = 0 ! No constrains are expected to be active, only useful for Warm Start
		call npoptn('Cold Start')
		call npoptn('Nolist')
		call npoptn('Iteration limit = 1000')
		call npoptn('Print file = 0')
		if (g_derOpt == DERIVATIVE_OPT) then
			call npoptn('Derivative level = 3')  ! 3: all gradients of objective and constraints are given
		else
			call npoptn('Derivative level = 0')  ! 0: no gradient given
		end if

	end subroutine OptimInit


    subroutine OptimMethod(nvar, x0, xsol, fvalue, inform)

        integer :: iter, nvar
        double precision, intent(IN) :: x0(nvar) 
        double precision, intent(OUT) :: xsol(nvar), fvalue
        integer, intent(OUT) :: inform
        ! double precision, external :: dF_k
        ! EXTERNAL NPSOLConFun, NPSOLObjFun, dU_c

        xsol = x0

        call NPSOL(nvar,g_nclin,g_ncnln,g_ldA,g_ldJ,g_ldR,g_linA,g_bl,g_bu,NPSOLConFun,NPSOLObjFun,inform,& 
            iter,g_istate,g_Con,g_cJac,g_clambda,fvalue,g_grad,g_Hess,xsol,g_iw,g_leniw,g_rw,g_lenrw)

        select case (inform)
            case (9)
                print *, "An input parameter was invalid"
            case (7)        
                print *, "The function derivatives returned by funcon or fcn is incorrect"
            case (6)
    !           g_WorkerInform = WORKER_NPSOL_1st_ORDER_FAIL
                print *, "x does not satisfy the 1st-order optimality conditions"
                ! write(*,'(2(I4), 5(f7.3, :, ", "))') DS(1), DS(2), xsol
            case (4)
                print *, "The Major iteration limit was reached"
            case (3)
                print *, "The nonlinear constraints could not be satisfied"
            case (2)
                print *, "The linear constraints and bounds could not be satisfied"
    !       case (0)
    !           print *, "Get an optimal solution x"
            case default
    !           g_WorkerInform = WORKER_NPSOL_OTHER_FAIL
    !           print *, "Terminated with other reasons"
        end select

        ! NPSOL solves the minimal problem, but we are dealing with the maximal problem
        ! so we change the objective function f(x) to -f(x). Now we change the optimal 
        ! objective value -f(x^*) backwards f(x^*)
        fvalue = -fvalue 

        if (abs(fvalue) > g_pInfty) then
            print *, "function value is NaN"
        end if

    end subroutine OptimMethod

	! nonlinear constraint function for npsol: g_nlnL <= NPSOLConFun(x) <= g_nlnU
	subroutine NPSOLConFun(mode,ncnln,n,ldJ,needc,x,cc,cJac,nstate)

		INTEGER :: mode,ncnln,n,ldJ,nstate, i
		INTEGER :: needc(ncnln)
		double precision :: x(n), cc(ncnln), cJac(ldJ,*)
		double precision :: control(Nplayers), value(Nplayers)
		double precision :: price(Nplayers), profits(Nplayers), transition_hazard_rates(Nplayers)
		double precision :: Bellman(Nplayers), FOC(Nplayers)

		control = x(1:Nplayers)
		value = x(Nplayers + 1:Nplayers*2)

		price = pricing(control)
		profits = revenue(price, control) - cost(control)
		transition_hazard_rates = conditional_transition_haz(control)
		Bellman = rho * value - (profits + matmul(transition_hazard_rates, vplusses - spread(value, 1, Nplayers)))
		FOC = price + d_c(control) +   &
			sum((vplusses - spread(value, 1, Nplayers)) * jac_transition(control), dim = 1)

		cc = (/Bellman, FOC/)

	end subroutine NPSOLConFun

	! the objective function for NPSOL method
	subroutine NPSOLObjFun(mode,n,x,v,grad,nstate)

		INTEGER :: mode,n,nstate
		double precision :: v
		double precision :: x(n), grad(n)	
		integer :: i

		v = sum(x(1:Nplayers))

		! ! NPSOL solves the minimal problem, but we are dealing with the maximal problem
		! ! so we change the objective function f(x) to -f(x)
		v = - v


	end subroutine NPSOLObjFun

	function revenue(p, q)

		double precision, dimension(Nplayers), intent(IN) :: p(:), q(:)
		double precision :: revenue(Nplayers)
		revenue = p * q

	end function revenue

	function pricing(control)

		double precision, intent(IN) :: control(Nplayers)
		double precision :: total_quantity, pricing(Nplayers)

		total_quantity = sum(control)

		pricing = A - total_quantity / Pop

	end function pricing

	function cost(x)

		double precision, intent(IN) :: x(Nplayers)
		double precision :: cost(Nplayers)
		cost = 0.d0

	end function cost

	function d_c(x)

		double precision, intent(IN) :: x(Nplayers)
		double precision :: d_c(Nplayers)
		d_c = 0.d0

	end function d_c

	function conditional_transition_haz(control)

		double precision, intent(in) :: control(Nplayers)
		double precision :: conditional_transition_haz(Nplayers)

		conditional_transition_haz = delta * control

	end function conditional_transition_haz


	function jac_transition(control)

		double precision, intent(in) :: control(Nplayers)
		double precision :: jac_transition(Nplayers, Nplayers)
		integer :: i

		do i = 1, Nplayers
			jac_transition(i, i) = delta
		end do

	end function jac_transition

	! subroutine DoInitialization()

	! 	integer :: i, j
	! 	double precision :: initialguess_value, initialguess_policy

	! 	do i = 1, max_DS
	! 		do j = 1, max_DS
	! 			valuearray(i,j) = 1.d0
	! 		end do
	! 	end do

	! end subroutine DoInitialization

	subroutine netcdfRead(filename, variable_name, nStates, variable_read)

		character (len = *), intent(in) :: filename
		integer, intent(in) :: nStates
		character (len = *), intent(in) :: variable_name
		integer :: ncid
		integer :: variable_ID
		logical :: file_exists

		double precision, intent(out) :: variable_read(nStates, nStates)

		inquire(file = filename, exist = file_exists)
		if (file_exists .eqv. .false.) then
			print *, "File ", trim(adjustl(filename)), " not found"
			stop
		end if

		call netcdfcheck( nf90_open(filename, nf90_nowrite, ncid) )

		call netcdfcheck( nf90_inq_varid(ncid, variable_name, variable_ID) )
		call netcdfcheck( nf90_get_var(ncid, variable_ID, variable_read) )  

		call netcdfcheck( nf90_close(ncid) )

	end subroutine netcdfRead


	subroutine netcdfWrite(filename, variable_name, array, max_DS, Nplayers)
		
		character (len = *), intent(in) :: filename
		character (len = *), intent(in) :: variable_name
		integer, intent(in) :: max_DS, Nplayers
		double precision, intent(in) :: array(max_DS, max_DS)

		integer :: ncid
		integer :: variable_ID

		character (len = *), parameter :: Xname = "Dim1"
		character (len = *), parameter :: Yname = "Dim2"
		integer :: X_dimid, Y_dimid

		integer :: dimids(Nplayers)


		! Create the file, delete it if it already exists (clobber)
		call netcdfcheck( nf90_create(filename, nf90_clobber, ncid) )
		
		call netcdfcheck( nf90_def_dim(ncid, Xname, max_DS, X_dimid) )
		call netcdfcheck( nf90_def_dim(ncid, Yname, max_DS, Y_dimid) )

		dimids = (/ X_dimid, Y_dimid /)
		call netcdfcheck( nf90_def_var(ncid, variable_name, NF90_DOUBLE, dimids, variable_ID) )

		! End define mode.
		call netcdfcheck( nf90_enddef(ncid) )

		call netcdfcheck( nf90_put_var(ncid, variable_ID, array) )

		call netcdfcheck( nf90_close(ncid) )

	end subroutine netcdfWrite

	subroutine netcdfcheck(status)
		integer, intent (in) :: status
		
		if(status /= nf90_noerr) then 
			print *, trim(nf90_strerror(status))
			stop 2
		end if
	end subroutine netcdfcheck

end program OilGame