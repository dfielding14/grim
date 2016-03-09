︠fbddea82-3103-42dd-bbf9-04c4a571eee5as︠
%auto
typeset_mode(True)
︡47b846cd-b4af-4afd-8518-b6796e4cfa8b︡{"auto":true}︡
︠0d1a5418-4348-4f75-9af5-685a94ee8334i︠
%md
## $\mathtt{balbusaur}$<br> ##
#### A framework for automated linear analysis

#### Initialize the operators $\mathtt{d\_dt()}$, $\mathtt{d\_dX1()}$ and $\mathtt{d\_dX2()}$ acting on variables:
1. Density $\rho$
2. Internal energy $u$
3. Velocity in $X^1$ direction $u^1$
4. Velocity in $X^2$ direction $u^2$
5. Velocity in $X^3$ direction $u^3$


* Mean variables      : $\rho_0$, $u_0$, $u^1_0$, $u^2_0$, $u^3_0$
* Perturbed variables : $\delta_\rho$, $\delta_u$, $\delta_{u^1}$, $\delta_{u^2}$, $\delta_{u^3}$


︡de554342-ffe2-4092-97f1-65f5f7e57804︡︡{"done":true,"md":"## $\\mathtt{balbusaur}$<br> ##\n#### A framework for automated linear analysis\n\n#### Initialize the operators $\\mathtt{d\\_dt()}$, $\\mathtt{d\\_dX1()}$ and $\\mathtt{d\\_dX2()}$ acting on variables:\n1. Density $\\rho$\n2. Internal energy $u$\n3. Velocity in $X^1$ direction $u^1$\n4. Velocity in $X^2$ direction $u^2$\n5. Velocity in $X^3$ direction $u^3$\n\n\n* Mean variables      : $\\rho_0$, $u_0$, $u^1_0$, $u^2_0$, $u^3_0$\n* Perturbed variables : $\\delta_\\rho$, $\\delta_u$, $\\delta_{u^1}$, $\\delta_{u^2}$, $\\delta_{u^3}$"}
︠81d0b45a-c4a1-4a7c-8b52-5a3198afcfbe︠
def linearize(term):
    return taylor(term, (delta_rho, 0), \
                        (delta_u, 0),   \
                        (delta_u1, 0),  \
                        (delta_u2, 0),  \
                        (delta_u3, 0),  \
                        (delta_rho_dt, 0), \
                        (delta_u_dt, 0),   \
                        (delta_u1_dt, 0),  \
                        (delta_u2_dt, 0),  \
                        (delta_u3_dt, 0), 1 \
                 ).simplify_full()

def d_dX1(term):
    term  = Expression(SR, linearize(term))

    expr =   term.coefficient(delta_rho) * I * k1 * delta_rho \
           + term.coefficient(delta_u)   * I * k1 * delta_u   \
           + term.coefficient(delta_u1)  * I * k1 * delta_u1  \
           + term.coefficient(delta_u2)  * I * k1 * delta_u2  \
           + term.coefficient(delta_u3)  * I * k1 * delta_u3

    return expr

def d_dX2(term):
    term  = Expression(SR, linearize(term))

    expr =   term.coefficient(delta_rho) * I * k2 * delta_rho \
           + term.coefficient(delta_u)   * I * k2 * delta_u   \
           + term.coefficient(delta_u1)  * I * k2 * delta_u1  \
           + term.coefficient(delta_u2)  * I * k2 * delta_u2  \
           + term.coefficient(delta_u3)  * I * k2 * delta_u3  


    return expr

def d_dX3(term):
    term  = Expression(SR, linearize(term))

    expr =   term.coefficient(delta_rho) * I * k3 * delta_rho \
           + term.coefficient(delta_u)   * I * k3 * delta_u   \
           + term.coefficient(delta_u1)  * I * k3 * delta_u1  \
           + term.coefficient(delta_u2)  * I * k3 * delta_u2  \
           + term.coefficient(delta_u3)  * I * k3 * delta_u3  

    return expr

def d_dt(term):
    term  = Expression(SR, linearize(term))

    expr =   term.coefficient(delta_rho) * delta_rho_dt \
           + term.coefficient(delta_u)   * delta_u_dt   \
           + term.coefficient(delta_u1)  * delta_u1_dt  \
           + term.coefficient(delta_u2)  * delta_u2_dt  \
           + term.coefficient(delta_u3)  * delta_u3_dt

    return expr

︡0d94fe8f-4ce7-4f6a-91c4-e2036727b9c2︡︡{"done":true}
︠90256b3f-1117-4e18-9737-e09ad6feec2ei︠
%md
#### Options:

1. $\mathtt{EVOLVE\_B\_FIELDS}$       : 0 or 1
2. $\mathtt{CONDUCTION}$              : 0 or 1
3. $\mathtt{VISCOSITY}$               : 0 or 1
4. $\mathtt{FAKE\_EMHD}$              : 0 or 1
5. $\mathtt{TURN\_OFF\_MEAN\_B2}$      : 0 or 1
6. $\mathtt{TURN\_OFF\_K2\_PERTURBATIONS}$ : 0 or 1
7. $\mathtt{PRINCIPAL\_COEFFICIENTS}$ : 0 or 1
︡2653d9cb-060b-4e81-bc75-c640b576d091︡︡{"done":true,"md":"#### Options:\n\n1. $\\mathtt{EVOLVE\\_B\\_FIELDS}$       : 0 or 1\n2. $\\mathtt{CONDUCTION}$              : 0 or 1\n3. $\\mathtt{VISCOSITY}$               : 0 or 1\n4. $\\mathtt{FAKE\\_EMHD}$              : 0 or 1\n5. $\\mathtt{TURN\\_OFF\\_MEAN\\_B2}$      : 0 or 1\n6. $\\mathtt{TURN\\_OFF\\_K2\\_PERTURBATIONS}$ : 0 or 1\n7. $\\mathtt{PRINCIPAL\\_COEFFICIENTS}$ : 0 or 1"}
︠293cd666-2808-49e0-8989-5449141bd297︠
# Spatiotemporal variables
t, omega, k1, k2, k3 = var('t, omega, k1, k2, k3')

# Constants:
# gamma : Adiabatic index
# kappa : Heat conductivity
gamma, kappa = var('gamma, kappa')

# Background mean values: Symbolic variables
rho0, u0, u10, u20, u30 = var('rho0, u0, u10, u20, u30')

# Perturbations in space
delta_rho, delta_u, delta_u1, delta_u2, delta_u3 = \
    var('delta_rho, delta_u, delta_u1, delta_u2, delta_u3')

# Perturbations in time
delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt = \
    var('delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt')

# Inputs:

TURN_OFF_K2_PERTURBATIONS = 0

# Equilibrium states:
# rho0 and u0 are NEVER set to zero
u10 = 0
u20 = 0
u30 = 0
if (TURN_OFF_K2_PERTURBATIONS):
    k2 = 0

rho = rho0 + delta_rho
u   = u0   + delta_u
u1  = u10  + delta_u1
u2  = u20  + delta_u2
u3  = u30  + delta_u3

# Introducing:
# cs  : sound speed
P        = (gamma - 1)*u
T        = P/rho
cs       = sqrt(gamma * P/ (rho + gamma*u) )
Lambda   = function('Lambda', rho, u)

# Inputs for numerical diagonalization for finite k modes
k1_num    = 2*pi
k2_num    = 4*pi

rho0_num = 1
u0_num   = 2
u10_num  = 0
u20_num  = 0
u30_num  = 0

gamma_num = 4/3
︡63288c65-742f-4cbe-b4e5-192628b5d086︡︡{"done":true}
︠c953ed04-b6ab-4a98-a57f-2711f7d0d62bi︠
# Spatiotemporal variables
t, omega, k1, k2, k3 = var('t, omega, k1, k2, k3')

# Constants:
# Gamma : Adiabatic index
# kappa : Heat conductivity
Gamma, kappa = var('Gamma, kappa')

# Background mean values: Symbolic variables
rho0, u0, u10, u20, u30, q0 = var('rho0, u0, u10, u20, u30, q0')

# Perturbations in space
delta_rho, delta_u, delta_u1, delta_u2, delta_u3, delta_q = \
    var('delta_rho, delta_u, delta_u1, delta_u2, delta_u3, delta_q')

# Perturbations in time
delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt, delta_q_dt = \
    var('delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt, delta_q_dt')

# Inputs:

CONDUCTION                = 1
TURN_OFF_K2_PERTURBATIONS = 0

# Equilibrium states:
# rho0 and u0 are NEVER set to zero
u10 = 0
u20 = 0
u30 = 0
if (TURN_OFF_K2_PERTURBATIONS):
    k2 = 0

rho = rho0 + delta_rho
u   = u0   + delta_u
u1  = u10  + delta_u1
u2  = u20  + delta_u2
u3  = u30  + delta_u3
q   = q0   + delta_q

# Introducing:
# chi : Thermal diffusivity
# nu  : kinematic viscosity
# cs  : sound speed
P   = (Gamma - 1)*u
T   = P/rho
cs  = sqrt(Gamma * P/ (rho + Gamma*u) )

# Inputs for numerical diagonalization for finite k modes
k1_num    = 2*pi
k2_num    = 4*pi

rho0_num = 1
u0_num   = 2
u10_num  = 0
u20_num  = 0
u30_num  = 0
q0_num   = 0 ##<---- what goes here?

Gamma_num = 4/3
%md

#### All the physics is below
︡55b7eb28-f2b9-4f0d-8626-7ea38438972d︡︡{"hide":"input"}︡{"md":"\n#### All the physics is below","done":false}︡{"done":true}
︠a7637a36-f643-4769-aed3-98f0e4e300b9︠
Eqn_rho = linearize(d_dt(rho) + d_dX1(rho*u1) + d_dX2(rho*u2))
Eqn_u1  = linearize(d_dt(u1) + u1*d_dX1(u1) + u2*d_dX2(u1) + d_dX1(P)/rho )
Eqn_u2  = linearize(d_dt(u2) + u1*d_dX1(u2) + u2*d_dX2(u2) + d_dX2(P)/rho )
Eqn_u3  = linearize(d_dt(u3) + u1*d_dX1(u3) + u2*d_dX2(u3) )

# The paper uses an evolution equation for entropy, not an energy equation
Eqn_entropy = linearize(P/(gamma-1) * (  d_dt(log(P/rho**gamma))   \
                                    + u1*d_dX1(log(P/rho**gamma))  \
                                    + u2*d_dX2(log(P/rho**gamma))) \
                        + rho**2 * Lambda)

Eqns          = [Eqn_rho==0, Eqn_u1==0, Eqn_u2==0, Eqn_u3==0, Eqn_entropy==0]
delta_vars    = [delta_rho, delta_u1, delta_u2, delta_u3, delta_u]
delta_vars_dt = [delta_rho_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt, delta_u_dt]

solutions = solve(Eqns, delta_vars_dt, solution_dict=True)

solns_delta_vars_dt = []
for dvar_dt in delta_vars_dt:
    solns_delta_vars_dt.append(solutions[0][dvar_dt])

M = jacobian(solns_delta_vars_dt, delta_vars)
M = M.apply_map(lambda x : x.simplify_full())

pretty_print("Linearized system : ", )
print("\n")
pretty_print(Matrix(delta_vars_dt).transpose(), " = ", M, Matrix(delta_vars).transpose())
print("\n\n")
#pretty_print("Analytic eigenvalues and eigenvectors in the $k_1, k_2 \\rightarrow 0$ limit : ", )
#M.subs(k1=0, k2=0).eigenvectors_right()

# Numerical diagonalization:

#M_numerical = M.subs(rho0=rho0_num, u0=u0_num, u10=u10_num, u20=u20_num, u30=u30_num, \
#                     Gamma=Gamma_num, \
#                     k1=k1_num, k2=k2_num \
#                    )

#M_numerical = M_numerical.change_ring(CDF)
#eigenvecs   = M_numerical.eigenvectors_right()

#pretty_print("Numerical eigenvalues and eigenvectors for $k_1 = $", k1_num, " , $k_2$ = ", k2_num, ":\n")


#for i in xrange(len(eigenvecs)):
#    print("--------------------------")
#    print("Eigenvalue   = ", eigenvecs[i][0])
#    print(delta_rho,  " = ", eigenvecs[i][1][0][0])
#    print(delta_u1,   " = ", eigenvecs[i][1][0][1])
#    print(delta_u2,   " = ", eigenvecs[i][1][0][2])
#    print(delta_u3,   " = ", eigenvecs[i][1][0][3])
#    print(delta_u,    " = ", eigenvecs[i][1][0][4])
︡d867c10d-df14-4377-bb8c-2dfa2133ca2b︡︡{"html":"<div align='center'>Linearized system : </div>","done":false}︡{"stdout":"\n\n","done":false}︡{"html":"<div align='center'>$\\displaystyle \\left(\\begin{array}{r}\n\\delta_{\\rho_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u1}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u2}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u3}_{\\mathit{dt}}} \\\\\n\\delta_{u_{\\mathit{dt}}}\n\\end{array}\\right)$  =  $\\displaystyle \\left(\\begin{array}{rrrrr}\n0 &amp; -i \\, k_{1} \\rho_{0} &amp; -i \\, k_{2} \\rho_{0} &amp; 0 &amp; 0 \\\\\n0 &amp; 0 &amp; 0 &amp; 0 &amp; \\frac{{\\left(-i \\, \\gamma + i\\right)} k_{1}}{\\rho_{0}} \\\\\n0 &amp; 0 &amp; 0 &amp; 0 &amp; \\frac{{\\left(-i \\, \\gamma + i\\right)} k_{2}}{\\rho_{0}} \\\\\n0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 \\\\\n-\\rho_{0}^{2} D[0]\\left(\\Lambda\\right)\\left(\\rho_{0}, u_{0}\\right) - 2 \\, \\rho_{0} \\Lambda\\left(\\rho_{0}, u_{0}\\right) &amp; -i \\, \\gamma k_{1} u_{0} &amp; -i \\, \\gamma k_{2} u_{0} &amp; 0 &amp; -\\rho_{0}^{2} D[1]\\left(\\Lambda\\right)\\left(\\rho_{0}, u_{0}\\right)\n\\end{array}\\right)$ $\\displaystyle \\left(\\begin{array}{r}\n\\delta_{\\rho} \\\\\n\\delta_{u_{1}} \\\\\n\\delta_{u_{2}} \\\\\n\\delta_{u_{3}} \\\\\n\\delta_{u}\n\\end{array}\\right)$</div>","done":false}︡{"stdout":"\n\n\n","done":false}︡{"done":true}
︠d39e6e67-8c79-4c5a-8606-03aa2a4281dc︠
M.eigenvectors_right()
︡178e1c57-f37f-4811-a3b8-c84859544e66︡{"done":false,"stderr":"Error in lines 1-1\n"}︡{"done":false,"stderr":"Traceback (most recent call last):\n  File \"/projects/sage/sage-6.10/local/lib/python2.7/site-packages/smc_sagews/sage_server.py\", line 905, in execute\n    exec compile(block+'\\n', '', 'single') in namespace, locals\n  File \"\", line 1, in <module>\n  File \"sage/matrix/matrix_symbolic_dense.pyx\", line 287, in sage.matrix.matrix_symbolic_dense.Matrix_symbolic_dense.eigenvectors_right (/projects/sage/sage-6.10/src/build/cythonized/sage/matrix/matrix_symbolic_dense.c:3424)\n    return self.transpose().eigenvectors_left()\n  File \"sage/matrix/matrix_symbolic_dense.pyx\", line 255, in sage.matrix.matrix_symbolic_dense.Matrix_symbolic_dense.eigenvectors_left (/projects/sage/sage-6.10/src/build/cythonized/sage/matrix/matrix_symbolic_dense.c:2899)\n    [evals,mults],evecs=self.transpose()._maxima_(maxima).eigenvectors()._sage_()\n  File \"/projects/sage/sage-6.10/local/lib/python2.7/site-packages/sage/interfaces/interface.py\", line 632, in __call__\n    return self._obj.parent().function_call(self._name, [self._obj] + list(args), kwds)\n  File \"/projects/sage/sage-6.10/local/lib/python2.7/site-packages/sage/interfaces/interface.py\", line 533, in function_call\n    return self.new(s)\n  File \"/projects/sage/sage-6.10/local/lib/python2.7/site-packages/sage/interfaces/interface.py\", line 308, in new\n    return self(code)\n  File \"/projects/sage/sage-6.10/local/lib/python2.7/site-packages/sage/interfaces/interface.py\", line 243, in __call__\n    return cls(self, x, name=name)\n  File \"/projects/sage/sage-6.10/local/lib/python2.7/site-packages/sage/interfaces/interface.py\", line 670, in __init__\n    raise TypeError(x)\nTypeError: ECL says: Console interrupt.\n"}︡{"done":true}︡
︠d03a40b1-fe56-44dc-aaf1-d8875137733b︠
︠27016e95-2474-46f3-aa37-84d0d58de8be︠

linearize(q)
︡908a0f58-8bc8-445a-a312-cd7982ebc2b3︡︡{"html":"<div align='center'>$\\displaystyle \\delta_{\\rho} D[0]\\left(q\\right)\\left(\\rho_{0}, u_{0}, 0, 0\\right) + \\delta_{u} D[1]\\left(q\\right)\\left(\\rho_{0}, u_{0}, 0, 0\\right) + \\delta_{\\rho_{\\mathit{dt}}} D[2]\\left(q\\right)\\left(\\rho_{0}, u_{0}, 0, 0\\right) + \\delta_{u_{\\mathit{dt}}} D[3]\\left(q\\right)\\left(\\rho_{0}, u_{0}, 0, 0\\right) + q\\left(\\rho_{0}, u_{0}, 0, 0\\right)$</div>","done":false}︡{"done":true}
︠9b781967-b0e9-48ca-a619-d4fbf5506cdb︠










