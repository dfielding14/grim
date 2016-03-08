︠fbddea82-3103-42dd-bbf9-04c4a571eee5as︠
%auto
typeset_mode(True)
︡8835d497-f924-4f48-b1eb-fc094bd03276︡︡{"auto":true}︡{"done":true}
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
6. Heat flux magnitude $q$

* Mean variables      : $\rho_0$, $u_0$, $u^1_0$, $u^2_0$, $u^3_0$, $q_0$, 
* Perturbed variables : $\delta_\rho$, $\delta_u$, $\delta_{u^1}$, $\delta_{u^2}$, $\delta_{u^3}$, $\delta_q$


︡de554342-ffe2-4092-97f1-65f5f7e57804︡︡{"done":true,"md":"## $\\mathtt{balbusaur}$<br> ##\n#### A framework for automated linear analysis\n\n#### Initialize the operators $\\mathtt{d\\_dt()}$, $\\mathtt{d\\_dX1()}$ and $\\mathtt{d\\_dX2()}$ acting on variables:\n1. Density $\\rho$\n2. Internal energy $u$\n3. Velocity in $X^1$ direction $u^1$\n4. Velocity in $X^2$ direction $u^2$\n5. Velocity in $X^3$ direction $u^3$\n6. Heat flux magnitude $q$\n\n* Mean variables      : $\\rho_0$, $u_0$, $u^1_0$, $u^2_0$, $u^3_0$, $q_0$, \n* Perturbed variables : $\\delta_\\rho$, $\\delta_u$, $\\delta_{u^1}$, $\\delta_{u^2}$, $\\delta_{u^3}$, $\\delta_q$"}
︠81d0b45a-c4a1-4a7c-8b52-5a3198afcfbes︠
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

︡fd497d89-0cd8-4577-bd51-fc7aa62b361c︡︡{"done":true}
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
︠c953ed04-b6ab-4a98-a57f-2711f7d0d62b︠
# Spatiotemporal variables
t, omega, k1, k2, k3 = var('t, omega, k1, k2, k3')

# Constants:
# Gamma : Adiabatic index
# kappa : Heat conductivity
Gamma, kappa = var('Gamma, kappa')

# Background mean values: Symbolic variables
rho0, u0, u10, u20, u30 = var('rho0, u0, u10, u20, u30')

# Perturbations in space
delta_rho, delta_u, delta_u1, delta_u2, delta_u3 = \
    var('delta_rho, delta_u, delta_u1, delta_u2, delta_u3')

# Perturbations in time
delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt = \
    var('delta_rho_dt, delta_u_dt, delta_u1_dt, delta_u2_dt, delta_u3_dt')

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

# Introducing:
# chi : Thermal diffusivity
# nu  : kinematic viscosity
# cs  : sound speed
P   = (Gamma - 1)*u
T   = P/rho
cs  = sqrt(Gamma * P/ (rho + Gamma*u) )
chi   = phi * cs**2 * tau
kappa = rho * chi

# Inputs for numerical diagonalization for finite k modes
k1_num    = 2*pi
k2_num    = 4*pi

rho0_num = 1
u0_num   = 2
u10_num  = 0
u20_num  = 0
u30_num  = 0

Gamma_num = 4/3
︡8c9e36b4-17b4-4cc1-b74a-d3bb91050a64︡︡{"done":true}
︠c3616227-8585-41a4-ba9d-fe429350abe3si︠
%md

#### All the physics is below
︡e875cde8-d3ac-4dfb-8f25-df7e74479b02︡{"md":"\n#### All the physics is below\n"}︡
︠a7637a36-f643-4769-aed3-98f0e4e300b9s︠
Eqn_rho = linearize(d_dt(rho) + d_dX1(rho*u1) + d_dX2(rho*u2))
Eqn_u1  = linearize(d_dt(u1) + u1*d_dX1(u1) + u2*d_dX2(u1) + d_dX1(P)/rho )
Eqn_u2  = linearize(d_dt(u2) + u1*d_dX1(u2) + u2*d_dX2(u2) + d_dX2(P)/rho )
Eqn_u3  = linearize(d_dt(u3) + u1*d_dX1(u3) + u2*d_dX2(u3) )

# The paper uses an evolution equation for entropy, not an energy equation
Eqn_entropy = linearize(P/(Gamma-1) * (  d_dt(log(P/rho**Gamma))   \
                                    + u1*d_dX1(log(P/rho**Gamma))  \
                                    + u2*d_dX2(log(P/rho**Gamma))) \
                        + 0)

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
pretty_print("Analytic eigenvalues and eigenvectors in the $k_1, k_2 \\rightarrow 0$ limit : ", )
M.subs(k1=0, k2=0).eigenvectors_right()

# Numerical diagonalization:

M_numerical = M.subs(rho0=rho0_num, u0=u0_num, u10=u10_num, u20=u20_num, u30=u30_num, \
                     Gamma=Gamma_num, \
                     k1=k1_num, k2=k2_num \
                    )

M_numerical = M_numerical.change_ring(CDF)
eigenvecs   = M_numerical.eigenvectors_right()

pretty_print("Numerical eigenvalues and eigenvectors for $k_1 = $", k1_num, " , $k_2$ = ", k2_num, ":\n")


for i in xrange(len(eigenvecs)):
    print("--------------------------")
    print("Eigenvalue   = ", eigenvecs[i][0])
    print(delta_rho,  " = ", eigenvecs[i][1][0][0])
    print(delta_u1,   " = ", eigenvecs[i][1][0][1])
    print(delta_u2,   " = ", eigenvecs[i][1][0][2])
    print(delta_u3,   " = ", eigenvecs[i][1][0][3])
    print(delta_u,    " = ", eigenvecs[i][1][0][4])
︡57c4f275-b7b7-4ab1-80cb-95730b92cac1︡︡{"html":"<div align='center'>Linearized system : </div>","done":false}︡{"stdout":"\n\n","done":false}︡{"html":"<div align='center'>$\\displaystyle \\left(\\begin{array}{r}\n\\delta_{\\rho_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u1}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u2}_{\\mathit{dt}}} \\\\\n\\delta_{\\mathit{u3}_{\\mathit{dt}}} \\\\\n\\delta_{u_{\\mathit{dt}}}\n\\end{array}\\right)$  =  $\\displaystyle \\left(\\begin{array}{rrrrr}\n0 &amp; -i \\, k_{1} \\rho_{0} &amp; -i \\, k_{2} \\rho_{0} &amp; 0 &amp; 0 \\\\\n0 &amp; 0 &amp; 0 &amp; 0 &amp; \\frac{{\\left(-i \\, \\Gamma + i\\right)} k_{1}}{\\rho_{0}} \\\\\n0 &amp; 0 &amp; 0 &amp; 0 &amp; \\frac{{\\left(-i \\, \\Gamma + i\\right)} k_{2}}{\\rho_{0}} \\\\\n0 &amp; 0 &amp; 0 &amp; 0 &amp; 0 \\\\\n0 &amp; -i \\, \\Gamma k_{1} u_{0} &amp; -i \\, \\Gamma k_{2} u_{0} &amp; 0 &amp; 0\n\\end{array}\\right)$ $\\displaystyle \\left(\\begin{array}{r}\n\\delta_{\\rho} \\\\\n\\delta_{u_{1}} \\\\\n\\delta_{u_{2}} \\\\\n\\delta_{u_{3}} \\\\\n\\delta_{u}\n\\end{array}\\right)$</div>","done":false}︡{"stdout":"\n\n\n","done":false}︡{"html":"<div align='center'>Analytic eigenvalues and eigenvectors in the $k_1, k_2 \\rightarrow 0$ limit : </div>","done":false}︡{"html":"<div align='center'>[($\\displaystyle 0$, [$\\displaystyle \\left(1,\\,0,\\,0,\\,0,\\,0\\right)$, $\\displaystyle \\left(0,\\,1,\\,0,\\,0,\\,0\\right)$, $\\displaystyle \\left(0,\\,0,\\,1,\\,0,\\,0\\right)$, $\\displaystyle \\left(0,\\,0,\\,0,\\,1,\\,0\\right)$, $\\displaystyle \\left(0,\\,0,\\,0,\\,0,\\,1\\right)$], $\\displaystyle 5$)]</div>","done":false}︡{"html":"<div align='center'>Numerical eigenvalues and eigenvectors for $k_1 = $ $\\displaystyle 2 \\, \\pi$  , $k_2$ =  $\\displaystyle 4 \\, \\pi$ :\n</div>","done":false}︡{"stdout":"--------------------------\n('Eigenvalue   = ', 0.0)\n(delta_rho, ' = ', 1.0)\n(delta_u1, ' = ', 0.0)\n(delta_u2, ' = ', 0.0)\n(delta_u3, ' = ', 0.0)\n(delta_u, ' = ', 0.0)\n--------------------------\n('Eigenvalue   = ', 13.246117687728134*I)\n(delta_rho, ' = ', 0.33333333333333326 - 6.409875621278545e-17*I)\n(delta_u1, ' = ', -0.14054567378526128 + 2.7026408642169405e-17*I)\n(delta_u2, ' = ', -0.2810913475705225 + 5.40528172843388e-17*I)\n(delta_u3, ' = ', 0.0)\n(delta_u, ' = ', 0.8888888888888888)\n--------------------------\n('Eigenvalue   = ', -5.5477369469122726e-17 + 7.591448211103332e-16*I)\n(delta_rho, ' = ', 0.9831264779760727)\n(delta_u1, ' = ', 0.0005882405159709954 + 0.16361392549225612*I)\n(delta_u2, ' = ', -0.00029412025798556014 - 0.0818069627461281*I)\n(delta_u3, ' = ', 0.0)\n(delta_u, ' = ', 5.1165785135882766e-17 - 1.168379630072576e-16*I)\n--------------------------\n('Eigenvalue   = ', -6.921127172674474e-16 - 13.246117687728134*I)\n(delta_rho, ' = ', 0.3333333333333334 + 6.715343395427227e-17*I)\n(delta_u1, ' = ', 0.1405456737852614 - 3.030838367273343e-17*I)\n(delta_u2, ' = ', 0.2810913475705226 + 8.137810935710134e-17*I)\n(delta_u3, ' = ', 0.0)\n(delta_u, ' = ', 0.8888888888888888)\n--------------------------\n('Eigenvalue   = ', 0.0)\n(delta_rho, ' = ', 0.0)\n(delta_u1, ' = ', 0.0)\n(delta_u2, ' = ', 0.0)\n(delta_u3, ' = ', 1.0)\n(delta_u, ' = ', 0.0)\n","done":false}︡{"done":true}
︠d39e6e67-8c79-4c5a-8606-03aa2a4281dcs︠
M.eigenvectors_right()
︡065016cc-abfd-4e1c-ab25-fc8f4b6285aa︡︡{"html":"<div align='center'>[($\\displaystyle -\\sqrt{-\\frac{\\Gamma^{2} k_{1}^{2} u_{0}}{\\rho_{0}} - \\frac{\\Gamma^{2} k_{2}^{2} u_{0}}{\\rho_{0}} + \\frac{\\Gamma k_{1}^{2} u_{0}}{\\rho_{0}} + \\frac{\\Gamma k_{2}^{2} u_{0}}{\\rho_{0}}}$, [$\\displaystyle \\left(1,\\,\\frac{{\\left(i \\, \\Gamma^{2} - i \\, \\Gamma\\right)} k_{1} u_{0}}{\\sqrt{-\\frac{\\Gamma^{2} k_{1}^{2} u_{0}}{\\rho_{0}} - \\frac{\\Gamma^{2} k_{2}^{2} u_{0}}{\\rho_{0}} + \\frac{\\Gamma k_{1}^{2} u_{0}}{\\rho_{0}} + \\frac{\\Gamma k_{2}^{2} u_{0}}{\\rho_{0}}} \\rho_{0}^{2}},\\,-\\frac{i \\, \\sqrt{\\Gamma} \\sqrt{-\\Gamma + 1} k_{2} \\sqrt{u_{0}}}{\\sqrt{k_{1}^{2} + k_{2}^{2}} \\rho_{0}^{\\frac{3}{2}}},\\,0,\\,\\frac{\\Gamma u_{0}}{\\rho_{0}}\\right)$], $\\displaystyle 1$), ($\\displaystyle \\sqrt{-\\frac{\\Gamma^{2} k_{1}^{2} u_{0}}{\\rho_{0}} - \\frac{\\Gamma^{2} k_{2}^{2} u_{0}}{\\rho_{0}} + \\frac{\\Gamma k_{1}^{2} u_{0}}{\\rho_{0}} + \\frac{\\Gamma k_{2}^{2} u_{0}}{\\rho_{0}}}$, [$\\displaystyle \\left(1,\\,\\frac{{\\left(-i \\, \\Gamma^{2} + i \\, \\Gamma\\right)} k_{1} u_{0}}{\\sqrt{-\\frac{\\Gamma^{2} k_{1}^{2} u_{0}}{\\rho_{0}} - \\frac{\\Gamma^{2} k_{2}^{2} u_{0}}{\\rho_{0}} + \\frac{\\Gamma k_{1}^{2} u_{0}}{\\rho_{0}} + \\frac{\\Gamma k_{2}^{2} u_{0}}{\\rho_{0}}} \\rho_{0}^{2}},\\,\\frac{i \\, \\sqrt{\\Gamma} \\sqrt{-\\Gamma + 1} k_{2} \\sqrt{u_{0}}}{\\sqrt{k_{1}^{2} + k_{2}^{2}} \\rho_{0}^{\\frac{3}{2}}},\\,0,\\,\\frac{\\Gamma u_{0}}{\\rho_{0}}\\right)$], $\\displaystyle 1$), ($\\displaystyle 0$, [$\\displaystyle \\left(1,\\,0,\\,0,\\,0,\\,0\\right)$, $\\displaystyle \\left(0,\\,1,\\,-\\frac{k_{1}}{k_{2}},\\,0,\\,0\\right)$, $\\displaystyle \\left(0,\\,0,\\,0,\\,1,\\,0\\right)$], $\\displaystyle 3$)]</div>","done":false}︡{"done":true}
︠ff95f3fd-65c9-4405-ac1a-18ab711c271c︠









