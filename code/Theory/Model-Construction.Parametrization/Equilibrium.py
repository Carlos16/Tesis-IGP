from sympy import Matrix,lambdify,symbols,simplify
# Definitions of variables in a general setting 

r = symbols('r')
K = symbols('K')
a_1 = symbols('a_1')
a_2 = symbols('a_2')
a_3 = symbols('a_3')
e_1 = symbols('e_1')
e_2 = symbols('e_2')
e_3 = symbols('e_3')
q_1 = symbols('q_1')
q_2 = symbols('q_2')
m_P = symbols('m_P')
m_C = symbols('m_C')
m_R = symbols('m_R')

R = symbols('R')
C = symbols('C')
P = symbols('P')
t_hp = symbols('t_hp')
t_hc = symbols('t_hc')
ha1 = symbols('ha1')
ha2 = symbols('ha2')
ha3 = symbols('ha3')


# Matrix for solving to equilibrium
E_matrix= Matrix([[r/K,a_1/m_C,a_2/m_P],[e_1*a_1/m_C,0,-a_3/m_P],[a_2*e_2,a_3*e_3,0]])
B = Matrix([r,q_1,q_2*m_P])

# little code for computing the equilibrium values using cramer's rule
def Cramers_L_solver(E,B,j):
    E_prime= E[:,:]
    for i in range(E.rows):
        E_prime[i,j]=B[i]
    x_j=E_prime.det()/E.det()
    x_j=simplify(x_j)
    return x_j
#################

#Finding equilibrium values
R_eq = Cramers_L_solver(E_matrix,B,0)
C_eq = Cramers_L_solver(E_matrix,B,1)
P_eq = Cramers_L_solver(E_matrix,B,2)


##Equilibrium functions
set_R_eq = lambdify((m_P,m_C,K,q_1,q_2,r,a_1,a_2,a_3,e_1,e_2,e_3),R_eq)

set_C_eq = lambdify((m_P,m_C,K,q_1,q_2,r,a_1,a_2,a_3,e_1,e_2,e_3),C_eq)

set_P_eq = lambdify((m_P,m_C,K,q_1,q_2,r,a_1,a_2,a_3,e_1,e_2,e_3),P_eq)



#RM case


##Expressions
EqR_equation = (1+t_hc*R*ha1) *(1+t_hp*R*ha2)*r*(1-R/K)  - (ha1*C*(1+t_hp*R*ha2) + ha2*P*(1+t_hc*R*ha1)) 

EqP_equation = (1+t_hp*ha3*C)*e_1*ha1*R - (1+t_hc*ha1*R)*ha3*P  - q_1 *(1+t_hp*ha3*C)*(1+t_hc*ha1*R)
EqC_equation = (1+t_hp*ha3*C) *e_2*ha2*R + (1+t_hp*ha2*R)*e_3*ha3*C- q_2 *(1+t_hp*ha2*R)*(1+t_hp*ha3*C) 


##Convert to functions

EqR_equation = lambdify((K,q_1,q_2,r,ha1,ha2,ha3,e_1,e_2,e_3,t_hc,t_hp,R,C,P),EqR_equation)
EqC_equation = lambdify((K,q_1,q_2,r,ha1,ha2,ha3,e_1,e_2,e_3,t_hc,t_hp,R,C,P),EqC_equation)
EqP_equation = lambdify((K,q_1,q_2,r,ha1,ha2,ha3,e_1,e_2,e_3,t_hc,t_hp,R,C,P),EqP_equation)


