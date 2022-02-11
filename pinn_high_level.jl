using NeuralPDE, Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Quadrature,Cubature
import ModelingToolkit: Interval, infimum, supremum

# Parameters and variables
@parameters t, r
@variables q_s_bar(..), q_p_bar(..), q_n_bar(..), c_sn_bar(..), c_sp_bar(..), c_sp(..), c_sn(..)
Dt = Differential(t)
Dr = Differential(r)

# Constants
D_sp = 1e-14
D_sn = 5.3125e-14
F = 96487
S_p = 1
S_n = 0.5918
R_p = 1.7183e-5
R_n = 1.6192e-5

I = 1.0 #(A)

# Equations
eqs = [
    Dt(c_sp_bar(r,t)) ~ (-3)*(I/(F*S_p))/R_p,
    Dt(c_sn_bar(r,t)) ~ (-3)*(I/(F*S_n))/R_n,
    Dt(q_p_bar(r,t)) + 30*(D_sp/(R_p)^2)*q_s_bar(r,t) + (45/2)*((I/(F*S_p))/R_p^2) ~ 0,
    Dt(q_n_bar(r,t)) + 30*(D_sn/(R_n)^2)*q_s_bar(r,t) + (45/2)*((I/(F*S_n))/R_n^2) ~ 0, 
]

# Conditions
bcs = [
    Dr(c_sp(0,t)) ~ 0,
    Dr(c_sn(0,t)) ~ 0,
    D_sp*Dr(c_sp(R_p,t)) ~ -I/(F*S_p),
    D_sn*Dr(c_sn(R_p,t)) ~ -I/(F*S_n),  
    c_sp(r,0) ~ (-1/(3*R_p))*(I/(F*S_p))^2,
    c_sn(r,0) ~ (-1/(3*R_n))*(I/(F*S_n))^2,
    q_p_bar(r,t) ~ (-15*D_sp/R_p^2)*q_s_bar(r,t)^2-(45/4)*(1/R_p^2)*(1/(F*S_p))^2, 
    q_n_bar(r,t) ~ (-15*D_sn/R_n^2)*q_n_bar(r,t)^2-(45/4)*(1/R_n^2)*(1/(F*S_n))^2
]


# Domains
domains = [
    t ∈ Interval(0.0, 1.0),
    r ∈ Interval(0.0, R_n)
]

# Neural network
input_ = length(domains)
n = 15
chain = [FastChain(FastDense(input_, n, Flux.σ), FastDense(n,n,Flux.σ), FastDense(n,1)) for _ in 1:7]
initθ = map(c -> Float64.(c), DiffEqFlux.initial_params.(chain))

_strategy = QuadratureTraining()
discretization = PhysicsInformedNN(chain, _strategy, init_params=initθ)

@named pde_system = PDESystem(eqs,bcs, domains, [r,t], [q_s_bar(r,t), q_p_bar(r,t), q_n_bar(r,t), c_sn_bar(r,t), c_sp_bar(r,t), c_sp(r,t), c_sn(r,t)])
prob = discretize(pde_system, discretization)
sys_prob = symbolic_discretize(pde_system, discretization)

pde_inner_loss_functions = prob.f.f.loss_function.pde_loss_function.pde_loss_functions.contents
bcs_inner_loss_functions = prob.f.f.loss_function.bcs_loss_function.bc_loss_functions.contents

cb = function (p,l)
    println("loss: ", l)
    println("pde_losses: ", map(l_ -> l_(p), pde_inner_loss_functions))
    println("bcs_losses: ", map(l_ -> l_(p), bcs_inner_loss_functions))
    return false
end

res = GalacticOptim.solve(prob,BFGS(); cb = cb, maxiters=10)

phi = discretization.phi

