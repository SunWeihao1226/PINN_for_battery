# Reference: https://neuralpde.sciml.ai/stable/pinn/system/
#            Marquis, S. G., Sulzer, V., Timms, R., Please, C. P., & Chapman, S. J. (2019). An asymptotic derivation of a single particle model with electrolyte. arXiv [physics.chem-ph]. Opgehaal van http://arxiv.org/abs/1905.12553

using NeuralPDE, Flux, ModelingToolkit, GalacticOptim, Optim, DiffEqFlux
using Quadrature,Cubature
import ModelingToolkit: Interval, infimum, supremum

@parameters t, r
@variables c_sp(..), c_sn(..)
Dt = Differential(t)
Dr = Differential(r)

# Constants
C_p = 0.0442
C_n = 0.1134
a_n = 1.8
a_p = 1.5
y_n = 1
y_p = 2.0501
L_p = 0.4444
L_n = 0.4444
I = 1.0


# Equations
eqs = [
    C_p*Dt(c_sp(r,t)) ~ (1/(r)^2)*Dr(((r)^2)*Dr(c_sp(r,t))),
    C_n*Dt(c_sp(r,t)) ~ (1/(r)^2)*Dr(((r)^2)*Dr(c_sn(r,t)))

]

# Conditions
bcs = [
    Dr(c_sp(0,t)) ~ 0,
    Dr(c_sn(0,t)) ~ 0,
    ((a_p*y_p)/(C_p))*Dr(c_sp(1,t)) ~ I/(L_p),
    ((a_n*y_n)/(C_n))*Dr(c_sn(1,t)) ~ -I/(L_n), 
    c_sp(r,0) ~ 0.6,
    c_sn(r,0) ~ 0.8
]


# Domains
domains = [
    t ∈ Interval(0.0, 1.0),
    r ∈ Interval(0.0, L_n)
]


# Neural network
input_ = length(domains)
n = 15
chain = [FastChain(FastDense(input_, n, Flux.σ), FastDense(n,n,Flux.σ), FastDense(n,1)) for _ in 1:2]
initθ = map(c -> Float64.(c), DiffEqFlux.initial_params.(chain))

_strategy = QuadratureTraining()
discretization = PhysicsInformedNN(chain, _strategy, init_params=initθ)

@named pde_system = PDESystem(eqs,bcs, domains, [r,t], [c_sp(r,t), c_sn(r,t)])
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

# phi = discretization.phi