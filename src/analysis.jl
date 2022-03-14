include("ode_spm.jl")
@time include("pinn_spm.jl")

# Create plots for PINN predictions
c_sp_plt = plot(c_sp_predict, fmt = :png, title="c_sp predict", label = "c_sp")
c_sn_plt = plot(c_sn_predict, fmt = :png, title="c_sn predict", label = "c_sn")

savefig(c_sp_plt,"../plots/c_sp_predict.png")
savefig(c_sn_plt,"../plots/c_sn_predict.png")

# println(c_sp_predict)
# println(c_sn_predict)
# println(sol)

