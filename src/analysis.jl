using Plots

include("ode_spm.jl")
@time include("pinn_spm.jl")

# Create plots for PINN predictions
# Posivive and Negative Concentrations vs dt 

plt_c_sp = plot(title = "Positive Concentration vs. dt", xlabel = "dt", ylabel = "Concentration")
for i in (1:11)
    lab = string("dr= ", string(i/10))
    plot!(x_axis, pred_c_sp_dt[i], label = lab)
end

plt_c_sn = plot(title = "Negative Concentration vs. dt", xlabel = "dt", ylabel = "Concentration")
for i in (1:11)
    lab = string("dr= ", string(i/10))
    plot!(x_axis, pred_c_sn_dt[i], label = lab)
end

savefig(plt_c_sp,"../plots/c_sp_predict.png")
savefig(plt_c_sn,"../plots/c_sn_predict.png")