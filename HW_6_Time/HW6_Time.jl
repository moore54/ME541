using PyPlot, Dierckx


function gauss_sidel(aw,ap,ae,bstar,Tstar,relax)

    T = zeros(Tstar)
    ERR = 1E20
    j=0
    while ERR>1E-6
        i = 1
        T[1] = (0+ae[i]*Tstar[i+1]+bstar[i])/ap[i]
        for i = 2:length(Tstar)-1
            T[i] = (aw[i]*T[i-1]+ae[i]*Tstar[i+1]+bstar[i])/ap[i]
        end
        T[end] = (aw[end]*T[end-1]+0+bstar[end])/ap[end]

        #Apply Relaxation
        T = Tstar+relax*(T-Tstar)
        ERR = maximum(abs.(T-Tstar))
        Tstar = copy(T)
        j+=1

    end

    return T
end

function converge_one_timestep(N_node,Tstar,Told;ii=1,
    scheme = "central",
    h_conv = 1.0,
    u_vel = 10,
    Tinf = 273,
    Tsur = 293,
    Tb = 400,
    rho = 1,
    C_heat = 875,
    rod_L = 0.02,
    b_thick = 0.0003,
    gamma_conv = 1,
    epsilon = 1.0,
    theta = 5.67E-8,
    k_in = 401,
    TOL = 1e-6,
    relax = 1.0,
    delta_t = 0.01)

    T = []
    node_x = []
    qw = []
    qe = []
    dxdt = []

    qestar = ones(length(Tstar)-1)*1e20
    qwstar = ones(length(Tstar)-1)*1e20

    while true

        T,node_x,qw,qe,dxdt = inner_iteration(N_node,Tstar,Told;
        scheme = scheme,
        h_conv = h_conv,
        u_vel = u_vel,
        Tinf = Tinf,
        Tsur = Tsur,
        Tb = Tb,
        rho = rho,
        C_heat = C_heat,
        rod_L = rod_L,
        b_thick = b_thick,
        gamma_conv = gamma_conv,
        epsilon = epsilon,
        theta = theta,
        k_in = k_in,
        relax = relax,
        delta_t = delta_t)

        # diff = maximum(abs.(T-Tstar)) #max temperature error
        # diff = (abs.(T[end]-Tstar[end])) #tip temperature
        # diff = (abs.(T[1]-Tstar[1])) #base temperature
        # diff = (abs.(qe[end-1]-qestar[end-1])) #tip heat flux
        diff = (abs.(qw[1]-qwstar[1])) #base heat flux

        Tstar = copy(T)
        qwstar = copy(qw)
        qestar = copy(qe)

        # println("i: $i diff $diff T: $Tstar")
        ii+=1
        # println(ii)
        if ii>1 && diff<TOL
            break
        end

        if false#(ii%100==0)
            PyPlot.plot(node_x,T,".")
            PyPlot.pause(0.1)
        end
    end
    return T,node_x,ii,qw[1],dxdt
end

function TDMA(aw,ap,ae,b_in)
    a = ap
    b = ae
    c = aw
    d = b_in

    PHI = zeros(b)
    P = zeros(b)
    Q = zeros(b)

    #1 find P[1]
    P[1] = b[1]/a[1]
    Q[1] = d[1]/a[1]

    #2-3 go through all interior
    for i = 2:length(b)
        P[i] = b[i]/(a[i]-c[i]*P[i-1])
        Q[i] = (d[i]+c[i]*Q[i-1])/(a[i]-c[i]*P[i-1])
    end
    #4 last point
    PHI[end] = Q[end]

    #5 Back sub to find PHIs
    for i = reverse(2:length(b))
        PHI[i-1] = P[i-1]*PHI[i]+Q[i-1]
    end
    return PHI
end

function A_fun(pec,scheme = "central")

    if scheme == "central"
        res = 1-0.5*abs(pec)
    elseif scheme == "upwind"
        res = 1
    elseif scheme == "hybrid"
        res = max(0,1-0.5*abs(pec))
    elseif scheme == "power"
        res = max(0,(1-0.1*abs(pec))^5)
    end

    return res
end

function inner_iteration(N_node,Tstar,Told;
    scheme = "central",
    h_conv = 1.0,
    u_vel = 10,
    Tinf = 273,
    Tsur = 293,
    Tb = 400,
    rho = 1,
    C_heat = 875,
    rod_L = 0.02,
    b_thick = 0.0003,
    gamma_conv = 1,
    epsilon = 1.0,
    theta = 5.67E-8,
    k_in = 401,
    relax = 1.0,
    delta_t = 0.01)

    aw = zeros(N_node)
    ap = zeros(aw)
    # ap0 = zeros(aw)
    ae = zeros(aw)
    b = zeros(aw)
    k = zeros(aw)
    kw_h = zeros(aw) #harmonic mean
    ke_h = zeros(aw)
    Sc = zeros(N_node) #constant source term
    Sp = zeros(N_node) #slope of source term, multiplied by temp

    node_x = collect(linspace(0,rod_L,N_node))
    node_dx = node_x[2:end] - node_x[1:end-1]
    cv_x = (node_x[2:end]+node_x[1:end-1])/2
    cv_dx = [cv_x[1]-node_x[1];cv_x[2:end] - cv_x[1:end-1];node_x[end]-cv_x[end]]

    for i = 1:length(k)
        k[i] = k_in
    end

    for i = 1:N_node
        Sp[i] =-2/b_thick*(h_conv+4*epsilon*theta*Tstar[i]^3)
        Sc[i] =2/b_thick*(h_conv*Tinf+epsilon*theta*(3*Tstar[i]^4+Tsur^4))
    end

    #-------- Apply harmonic mean --------#
    node_dx_minus = cv_x-node_x[1:end-1]
    node_dx_plus = node_x[2:end]-cv_x

    for i = 1:length(kw_h)-1
        ke_h[i] = node_dx[i]*k[i+1]*k[i]/(k[i+1]*node_dx_minus[i]+k[i]*node_dx_plus[i])
        kw_h[i+1] = node_dx[i]*k[i]*k[i+1]/(k[i]*node_dx_plus[i]+k[i+1]*node_dx_minus[i])
    end

    #-------- Assemble the coefficients --------#
    # Apply Left BC
    aw[1] = 0
    ae[1] = 0
    ap[1] = 1
    b[1] = Tb

    # Apply Interior Points
    for i = 2:N_node-1
        #note D indexing
        ap0 = rho*C_heat*cv_dx[i]/delta_t #TODO: check if non-uniform spacing
        aw[i] = kw_h[i]/node_dx[i-1]
        ae[i] = ke_h[i]/node_dx[i]
        ap[i] = aw[i]+ae[i]+ap0-Sp[i-1]*cv_dx[i-1]
        b[i] = Sc[i-1]*cv_dx[i-1]+ap0*Told[i]
    end

    # Apply Right BC
    ap0 = rho*C_heat*cv_dx[end]/delta_t
    aw[end] = kw_h[end]/node_dx[end-1]
    ae[end] = 0
    ap[end] = aw[end]+ap0-Sp[end-1]*cv_dx[end-1]
    b[end] = Sc[end-1]*cv_dx[end-1]+ap0*Told[end]

    #-------- Solve the System --------#
    T = TDMA(aw,ap,ae,b)
    # T = gauss_sidel(aw,ap,ae,b,Tstar,relax)

    #------- Calculate Heat Flux -------#

    qw = -kw_h[2:end].*(T[2:end]-T[1:end-1])./node_dx - (Sc[2:end]+Sp[2:end].*T[2:end]).*cv_dx[2:end]
    qe = -ke_h[1:end-1].*(T[2:end]-T[1:end-1])./node_dx # - (Sc+Sp*Tp)*cv_dx

    dxdt = (cv_dx[end]/delta_t)


    return T,node_x,qw,qe,dxdt
end

rc("figure", figsize=(4.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.92)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
color_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7", "#000486","#700000","#006907","#4C0099"]

#Initialize solution
h_conv = [5.0,30.0]
Tinf = 293
Tsur = 293
T0 = 293
Tb_real = 353
rod_L = 0.01
b_thick = 0.0003
epsilon = [0.0, 1.0]
theta = 5.67E-8
k_in = 177
rho = 2770
C_heat = 875
relax = 1.2
delta_t = logspace(-4,-6,5)
gamma_conv = 1.0
scheme = "upwind"
print_time = [0.001,0.5,1.0,2.0,4.0,8.0,16.0,32.0,64.0]

a_cross = 1.0*b_thick
perim = 2*1.0+2*b_thick



# N_node = zeros(4)
# delta_x = zeros(N_node)
# node0 = 10.0
# for i = 1:length(N_node)
#     N_node[i] = node0*2
#     node_x = collect(linspace(0.0,rod_L,N_node[i]))
#     node_dx = node_x[2:end] - node_x[1:end-1]
#     delta_x[i] = node_dx[2]
# end
#
# N_node = round.(Int,N_node)


iters = []
PyPlot.close("all")
figname = "HW6_2_grid_conv"
Tsave = []
node_x = []
names = []

# Do grid convergence on h 30 e 1

delta_t = 0.008
N_node = 10
N_node_old = N_node
T = ones(N_node)*T0 # warmstart by placing this here

#Percentage change/ tolerance
tol_d_loop = 1E-2
tol_N_node = 1E-1
tol_delta_t = 1E-6
tol_time = 1E-6 #Absolute change

dxdtOgoal = 1

# tic()
# qw_save = []
# N_node_save = []
# delta_t_save = []
# qw = []
# dxdt = []
# d_loop = 1E20
# qw_loop_old = 1E20
# # while d_loop>tol_d_loop #keep looping
# println("rerun \n \n")
# d_N_node = 1E20
# qw_N_node_old = 1E20
# # while d_N_node>tol_N_node # N_node
# for jjj = 1:8
#
#     d_delta_t = 1E20
#     qw_delta_t_old = 1E20
#     # while d_delta_t>tol_delta_t # delta_t
#     for iii = logspace(5,-3,10)
#         delta_t = dxdtOgoal*rod_L/N_node*iii
#         time = 0
#         time_diff = 1E20
#         i_print = 1
#         Tinterp = Dierckx.Spline1D(linspace(0,rod_L,N_node_old),T)
#         T = Tinterp(linspace(0,rod_L,N_node))
#         while time_diff > tol_time
#
#             time+=delta_t
#
#             Tstar = copy(T) #break the reference
#             Told = copy(T)
#             T,node_x,iters,qw,dxdt = converge_one_timestep(N_node,Tstar,Told;
#             scheme = scheme,
#             h_conv = h_conv[2],
#             u_vel = h_conv[2],
#             Tinf = Tinf,
#             Tsur = Tsur,
#             Tb = Tb_real,
#             rho = rho,
#             C_heat = C_heat,
#             rod_L = rod_L,
#             b_thick = b_thick,
#             gamma_conv = gamma_conv,
#             epsilon = epsilon[2],
#             theta = theta,
#             k_in = k_in,
#             relax = relax,
#             delta_t = delta_t)
#
#             push!(Tsave,T)
#
#
#             # PyPlot.savefig("$(figname)_$(h_conv[j]).pdf",transparent = true)
#             time_diff = sum(abs.(Told-T))
#
#         end
#
#         push!(qw_save,qw)
#         push!(N_node_save,N_node)
#         push!(delta_t_save,delta_t)
#
#         N_node_old = N_node
#         d_delta_t = abs.((qw-qw_delta_t_old)/qw)*100
#         qw_delta_t_old = qw
#
#         # println(qw)
#         println(dxdt)
#
#         # if d_delta_t > tol_delta_t
#         #     delta_t = delta_t/1.25
#         #     # println("delta_t $delta_t")
#         # else
#         # println("move on")
#     end
#
#
#
#     d_N_node = abs.((qw-qw_N_node_old)/qw)*100
#     qw_N_node_old = qw
#     # if d_N_node > tol_N_node
#     N_node = N_node*2
#     N_node = round(Int,N_node)
#     # println("Node $N_node \n")
#     # else
#     # println("move on")
#     # end
#     println()
# end
#
# d_loop = abs.((qw-qw_loop_old)/qw)*100
# qw_loop_old = qw
# # end
# toc()
#
# figname = "qw_base"
# PyPlot.figure("qw_base")
# PyPlot.plot3D(delta_t_save,N_node_save,round.(qw_save/1e6,3),".")
# PyPlot.xlabel("delta t (s)")
# PyPlot.ylabel("Nodes (#)")
# PyPlot.zlabel("q\" (1E6)")
# PyPlot.savefig(figname,transparent = true)
#
# figname = "qw_base2D_nodes"
# PyPlot.figure("qw_base2D_nodes")
# PyPlot.plot(N_node_save,qw_save/1e6,".")
# PyPlot.xlabel("Nodes (#)")
# PyPlot.ylabel("q\" (1E6)")
# PyPlot.savefig(figname,transparent = true)
#
#
# figname = "qw_base2D_delta_t"
# PyPlot.figure("qw_base2D_delta_t")
# PyPlot.semilogx(delta_t_save,qw_save/1e6,".")
# PyPlot.xlabel("delta t (s)")
# PyPlot.ylabel("q\" (1E6)")
# PyPlot.savefig(figname,transparent = true)


N_node = 320
dx = rod_L/N_node
delta_t = dxdtOgoal*rod_L/N_node

# PyPlot.zlim(minimum(qw_save/1e6), minimum(qw_save/1e6)*1.01)

# PyPlot.legend(loc = "best")

Tsave = []
node_x = []
names = []

for k = 1:length(epsilon)
    for j = 1:length(h_conv)

        analy_x = linspace(0,rod_L,N_node)
        mcoeff = sqrt(h_conv[j]*perim/(k_in*a_cross))
        Tanalyitical = cosh.(mcoeff*(rod_L-analy_x))/cosh.(mcoeff*rod_L)*(Tb_real-Tinf)+Tinf
        PyPlot.figure("$(figname)_$(h_conv[j])_$(epsilon[k])")
        PyPlot.plot(analy_x,Tanalyitical,color = color_cycle[j],label = "Analytical $(h_conv[j]) m/s")
        time = 0
        Esum = []
        qw_save = []
        T = ones(N_node)*T0

        # PyPlot.plot(linspace(0,rod_L,length(T)),T,".",label = "t 0.0-")
        max_diff = 1E20
        i_print = 1
        while max_diff > 1e-6
        #for i = 1:length(time)
        time+=delta_t

            Tstar = copy(T) #break the reference
            Told = copy(T)
            T,node_x,iters,qw = converge_one_timestep(N_node,Tstar,Told;
            scheme = scheme,
            h_conv = h_conv[j],
            u_vel = h_conv[j],
            Tinf = Tinf,
            Tsur = Tsur,
            Tb = Tb_real,
            rho = rho,
            C_heat = C_heat,
            rod_L = rod_L,
            b_thick = b_thick,
            gamma_conv = gamma_conv,
            epsilon = epsilon[k],
            theta = theta,
            k_in = k_in,
            relax = relax,
            delta_t = delta_t)

            push!(Tsave,T)
            push!(qw_save,qw)

            # PyPlot.savefig("$(figname)_$(h_conv[j]).pdf",transparent = true)
            max_diff = maximum(abs.(Told-T))
            push!(Esum, maximum(abs.(Tanalyitical-T)))

            if max_diff < 1e-6#time%print_time[i_print]<delta_t && time%print_time[i_print]>-delta_t
                i_print+=1
                PyPlot.plot(node_x,T,".",label = "t $(round(time,3))")
                PyPlot.xlabel("Y-location")
                PyPlot.ylabel("T")
                PyPlot.legend(loc = "best")
            end


        end

        rc("figure", figsize=(6.5, 2.6))
        rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.72)
        PyPlot.plot(node_x,T,".",label = "t $(round(time,3))")
        PyPlot.xlabel("Y-location")
        PyPlot.ylabel("T")
        PyPlot.legend(loc="lower left", bbox_to_anchor=(1, 0.5))

        rc("figure", figsize=(4.5, 2.6))
        rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.92)
        PyPlot.figure("error")
        PyPlot.semilogy(0:delta_t:time,Esum,"-",label = "$(h_conv[j])_$(epsilon[k])")
        PyPlot.xlabel("t (s)")
        PyPlot.ylabel("Maximum Steady State Error")
        PyPlot.legend(loc = "best")

        PyPlot.figure("qw_base")
        PyPlot.semilogy(0:delta_t:time,qw_save,"-",label = "$(h_conv[j])_$(epsilon[k])")
        PyPlot.xlabel("t (s)")
        PyPlot.ylabel("Base Heat Transfer (q\")")
        PyPlot.legend(loc = "best")

        # push!(names,"V_$(h_conv[j])m/s_$delta_x[i]")
    end
end
