using PyPlot


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

function itersolve(N_node,Tstar;ii=1,
    scheme = "central",
    u_vel = 10,
    Tinf = 273,
    Tb = 400,
    rho = 1,
    rod_L = 0.02,
    gamma = 1,
    epsilon = 1.0,
    theta = 5.67E-8,
    k_in = 401,
    TOL = 1e-6,
    relax = 1.0)
    T = []
    node_x = []
    qw = []
    qe = []
    qestar = ones(length(Tstar)-1)*1e20
    qwstar = ones(length(Tstar)-1)*1e20
    while true

        T,node_x,qw,qe = runsolver_iterative(N_node,Tstar;
        scheme = scheme,
        u_vel = u_vel,
        Tinf = Tinf,
        Tb = Tb,
        rho = rho,
        rod_L = rod_L,
        gamma = gamma,
        epsilon = epsilon,
        theta = theta,
        k_in = k_in,
        relax = relax)

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
    return T,node_x,ii
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


function runsolver_iterative(N_node,Tstar;
    scheme = "central",
    u_vel = 10,
    Tinf = 273,
    Tb = 400,
    rho = 1,
    rod_L = 0.02,
    gamma = 1,
    epsilon = 1.0,
    theta = 5.67E-8,
    k_in = 401,
    relax = 1.0)

    aw = zeros(N_node)
    ap = zeros(aw)
    ae = zeros(aw)
    b = zeros(aw)
    k = zeros(aw)
    kw_h = zeros(aw) #harmonic mean
    ke_h = zeros(aw)
    Sc = zeros(N_node) #constant source term
    Sp = zeros(N_node) #slope of source term, multiplied by temp

    # cv_x = collect(linspace(0,rod_L,N_node+1))
    # cv_dx = cv_x[2:end] - cv_x[1:end-1]
    # node_x = [0;(cv_x[2:end]+cv_x[1:end-1])/2;cv_x[end]]
    # node_dx = node_x[2:end] - node_x[1:end-1]

    # node_x = collect(linspace(rod_L/(N_node*2),rod_L-rod_L/(N_node*2),N_node))
    # node_dx = node_x[2:end] - node_x[1:end-1]
    # cv_x = [0;(node_x[2:end]+node_x[1:end-1])/2;rod_L]
    # cv_dx = cv_x[2:end] - cv_x[1:end-1]

    node_x = collect(linspace(0,rod_L,N_node))
    node_dx = node_x[2:end] - node_x[1:end-1]
    cv_x = (node_x[2:end]+node_x[1:end-1])/2
    cv_dx = [cv_x[1]-node_x[1];cv_x[2:end] - cv_x[1:end-1];node_x[end]-cv_x[end]]

    for i = 1:length(k)
        k[i] = k_in
        # if node_x[i]<boundary
            # k[i] = k_in
        # else
        #     # k[i] = 137*exp(25*node_x[i]-2)
        # end
    end

    for i = 1:N_node

        Sp[i] =0#-4/D*h-16*epsilon*theta/D*Tstar[i]^3
        Sc[i] =0#4/D*h*Tinf+12*epsilon*theta/D*Tstar[i]^4+4*epsilon*theta/D*Tinf^4

    end

    #-------- Apply harmonic mean --------#
    node_dx_minus = cv_x-node_x[1:end-1]
    node_dx_plus = node_x[2:end]-cv_x

    for i = 1:length(kw_h)-1
        ke_h[i] = node_dx[i]*k[i+1]*k[i]/(k[i+1]*node_dx_minus[i]+k[i]*node_dx_plus[i])
        kw_h[i+1] = node_dx[i]*k[i]*k[i+1]/(k[i]*node_dx_plus[i]+k[i+1]*node_dx_minus[i])
    end

    F = rho*u_vel #TODO will likely change
    D = gamma./node_dx

    pec = F./D
    #-------- Assemble the coefficients --------#
    # Apply Left BC
    aw[1] = 0
    ae[1] = 0
    ap[1] = 1
    b[1] = 1

    # Apply Interior Points
    for i = 2:N_node-1
        #note D indexing
        aw[i] = D[i-1]*A_fun(pec[i-1],scheme)+max(F,0)#kw_h[i]/node_dx[i-1]
        ae[i] = D[i+0]*A_fun(pec[i],scheme)+max(-F,0)#ke_h[i]/node_dx[i]
        ap[i] = aw[i]+ae[i]+(F-F)-Sp[i-1]*cv_dx[i-1]
        b[i] = Sc[i-1]*cv_dx[i-1]
    end

    # Apply Right BC
    aw[end] = 0
    ae[end] = 0
    ap[end] = 1
    b[end] = 0

    #-------- Solve the System --------#
    T = TDMA(aw,ap,ae,b)
    # T = gauss_sidel(aw,ap,ae,b,Tstar,relax)

    #------- Calculate Heat Flux -------#

    qw = -kw_h[2:end].*(T[2:end]-T[1:end-1])./node_dx
    qe = -ke_h[1:end-1].*(T[2:end]-T[1:end-1])./node_dx

    return T,node_x,qw,qe
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
h = 10
Tinf = 273
Tb = 0#400
boundary = 1E20
rod_L = 1.0
D = 0.003
epsilon = 0.0
theta = 5.67E-8
k_in = 401
T0 = 0400
relax = 1.2
rho = 1.0
gamma = 1.0
scheme = "upwind"
scheme_a = ["central","upwind", "hybrid", "power"]


N_node = zeros(4)
N_node = 21
u_vela = [1.0,20.0,75.0]
iters_array = zeros(length(N_node))
PyPlot.close("all")
figname = "HW3_2_grid_conv"
Tsave = []
node_x = []
names = []
Esum = []
for j = 1:length(u_vela)
    PyPlot.figure("$(figname)_$(u_vela[j])")
    analy_x = linspace(0,rod_L,21)
    Pec = rho*u_vela[j]*rod_L/gamma
    Tanalyitical = 1-(exp.(Pec*analy_x/rod_L)-1)/(exp.(Pec)-1)
    PyPlot.plot(analy_x,Tanalyitical,color = color_cycle[j],label = "Analytical $(u_vela[j]) m/s")
    for i = 1:length(scheme_a)
    scheme = scheme_a[i]

    # j = 1

        T = ones(N_node+2)*T0

        Tstar = copy(T)
        T,node_x,iters = itersolve(N_node,Tstar;
        scheme = scheme,
        u_vel = u_vela[j],
        Tinf = Tinf,
        Tb = Tb,
        rho = rho,
        rod_L = rod_L,
        gamma = gamma,
        epsilon = epsilon,
        theta = theta,
        k_in = k_in,
        relax = relax)

        push!(Tsave,T)
        push!(names,"V_$(u_vela[j])m/s_$scheme")
        # printout = [node_x,T]
        # println("
        # $(u_vela[j]) m/s, $(scheme_a[i])
        # $(printout)
        #
        # ")

        PyPlot.plot(node_x,T,".",color = color_cycle[i],label = "$(scheme_a[i])")
        PyPlot.xlabel("Y-location")
        PyPlot.ylabel(L"\phi")
        PyPlot.legend(loc = "best")
        PyPlot.savefig("$(figname)_$(u_vela[j]).pdf",transparent = true)

        push!(Esum, sum(abs.(Tanalyitical-T)))
    end
end


nameout = "./hw5.txt"

write(nameout,"x\t$(names[1])\t$(names[2])\t$(names[3])\t$(names[4])\t$(names[5])\t$(names[6])\t$(names[7])\t$(names[8])\t$(names[9])\t$(names[10])\t$(names[11])\t$(names[12])\n")

open(nameout,"a") do x
    for i = 1:length(node_x)
        write(x,"$(node_x[i])\t$(Tsave[1][i])\t$(Tsave[2][i])\t$(Tsave[3][i])\t$(Tsave[4][i])\t$(Tsave[5][i])\t$(Tsave[6][i])\t$(Tsave[7][i])\t$(Tsave[8][i])\t$(Tsave[9][i])\t$(Tsave[10][i])\t$(Tsave[11][i])\t$(Tsave[12][i])\n")
    end
end

# PyPlot.savefig("$figname.pdf",transparent = true)
