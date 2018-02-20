using PyPlot


function gauss_sidel(aw,ap,ae,bstar,Tstar,relax)

    T = zeros(Tstar)
    i = 1
    T[1] = (0+ae[i]*Tstar[i+1]+bstar[i])/ap[i]
    for i = 2:length(Tstar)-1
        T[i] = (aw[i]*T[i-1]+ae[i]*Tstar[i+1]+bstar[i])/ap[i]
    end
    T[end] = (aw[end]*T[end-1]+0+bstar[end])/ap[end]

    #Apply Relaxation
    T = Tstar+relax*(T-Tstar)

    # println(Tstar)
    return T
end

function itersolve(N_cv,Tstar,Sstar;ii=1,
    h = 10,
    Tinf = 273,
    Tb = 400,
    boundary = 1E20,
    rod_L = 0.02,
    D = 0.003,
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

        T,node_x,Sstar,qw,qe = runsolver_iterative(N_cv,Tstar,Sstar;
        h = h,
        Tinf = Tinf,
        Tb = Tb,
        boundary = boundary,
        rod_L = rod_L,
        D = D,
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
            # println("$T")
        end
    end
    return T,node_x,Sstar,ii
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
        # println(i)
        PHI[i-1] = P[i-1]*PHI[i]+Q[i-1]
    end
    return PHI
end


function runsolver_iterative(N_cv,Tstar,Sstar;
    h = 10,
    Tinf = 273,
    Tb = 400,
    boundary = 1E20,
    rod_L = 0.02,
    D = 0.003,
    epsilon = 1.0,
    theta = 5.67E-8,
    k_in = 401,
    relax = 1.0)

    aw = zeros(N_cv+2)
    ap = zeros(aw)
    ae = zeros(aw)
    b = zeros(aw)
    k = zeros(aw)
    kw_h = zeros(aw) #harmonic mean
    ke_h = zeros(aw)
    Sc = zeros(N_cv) #constant source term
    Sp = zeros(N_cv) #slope of source term, multiplied by temp

    cv_x = collect(linspace(0,rod_L,N_cv+1))
    cv_dx = cv_x[2:end] - cv_x[1:end-1]
    node_x = [0;(cv_x[2:end]+cv_x[1:end-1])/2;cv_x[end]]
    node_dx = node_x[2:end] - node_x[1:end-1]

    for i = 1:length(k)
        if node_x[i]<boundary
            k[i] = k_in
        else
            # k[i] = 137*exp(25*node_x[i]-2)
        end
    end

    for i = 1:N_cv

        Sp[i] =-4/D*h-16*epsilon*theta/D*Tstar[i+1]^3
        Sc[i] =4/D*h*Tinf#Sstar[i]-Sp[i]*Tstar[i+1]

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
    for i = 2:N_cv+1
        aw[i] = kw_h[i]/node_dx[i-1]
        ae[i] = ke_h[i]/node_dx[i]
        ap[i] = aw[i]+ae[i]-Sp[i-1]*cv_dx[i-1]
        b[i] = Sc[i-1]*cv_dx[i-1]
    end

    # Apply Right BC
    aw[end] = kw_h[end]/node_dx[end]
    ae[end] = 0
    ap[end] = kw_h[end]/node_dx[end]
    b[end] = 0

    #-------- Solve the System --------#
    # T = TDMA(aw,ap,ae,b)
    T = gauss_sidel(aw,ap,ae,b,Tstar,relax)
    S = Sc + Sp.*T[2:end-1]

    #------- Calculate Heat Flux -------#

    qw = -kw_h[2:end].*(T[2:end]-T[1:end-1])./node_dx
    qe = -ke_h[1:end-1].*(T[2:end]-T[1:end-1])./node_dx

    return T,node_x,S,qw,qe
end

rc("figure", figsize=(4.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.92)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])


#Initialize solution
N_cv = 100

h = 10
Tinf = 273
Tb = 400
boundary = 1E20
rod_L = 0.02
D = 0.003
epsilon = 1.0
theta = 5.67E-8
k_in = 401
T0 = 400
relax = 1.0

analy_x = linspace(0,rod_L,100)
n = sqrt(4*h/(D*k_in))
Tanalyitical = Tinf+(Tb-Tinf)*(cosh.(n*(rod_L-analy_x)))/(cosh.(n*rod_L))


N_cv = zeros(4)
N_cv[1] = 10
for i = 2:length(N_cv)
    N_cv[i] = N_cv[i-1]*2
end
N_cv = round.(Int,N_cv)
iters_array = zeros(length(N_cv))
PyPlot.close("all")
figname = "HW3_2_grid_conv"
Tsave = []
PyPlot.figure(figname)
for j = 1:length(N_cv)
    analy_x = linspace(0,rod_L,N_cv[j]+2)
    n = sqrt(4*h/(D*k_in))
    T = Tinf+(Tb-Tinf)*(cosh.(n*(rod_L-analy_x)))/(cosh.(n*rod_L))
    # T = ones(N_cv[j]+2)*T0
    Sstar = zeros(N_cv[j])
    Tstar = copy(T)
    T,node_x,Sstar,iters = itersolve(N_cv[j],Tstar,Sstar;
    h = h,
    Tinf = Tinf,
    Tb = Tb,
    boundary = boundary,
    rod_L = rod_L,
    D = D,
    epsilon = epsilon,
    theta = theta,
    k_in = k_in,
    relax = relax)

    push!(Tsave,T)


    PyPlot.plot(node_x,T,label = "$(N_cv[j]) CV")

    iters_array[j] = iters
end
alay_x = linspace(0,rod_L,100)
Tanalyitical = Tinf+(Tb-Tinf)*(cosh.(n*(rod_L-alay_x)))/(cosh.(n*rod_L))
PyPlot.plot(alay_x,Tanalyitical,label = "Analytical")
PyPlot.xlabel("Y-location")
PyPlot.ylabel("T (K)")
PyPlot.legend(loc = "best")
PyPlot.savefig(figname,transparent = true)

PyPlot.figure()
for i = 1:length(N_cv)
PyPlot.plot(N_cv[i],Tsave[i][end],".")
end
PyPlot.plot(100,Tanalyitical[end],".")
