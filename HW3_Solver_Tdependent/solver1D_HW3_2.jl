using PyPlot

#-------- Input Parameters --------#

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

# println("
# a $a
# b $b
# c $c
# d $d
# P $P
# Q $Q")

    return PHI
end


function runsolver_warmstart(N_cv;
    h = 10,
    Tinf = 273,
    Tb = 400,
    boundary = 1E20,
    rod_L = 0.02,
    D = 0.003,
    epsilon = 0.0,
    theta = 5.67E-8,
    k_in = 401)

    aw = zeros(N_cv+2)
    ap = zeros(aw)
    ae = zeros(aw)
    b = zeros(aw)
    k = zeros(aw)
    kw_h = zeros(aw) #harmonic mean
    ke_h = zeros(aw)
    Sc = zeros(N_cv) #constant source term
    Sp = zeros(N_cv) #slope of source term, multiplied by temp


    Sp_spec = -4/D*h
    Sc_spec = 4/D*h*Tinf


    cv_x = collect(linspace(0,rod_L,N_cv+1))
    cv_dx = cv_x[2:end] - cv_x[1:end-1]
    node_x = [0;(cv_x[2:end]+cv_x[1:end-1])/2;cv_x[end]]
    node_dx = node_x[2:end] - node_x[1:end-1]

    for i = 1:length(k)
        if node_x[i]<boundary
            k[i] = 401
        else
            # k[i] = 137*exp(25*node_x[i]-2)
        end
    end

    for i = 1:length(Sc)
        if cv_x[i]<boundary
            Sp[i] =Sp_spec
            Sc[i] =Sc_spec
        else
            # Sc[i] = 0
            # Sp[i] = 0
        end
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

    # tic()
    # T = A\b
    # toc()
    # tic()
    T = TDMA(aw,ap,ae,b)
    # toc()
    # println("A")
    # for i = 1:length(A[:,1])
    #     println(A[i,:])
    # end
    #
    # println("
    # b: $b")
    #
    # println("
    # T: $T
    # TestT: $TestT")

    Sstar = Sc + Sp.*T[2:end-1]
    # println("1: $Sstar")
    return T,node_x,Sstar
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
    k_in = 401)

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

    T = TDMA(aw,ap,ae,b)
    S = Sc + Sp.*T[2:end-1]
    # S = Sstar+Sp.*(T[2:end-1]-Tstar[2:end-1])
    # println("2: $S")
    return T,node_x,S
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
epsilon = 0.0
theta = 5.67E-8
k_in = 401

T_e0,node_x,S_e0 = runsolver_warmstart(N_cv;
    h = h,
    Tinf = Tinf,
    Tb = Tb,
    boundary = boundary,
    rod_L = rod_L,
    D = D,
    epsilon = epsilon,
    theta = theta,
    k_in = k_in)

T = T_e0#zeros(T_e0)
Sstar = zeros(N_cv)
# Run until convergence
Tstar = zeros(T)
i = 1

while true

    T,node_x,Sstar = runsolver_iterative(N_cv,Tstar,Sstar;
    h = h,
    Tinf = Tinf,
    Tb = Tb,
    boundary = boundary,
    rod_L = rod_L,
    D = D,
    epsilon = epsilon,
    theta = theta,
    k_in = k_in)

    diff = maximum(abs.(T-Tstar))

    Tstar = copy(T)

    # println("i: $i diff $diff T: $Tstar")
    i+=1
    if i>5 && diff<1E-4
        break
    end
end

analy_x = linspace(0,rod_L,100)
n = sqrt(4*h/(D*k_in))
Tanalyitical = Tinf+(Tb-Tinf)*(cosh.(n*(rod_L-analy_x)))/(cosh.(n*rod_L))

PyPlot.close("all")
figname = "HW3_2"
PyPlot.figure(figname)
PyPlot.plot(analy_x,Tanalyitical,label = L"Analytical $\epsilon$ = 0.0")
PyPlot.plot(node_x,T_e0,label = L"TDMA $\epsilon$ = 0.0")
PyPlot.plot(node_x,T,label = L"Iterative $T_0$=zeros, $\epsilon$ = 0.0")
PyPlot.xlabel("x")
PyPlot.ylabel("T (K)")
PyPlot.legend(loc = "best")
PyPlot.savefig(figname,transparent = true)


epsilon = 1.0
Sstar = S_e0
# Run until convergence
Tstar = T_e0
while true

    T,node_x,Sstar = runsolver_iterative(N_cv,Tstar,Sstar;
    h = h,
    Tinf = Tinf,
    Tb = Tb,
    boundary = boundary,
    rod_L = rod_L,
    D = D,
    epsilon = epsilon,
    theta = theta,
    k_in = k_in)

    diff = maximum(abs.(T-Tstar))

    Tstar = copy(T)

    # println("i: $i diff $diff T: $Tstar")
    i+=1
    if i>5 && diff<1E-4
        break
    end
end

figname = "HW3_3"
PyPlot.figure(figname)
PyPlot.plot(analy_x,Tanalyitical,label = L"Analytical $\epsilon$ = 0.0")
PyPlot.plot(node_x,T,label = L"Iterative $\epsilon$ = 1.0")
PyPlot.xlabel("x")
PyPlot.ylabel("T (K)")
PyPlot.legend(loc = "best")
PyPlot.savefig(figname,transparent = true)
