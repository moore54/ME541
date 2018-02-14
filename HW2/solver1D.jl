using PyPlot

#-------- Input Parameters --------#


function runsolver(N_cv)
    aw = zeros(N_cv+2)
    ap = zeros(aw)
    ae = zeros(aw)
    b = zeros(aw)
    k = zeros(aw)
    kw_h = zeros(aw) #harmonic mean
    ke_h = zeros(aw)
    Sc = zeros(N_cv) #constant source term
    Sp = zeros(N_cv) #slope of source term, multiplied by temp


    h = 1000
    Tinf = 300
    boundary = 0.03

    cv_x = collect(linspace(0,0.05,N_cv+1))
    cv_dx = cv_x[2:end] - cv_x[1:end-1]
    node_x = [0;(cv_x[2:end]+cv_x[1:end-1])/2;cv_x[end]]
    node_dx = node_x[2:end] - node_x[1:end-1]

    for i = 1:length(k)
        if node_x[i]<boundary
            k[i] = 15
        else
            k[i] = 60
        end
    end

    for i = 1:length(Sc)
        if cv_x[i]<boundary
            Sc[i] = 4E6
            Sp[i] = 0
        else
            Sc[i] = 0
            Sp[i] = 0
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
    ae[1] = ke_h[1]/node_dx[1]
    ap[1] = ke_h[1]/node_dx[1]
    b[1] = 0

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
    ap[end] = kw_h[end]/node_dx[end]+h
    b[end] = h*Tinf

    #-------- Assemble the A_matrix --------#

    A = zeros(N_cv+2,N_cv+2)

    # First Row
    A[1,1] = ap[1]
    A[1,2] = -ae[1]

    # Interior Rows
    for i = 2:N_cv+1 #loop through the interior rows
        A[i,i-1] = -aw[i]
        A[i,i] = ap[i]
        A[i,i+1] = -ae[i]
    end

    # Last Row
    A[end,end-1] = -aw[end]
    A[end,end] = ap[end]

    #-------- Solve the System --------#

    T = A\b
    return T,node_x,A,b
end

rc("figure", figsize=(4.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.72)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])

PyPlot.figure("3")
N_cv = [5,20]
A = []
b = []
for i = 1:length(N_cv)
    T,node_x,A,b = runsolver(N_cv[i])
    PyPlot.plot(node_x,T,label = "$(N_cv[i]) Volumes")
end
PyPlot.xlabel("x")
PyPlot.ylabel("T (K)")
PyPlot.legend(loc = "best")
