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


    h = 10
    Tinf = 273
    Tb = 400
    boundary = 1E20
    rod_L = 0.02
    D = 0.003
    Sc_spec = 4/D*h*Tinf
    Sp_spec = -4/D*h

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
            Sc[i] =Sc_spec
            Sp[i] =Sp_spec
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
    return T,node_x
end

rc("figure", figsize=(4.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.72)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])

PyPlot.close("all")
PyPlot.figure("5")
N_cv = [5]#linspace(5,100,7)
N_cv = round.(Int,N_cv)
T = []
for i = 1:length(N_cv)
    T,node_x = runsolver(N_cv[i])
    PyPlot.plot(node_x,T,label = "$(N_cv[i]) Volumes")
    PyPlot.pause(0.01)
end
PyPlot.xlabel("x")
PyPlot.ylabel("T (K)")
PyPlot.legend(loc = "best")

# PyPlot.figure("5c")
# N_cv = zeros(10)
# for i = 1:length(N_cv)
#     N_cv[i] = 2^i
# end
# N_cv2 = round.(Int,N_cv)
# T2 = zeros(length(N_cv))
# for i = 1:length(N_cv)
# # i = 1
#     T,node_x = runsolver(N_cv2[i])
#     T2[i] = T[1]
#     println(i)
# end
# PyPlot.semilogx(N_cv,T2)
# PyPlot.xlabel("Number Control Volumes")
# PyPlot.ylabel("T_max (K)")
