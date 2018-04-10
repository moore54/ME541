using PyPlot
using Distributions
close("all")

meshgrid(x,y) = (repmat(x',length(y),1),repmat(y,1,length(x)))


N_x = 4 #horizontal pressure nodes
N_y = 4 #vertical pressure nodes

Inlet_pos = 0.0
Outlet_pos = 0.05
Top_pos = 0.01
Bot_pos = 0.0

U0 = 0.001
V0 = 0.0001
P0 = 0.001
rho = 1000
mu = 0.001
u_relax = 0.5
v_relax = 0.5
p_relax = 1.0


# + 1 since we have ghost cells of 0 value
u_x = collect(linspace(Inlet_pos,Outlet_pos,N_x+1))
v_y = collect(linspace(Bot_pos,Top_pos,N_y+1))
u_x = [u_x[1]-(u_x[2]-u_x[1]);u_x;u_x[end]+(u_x[end]-u_x[end-1])]
v_y = [v_y[1]-(v_y[2]-v_y[1]);v_y;v_y[end]+(v_y[end]-v_y[end-1])]
P_x = (u_x[2:end]-u_x[1:end-1])./2+u_x[1:end-1]
P_y = (v_y[2:end]-v_y[1:end-1])./2+v_y[1:end-1]
u_y = P_y
v_x = [P_x;P_x[end]+(P_x[end]-P_x[end-1])]

X,Y=meshgrid(u_x,u_y)

# include ghost cells of 0 value
P_val = ones(length(P_y),length(P_x))*P0
# P_val[:,1] = 0.0 # outside domain
# P_val[:,end] = 0.0
# P_val[1,:] = 0.0
# P_val[end,:] = 0.0

u_val = ones(length(u_y),length(u_x))*U0
# u_val[:,1] = 0.0
# u_val[:,end] = 0.0
# u_val[1,:] = 0.0
# u_val[end,:] = 0.0

v_val = ones(length(v_y),length(v_x))*V0
# v_val[:,end] = 0.0
# v_val[:,1] = 0.0
# v_val[1,:] = 0.0
# v_val[end,:] = 0.0

plot(P_x,ones(P_x)*P_y[1],"x",label = "Pressure X Nodes")
plot(ones(P_y)*P_x[1],P_y,"x",label = "Pressure Y Nodes")
plot(u_x,ones(u_x)*u_y[1],".",label = "u X Nodes")
plot(ones(u_y)*u_x[1],u_y,".",label = "u Y Nodes")
plot(v_x,ones(v_x)*v_y[1],"*",label = "v X Nodes")
plot(ones(v_y)*v_x[1],v_y,"*",label = "v Y Nodes")
legend(loc = "best")

resid = 1E20
iter = 1
figure("u_contour")
# while resid>1E-3 || iter <1E3

#
#
#
#------------- 1) u_vel ------------#
#
#
#
# println("u_vel")

A_mat = zeros((length(u_x)-2)*(length(u_y)-2),(length(u_x)-2)*(length(u_y)-2))
b_arr = zeros(length(A_mat[:,1]))
Ap_u = zeros(b_arr)


Ae = 0.0
Aw = 0.0
An = 0.0
As = 0.0
Ap = 0.0
b_con = 0.0

j = 1+1 # x position
i = 1+1 # y position + 1 since we have ghost cells of 0 value
for k = 1:(length(u_val[:,1])-2)*(length(u_val[1,:])-2)
    # k = 2
    # j = 2+1

    if k==(length(u_val[1,:])-2)*(j-1)+1 || k==1
        #Apply Inlet Conditions
        Ae = 0.0
        Aw = 0.0
        An = 0.0
        As = 0.0
        Ap = 1.0
        b_con = U0

    elseif k==(length(u_val[1,:])-2)*(j-1)
        #Apply outlet conditions
        Ae = 0.0
        Aw = 1.0
        An = 0.0
        As = 0.0
        Ap = 1.0
        b_con = 0.0

    else

        De = mu/(u_x[i+1]-u_x[i])
        Fe = rho/2*(u_val[j,i+1]+u_val[j,i])
        Ae = (De + max(-Fe,0.0))*(v_y[j+1]-v_y[j])

        Dw = mu/(u_x[i]-u_x[i-1])
        Fw = rho/2*(u_val[j,i]+u_val[j,i-1])
        Aw = (Dw + max(Fw,0.0))*(v_y[j+1]-v_y[j])

        Dn = mu/(u_y[j+1]-u_y[j])
        Fn = rho/2*(v_val[j+1,i]+v_val[j+1,i-1])
        An = (Dn + max(-Fn,0.0))*(v_x[i+1]-v_x[i])

        Ds = mu/(u_y[j]-u_y[j-1])
        Fs = rho/2*(v_val[j,i]+v_val[j,i-1])
        As = (Ds + max(Fs,0.0))*(v_x[i+1]-v_x[i])

        Few = (Fe-Fw)*(v_y[j+1]-v_y[j])
        Fns = (Fn-Fs)*(v_x[i+1]-v_x[i])

        Ap = Ae+Aw+An+As+ Few + Fns


        #Apply bottom wall boundary conditions but not on inlet or outlet
        if j==2 && (k!=(length(u_val[1,:])-2)*(j-1)+1 && k!=1 && k!=(length(u_val[1,:])-2)*(j-1))
            Ap += mu/(u_y[j]-Bot_pos)*(v_x[i]-v_x[i-1]) - As
            As = 0.0
            # println("bottom wall")

        #Apply top wall boundary conditions but not on inlet or outlet
        elseif j == (length(u_val[:,1])-1) && (k!=(length(u_val[1,:])-2)*(j-1)+1 && k!=1 && k!=(length(u_val[1,:])-2)*(j-1))
            Ap += mu/(Top_pos-u_y[j])*(v_x[i]-v_x[i-1]) - An
            An = 0.0
            # println("top wall")
        end

        Ap = Ap/u_relax

        b_con = (P_val[j,i-1]-P_val[j,i])*(v_y[j+1]-v_y[j]) + 0.0 + (1-u_relax)*Ap*u_val[j,i]

    end


    #Assemble A matrix
    if k-(length(u_val[1,:])-2)>0
        A_mat[k,k-(length(u_val[1,:])-2)] = -As
    end

    if k-1>0
        A_mat[k,k-1] = -Aw
    end

    A_mat[k,k] = Ap

    if k+1<=(length(u_val[:,1])-2)*(length(u_val[1,:])-2)
        A_mat[k,k+1] = -Ae
    end

    if k+(length(u_val[1,:])-2)<=(length(u_val[:,1])-2)*(length(u_val[1,:])-2)
        A_mat[k,k+(length(u_val[1,:])-2)] = -An
    end

    #Record b and values that will be used in the pressure correction
    b_arr[k] = b_con

    Ap_u[k] = Ap

    # Increment indices
    i+=1
    if k==(length(u_val[1,:])-2)*(j-1)+1

        j+=1
        i=1+1
    end

    # println(A_mat[k,:])
end
# println("b_arr")
# println(b_arr)

# Solve for the u_velocity
u_val_new = A_mat\b_arr

# Reassemble u_val_new to correct dimensions
j = 1
for i = 1:length(u_val[:,1])-2
    u_val[i+1,2:end-1] = u_val_new[j:j+length(u_val[1,:])-3]
    j=j+length(u_val[1,:])-2
end


# Apply mass flowrate correction
Mout = sum(u_val[2:end-1,end-1])
Min = sum(u_val[2:end-1,2])
Min/Mout
u_val[2:end-1,end-1] = u_val[2:end-1,end-2]*Min/Mout

# println("v_vel")
#
#
#
#------------ 2) v_vel, use newest u_vel ----------#
#
#
#



A_mat = zeros((length(v_x)-1)*(length(v_y)-2),(length(v_x)-1)*(length(v_y)-2))
b_arr = zeros(length(A_mat[:,1]))
Ap_v = zeros(b_arr)

Ae = 0.0
Aw = 0.0
An = 0.0
As = 0.0
Ap = 0.0
b_con = 0.0

j = 1+1 # x position
i = 1+1 # y position + 1 since we have ghost cells of 0 value
for k = 1:(length(v_val[:,1])-2)*(length(v_val[1,:])-1)
    # k = 2
    # j = 2+1

    # Apply inlet condition
    if k==(length(v_val[1,:])-1)*(j-1)+1 || k==1

        Ae = 0.0
        Aw = 0.0
        An = 0.0
        As = 0.0
        Ap = 1.0
        b_con = 0.0

    elseif k==(length(v_val[1,:])-1)*(j-1)
        #Apply outlet conditions
        Ae = 0.0
        Aw = 1.0
        An = 0.0
        As = 0.0
        Ap = 1.0
        b_con = 0.0
        # println("outlet")

    else

        De = mu/(v_x[i+1]-v_x[i])
        Fe = rho/2*(u_val[j,i+1]+u_val[j-1,i+1])
        Ae = (De + max(-Fe,0.0))*(u_y[j]-u_y[j-1])

        Dw = mu/(v_x[i]-v_x[i-1])
        Fw = rho/2*(u_val[j,i]+u_val[j-1,i])
        Aw = (Dw + max(Fw,0.0))*(u_y[j]-u_y[j-1])

        Dn = mu/(v_y[j+1]-v_y[j])
        Fn = rho/2*(v_val[j,i]+v_val[j+1,i])
        An = (Dn + max(-Fn,0.0))*(u_x[i+1]-u_x[i])

        Ds = mu/(v_y[j]-v_y[j-1])
        Fs = rho/2*(v_val[j-1,i]+v_val[j,i])
        As = (Ds + max(Fs,0.0))*(u_x[i+1]-u_x[i])

        Few = (Fe-Fw)*(v_y[j+1]-v_y[j])
        Fns = (Fn-Fs)*(v_x[i+1]-v_x[i])

        Ap = Ae+Aw+An+As+ Few + Fns

        Ap = Ap/v_relax

        b_con = (P_val[j-1,i]-P_val[j,i])*(u_x[i+1]-u_x[i]) + 0.0 + (1-v_relax)*Ap*v_val[j,i]

        #Apply Top and Bottom Wall Conditions, but not on inlet or outlet
        if j==2 || j==(length(v_val[:,1])-1) && (k!=(length(v_val[1,:])-1)*(j-1)+1 && k!=1 && k!=(length(v_val[1,:])-1)*(j-1))
            Ae = 0.0
            Aw = 0.0
            An = 0.0
            As = 0.0
            Ap = 1.0
            b_con = 0.0
            # println("wall")
        end


    end


    if k-(length(v_val[1,:])-1)>0
        A_mat[k,k-(length(v_val[1,:])-1)] = -As
    end

    if k-1>0
        A_mat[k,k-1] = -Aw
    end

    A_mat[k,k] = Ap

    if k+1<=(length(v_val[:,1])-2)*(length(v_val[1,:])-2)
        A_mat[k,k+1] = -Ae
    end

    if k+(length(v_val[1,:])-1)<=(length(v_val[:,1])-2)*(length(v_val[1,:])-1)
        A_mat[k,k+(length(v_val[1,:])-1)] = -An
    end

    b_arr[k] = b_con

    Ap_v[k] = Ap

    # Increment indices
    i+=1
    if k==(length(v_val[1,:])-1)*(j-1)+1
        j+=1
        i=1+1
    end

    # println(A_mat[k,:])
end
# println("b_arr")
# println(b_arr)

# Solve for the v_velocity
v_val_new = A_mat\b_arr

# Reassemble u_val_new to correct dimensions
j = 1
for i = 1:length(v_val[:,1])-2
    v_val[i+1,2:end] = v_val_new[j:j+length(v_val[1,:])-2]
    j=j+length(v_val[1,:])-1
end

#
#
#
#--------- 3) Pressure --------#
# Set last column to 0 and don't solve for it, so if 4x4 pressure nodes, solving for 4x3
#
#
#
println("pressure")

A_mat = zeros((length(P_x)-3)*(length(P_y)-2),(length(P_x)-3)*(length(P_y)-2))
b_arr = zeros(length(A_mat[:,1]))
dw_save = zeros(length(A_mat[:,1]))
ds_save = zeros(length(A_mat[:,1]))

Ae = 0.0
Aw = 0.0
An = 0.0
As = 0.0
Ap = 0.0
b_con = 0.0

j = 1+1 # x position
i = 1+1 # y position + 1 since we have ghost cells of 0 value

for k = 1:(length(P_val[:,1])-2)*(length(P_val[1,:])-3)
    # k = 2
    # j = 2+1

    ae = v_y[j+1] - v_y[j]
    aw = v_y[j+1] - v_y[j]
    an = u_x[i+1] - u_x[i]
    as = u_x[i+1] - u_x[i]

# println(k)
# println(j)
# println(k+1+(j-2)*2)
    de = ae/Ap_u[k+1+(j-2)*2] # u_relax already applied, don't double apply
    dw = aw/Ap_u[k+(j-2)*2]
    dn = an/Ap_v[k+1+(length(v_val[1,:])-1)+(j-2)*3]
    ds = as/Ap_v[k+1+(j-2)*3]

    Ae = rho*de*ae
    Aw = rho*dw*aw
    An = rho*dn*an
    As = rho*ds*as



    Ap = Ae+Aw+An+As

    #Apply Top and Bottom Wall Conditions
    if j==2
        As = 0.0
        # println("bottom")
    elseif  j==(length(P_val[:,1])-1)
        An = 0.0
        # println("top")
    end


    if k==(length(P_val[1,:])-3)*(j-2)+1 || k==1
        # Apply inlet condition
        Aw = 0.0
        # println("inlet")
    elseif k==(length(P_val[1,:])-3)*(j-1)
        #Apply outlet conditions
        Ae = 0.0
        # println("outlet")
    end


    b_con = rho*u_val[j,i]*aw - rho*u_val[j,i+1]*ae + rho*v_val[j,i]*as - rho*v_val[j+1,i]*an


    # Assemble A matrix
    if k-(length(P_val[1,:])-3)>0
        A_mat[k,k-(length(P_val[1,:])-3)] = -As
    end

    if k-1>0
        A_mat[k,k-1] = -Aw
    end

    A_mat[k,k] = Ap

    if k+1<=(length(P_val[:,1])-2)*(length(P_val[1,:])-3)
        A_mat[k,k+1] = -Ae
    end

    if k+(length(P_val[1,:])-3)<=(length(P_val[:,1])-2)*(length(P_val[1,:])-3)
        A_mat[k,k+(length(P_val[1,:])-3)] = -An
    end

    b_arr[k] = b_con
    dw_save[k] = dw
    ds_save[k] = ds

    i+=1
    # Apply inlet index update
    if k==(length(P_val[1,:])-3)*(j-1) && k!=1
        # println("inlet")
        j+=1
        i=1+1
    end

    # println(round.(A_mat[k,:],3))
end
# println("b_arr")
# println(b_arr)
# Solve for the P correction
p_cor_new = A_mat\b_arr

# Reassemble p_cor_new to correct dimensions
p_cor = zeros(P_val)
dw_n = zeros(P_val)
ds_n = zeros(P_val)
j = 1
for i = 1:length(P_val[:,1])-2
    p_cor[i+1,2:end-2] = p_cor_new[j:j+length(p_cor[1,:])-4]
    dw_n[i+1,2:end-2] = dw_save[j:j+length(p_cor[1,:])-4]
    ds_n[i+1,2:end-2] = ds_save[j:j+length(p_cor[1,:])-4]
    j=j+length(p_cor[1,:])-3
end

P_val = P_val+p_relax*p_cor

for j = 1:length(u_val[:,1])-2
u_val[j,3:end-2] = u_val[j,3:end-2]+dw_n[j,2:end-2].*(p_cor[j,2:end-2]-p_cor[j,3:end-1])
end

for i = 1:length(v_val[1,:])-2
    v_val[3:end-2,i] = v_val[3:end-2,i]+ds_n[2:end-2,i].*(p_cor[2:end-2,i]-p_cor[3:end-1,i])
end
resid = norm(p_cor)
iter+=1

# if iter%100==0 || iter == 1
    # println(u_val)
    PyPlot.clf()
    CS = PyPlot.contour(X,Y,u_val,interpolation = "cubic",origin = "lower",cmap = PyPlot.cm[:gray])
    PyPlot.clabel(CS)
    # PyPlot.colorbar(orientation = "horizontal",extend = "both")
    PyPlot.pause(2.0)
# end
# end
