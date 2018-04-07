using PyPlot
close("all")

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

#TODO: add wall boundary conditions to x momentum, remove them from y, set y momentum walls to v = 0, do p boundaries

# + 1 since we have ghost cells of 0 value
u_x = collect(linspace(Inlet_pos,Outlet_pos,N_x+1))
v_y = collect(linspace(Bot_pos,Top_pos,N_y+1))
u_x = [u_x[1]-(u_x[2]-u_x[1]);u_x;u_x[end]+(u_x[end]-u_x[end-1])]
v_y = [v_y[1]-(v_y[2]-v_y[1]);v_y;v_y[end]+(v_y[end]-v_y[end-1])]
P_x = (u_x[2:end]-u_x[1:end-1])./2+u_x[1:end-1]
P_y = (v_y[2:end]-v_y[1:end-1])./2+v_y[1:end-1]
u_y = P_y
v_x = [P_x;P_x[end]+(P_x[end]-P_x[end-1])]

# include ghost cells of 0 value
P_val = ones(length(P_y),length(P_x))*P0
P_val[:,1] = 0.0 # outside domain
P_val[:,end] = 0.0
P_val[1,:] = 0.0
P_val[end,:] = 0.0

u_val = ones(length(u_y),length(u_x))*U0
u_val[:,1] = 0.0
u_val[:,end] = 0.0
u_val[1,:] = 0.0
u_val[end,:] = 0.0

v_val = ones(length(v_y),length(v_x))*V0
v_val[:,end] = 0.0
v_val[:,1] = 0.0
v_val[1,:] = 0.0
v_val[end,:] = 0.0

plot(P_x,ones(P_x)*P_y[1],"x",label = "Pressure X Nodes")
plot(ones(P_y)*P_x[1],P_y,"x",label = "Pressure Y Nodes")
plot(u_x,ones(u_x)*u_y[1],".",label = "u X Nodes")
plot(ones(u_y)*u_x[1],u_y,".",label = "u Y Nodes")
plot(v_x,ones(v_x)*v_y[1],"*",label = "v X Nodes")
plot(ones(v_y)*v_x[1],v_y,"*",label = "v Y Nodes")
legend(loc = "best")

#
#
#
#------------- 1) u_vel ------------#
#
#
#

A_mat = zeros((length(u_x)-2)*(length(u_y)-2),(length(u_x)-2)*(length(u_y)-2))
b_arr = zeros(length(A_mat[:,1]))
Ai1J = zeros(b_arr)
AiJ = zeros(b_arr)

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

    if k==(length(u_val[1,:])-2)*(i-1)+1 || k==1
        #Apply Inlet Conditions
        Ae = 0.0
        Aw = 0.0
        An = 0.0
        As = 0.0
        Ap = 1.0
        b_con = V0

        if k!=1
            i+=1
            j=1+1
        end

    elseif k==(length(u_val[1,:])-2)*(i-1)
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

        b_con = (P_val[j,i-1]-P_val[j,i])*(v_y[j+1]-v_y[j]) + 0.0 + (1-u_relax)*Ap/u_relax*u_val[j,i]

    end

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

    b_arr[k] = b_con

    Ai1J[k] = Ae
    AiJ[k] = Ap

    j+=1
    # println(A_mat[k,:])
end

# Solve for the u_velocity
u_val_new = A_mat\b_arr

# Reassemble u_val_new to correct dimensions
j = 1
for i = 1:length(u_val[:,1])-2
    u_val[i+1,2:end-1] = u_val_new[j:j+length(u_val[1,:])-3]
    j=j+length(u_val[1,:])-2
end

Mout = sum(u_val[:,end-1])
Min = sum(u_val[:,2])
u_val[:,end-1] = u_val[:,end-2]*Min/Mout


#
#
#
#------------ 2) v_vel, use newest u_vel ----------#
#
#
#



A_mat = zeros((length(v_x)-2)*(length(v_y)-2),(length(v_x)-2)*(length(v_y)-2))
b_arr = zeros(length(A_mat[:,1]))
AIj1 = zeros(b_arr)
AIj = zeros(b_arr)

Ae = 0.0
Aw = 0.0
An = 0.0
As = 0.0
Ap = 0.0
b_con = 0.0

j = 1+1 # x position
i = 1+1 # y position + 1 since we have ghost cells of 0 value
for k = 1:(length(v_val[:,1])-2)*(length(v_val[1,:])-2)
    # k = 2
    # j = 2+1

    if k==(length(v_val[1,:])-2)*(i-1)+1
        i+=1
        j=1+1
    end

    if j==2
        #Apply Bottom Wall Conditions
        De = mu/(v_x[i]-v_x[i-1])
        Fe = rho/2*(u_val[j,i+1]+u_val[j-1,i+1])
        Ae = (De + max(-Fe,0.0))*(u_y[j]-u_y[j-1])

        Dw = mu/(v_x[i+1]-v_x[i])
        Fw = rho/2*(u_val[j,i]+u_val[j-1,i])
        Aw = (Dw + max(Fw,0.0))*(u_y[j]-u_y[j-1])

        Dn = mu/(v_y[j+1]-v_y[j])
        Fn = rho/2*(v_val[j,i]+v_val[j+1,i])
        An = (Dn + max(-Fn,0.0))*(u_x[i+1]-u_x[i])

        Ds = mu/(v_y[j]-v_y[j-1])
        Fs = rho/2*(v_val[j-1,i]+v_val[j,i])
        As = 0.0#(Ds + max(Fs,0.0))*(u_x[i+1]-u_x[i])

        Few = (Fe-Fw)*(v_y[j+1]-v_y[j])
        Fns = (Fn-Fs)*(v_x[i+1]-v_x[i])

        Ap = Ae+Aw+An+As+ Few + Fns

        b_con = mu/(u_y[j]-Bot_pos)*(v_x[i+1]-v_x[i]) + (P_val[j-1,i]-P_val[j,i])*(u_x[i+1]-u_x[i]) + 0.0 + (1-v_relax)*Ap/v_relax*v_val[j,i]

    elseif j==(length(v_val[:,1])-2)
        #Apply Top wall conditions
        De = mu/(v_x[i]-v_x[i-1])
        Fe = rho/2*(u_val[j,i+1]+u_val[j-1,i+1])
        Ae = (De + max(-Fe,0.0))*(u_y[j]-u_y[j-1])

        Dw = mu/(v_x[i+1]-v_x[i])
        Fw = rho/2*(u_val[j,i]+u_val[j-1,i])
        Aw = (Dw + max(Fw,0.0))*(u_y[j]-u_y[j-1])

        Dn = mu/(v_y[j+1]-v_y[j])
        Fn = rho/2*(v_val[j,i]+v_val[j+1,i])
        An = 0.0#(Dn + max(-Fn,0.0))*(u_x[i+1]-u_x[i])

        Ds = mu/(v_y[j]-v_y[j-1])
        Fs = rho/2*(v_val[j-1,i]+v_val[j,i])
        As = (Ds + max(Fs,0.0))*(u_x[i+1]-u_x[i])

        Few = (Fe-Fw)*(v_y[j+1]-v_y[j])
        Fns = (Fn-Fs)*(v_x[i+1]-v_x[i])

        Ap = Ae+Aw+An+As+ Few + Fns

        b_con = mu/(u_y[j]-Bot_pos)*(v_x[i+1]-v_x[i]) + (P_val[j-1,i]-P_val[j,i])*(u_x[i+1]-u_x[i]) + 0.0 + (1-v_relax)*Ap/v_relax*v_val[j,i]


    else

        De = mu/(v_x[i]-v_x[i-1])
        Fe = rho/2*(u_val[j,i+1]+u_val[j-1,i+1])
        Ae = (De + max(-Fe,0.0))*(u_y[j]-u_y[j-1])

        Dw = mu/(v_x[i+1]-v_x[i])
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

        b_con = (P_val[j-1,i]-P_val[j,i])*(u_x[i+1]-u_x[i]) + 0.0 + (1-v_relax)*Ap/v_relax*v_val[j,i]

    end

    if k-(length(v_val[1,:])-2)>0
        A_mat[k,k-(length(v_val[1,:])-2)] = -As
    end

    if k-1>0
        A_mat[k,k-1] = -Aw
    end

    A_mat[k,k] = Ap

    if k+1<=(length(v_val[:,1])-2)*(length(v_val[1,:])-2)
        A_mat[k,k+1] = -Ae
    end

    if k+(length(v_val[1,:])-2)<=(length(v_val[:,1])-2)*(length(v_val[1,:])-2)
        A_mat[k,k+(length(v_val[1,:])-2)] = -An
    end

    b_arr[k] = b_con

    AIj1[k] = An
    AIj[k] = Ap

    j+=1
    # println(A_mat[k,:])
end

# Solve for the u_velocity
v_val_new = A_mat\b_arr

# Reassemble u_val_new to correct dimensions
j = 1
for i = 1:length(u_val[:,1])-2
    v_val[i+1,2:end-1] = v_val_new[j:j+length(v_val[1,:])-3]
    j=j+length(v_val[1,:])-2
end


# v_val_new = v_val
#
# # 3) Pressure
#
# A1 = v_y[j+1] - v_y[j]
# A2 = v_y[j+1] - v_y[j]
# A3 = u_x[i+1] - u_x[i]
# A4 = u_x[i+1] - u_x[i]
#
# d1 = A1*u_relax/Ai1J
# d2 = A2*u_relax/AiJ
# d3 = A3*u_relax/AIj1
# d4 = A4*u_relax/AIj
#
#
# Ae = rho*d1*A1
# Aw = rho*d2*A2
# An = rho*d3*A3
# As = rho*d4*A4
# Ap = Ae+Aw+An+As
# b_con = rho*u_val_new[j,i]*A2 - rho*u_val_new[j,i+1]*A1 - rho*v_val_new[j,i]*A4 - rho*v_val_new[j+1,i]*A3
