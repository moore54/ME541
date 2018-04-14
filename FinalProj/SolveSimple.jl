using PyPlot
using Distributions
close("all")

meshgrid(x,y) = (repmat(x',length(y),1),repmat(y,1,length(x)))


N_x = 4 #horizontal pressure nodes
N_y = 4 #vertical pressure nodes
second = false

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

Au = []
Av = []
AP = []

bu = []
bv = []
bP = []

# + 1 since we have ghost cells of 0 value
u_x = collect(linspace(Inlet_pos,Outlet_pos,N_x+1))
v_y = collect(linspace(Bot_pos,Top_pos,N_y+1))
u_x = [u_x[1]-(u_x[2]-u_x[1]);u_x;u_x[end]+(u_x[end]-u_x[end-1])]
v_y = [v_y[1]-(v_y[2]-v_y[1]);v_y;v_y[end]+(v_y[end]-v_y[end-1])]
P_x = (u_x[2:end]-u_x[1:end-1])./2+u_x[1:end-1]
P_y = (v_y[2:end]-v_y[1:end-1])./2+v_y[1:end-1]
u_y = P_y
v_x = [P_x;P_x[end]+(P_x[end]-P_x[end-1])]

Xu,Yu=meshgrid(u_x[2:end-1],u_y[2:end-1])
Xv,Yv=meshgrid(v_x[2:end-2],v_y[2:end-1])
Xp,Yp=meshgrid(P_x[2:end-1],P_y[2:end-1])

# include ghost cells of 0 value
P_val = ones(length(P_y),length(P_x))*P0
if second
    P_val_end=[0.00454,0.00303,0.00181,0.001,0.00465,0.0032,0.00196,0.001,0.00481,0.00339,0.00215,0.001,0.00499,0.00354,0.0023,0.001]
    j = 1
    for i = 1:length(P_val[:,1])-2
        P_val[i+1,2:end-1] = P_val_end[j:j+3]
        j=j+4
    end
end
# P_val[:,1] = 0.0 # outside domain
# P_val[:,end] = 0.0
# P_val[1,:] = 0.0
# P_val[end,:] = 0.0
P_val_star = copy(P_val)
P_val_new = copy(P_val)
#
# P_val[:,end-1] = 0.0

u_val = ones(length(u_y),length(u_x))*U0
if second
    uval_end = [0.001,0.00081,0.00077,0.00074,0.00085,0.001,0.00105,0.00102,0.00099,0.00113,0.001,0.00106,0.00103,0.00102,0.00114,0.001,0.00082,0.00079,0.00079,0.00088]
    j = 1
    for i = 1:length(u_val[:,1])-2
        u_val[i+1,2:end-1] = uval_end[j:j+length(u_val[1,:])-3]
        j=j+length(u_val[1,:])-2
    end
end
# u_val[:,1] = 0.0
# u_val[:,end] = 0.0
# u_val[1,:] = 0.0
# u_val[end,:] = 0.0
u_val_star = copy(u_val)
u_val_new = copy(u_val)

v_val = ones(length(v_y),length(v_x))*V0
if second
    v_val_end = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.43E-05,-3.98E-06,1.19E-06,7.15E-05,0.0,0.0,3.16E-06,2.82E-06,6.43E-06,8.66E-05,0.0,0.0,-9.83E-06,8.79E-06,8.50E-06,7.57E-05,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
    j = 1
    for i = 1:length(v_val[:,1])-2
        v_val[i+1,2:end] = v_val_end[j:j+length(v_val[1,:])-2]
        j=j+length(v_val[1,:])-1
    end
end
# v_val[:,end] = 0.0
# v_val[:,1] = 0.0
# v_val[1,:] = 0.0
# v_val[end,:] = 0.0
v_val_star = copy(v_val)
v_val_new = copy(v_val)
P_cor = zeros(P_val)

plot(P_x,ones(P_x)*P_y[1],"x",label = "Pressure X Nodes")
plot(ones(P_y)*P_x[1],P_y,"x",label = "Pressure Y Nodes")
plot(u_x,ones(u_x)*u_y[1],".",label = "u X Nodes")
plot(ones(u_y)*u_x[1],u_y,".",label = "u Y Nodes")
plot(v_x,ones(v_x)*v_y[1],"*",label = "v X Nodes")
plot(ones(v_y)*v_x[1],v_y,"*",label = "v Y Nodes")
legend(loc = "best")
Ae = 0.0
Aw = 0.0
An = 0.0
As = 0.0
Ap = 0.0
b_con = 0.0

Ae_matu = zeros(u_val)
Aw_matu = zeros(u_val)
An_matu = zeros(u_val)
As_matu = zeros(u_val)
Ap_matu = zeros(u_val)
b_matu = zeros(u_val)

A_matu = zeros((length(u_x)-2)*(length(u_y)-2),(length(u_x)-2)*(length(u_y)-2))
b_arru = zeros(length(A_matu[:,1]))

Ae_matv = zeros(v_val)
Aw_matv = zeros(v_val)
An_matv = zeros(v_val)
As_matv = zeros(v_val)
Ap_matv = zeros(v_val)
b_matv = zeros(v_val)

A_matv = zeros((length(v_x)-1)*(length(v_y)-2),(length(v_x)-1)*(length(v_y)-2))
b_arrv = zeros(length(A_matv[:,1]))

Ae_matP = zeros(P_val)
Aw_matP = zeros(P_val)
An_matP = zeros(P_val)
As_matP = zeros(P_val)
Ap_matP = zeros(P_val)
b_matP = zeros(P_val)

dw_n = zeros(P_val)
ds_n = zeros(P_val)

A_matP = zeros((length(P_x)-2)*(length(P_y)-3),(length(P_x)-2)*(length(P_y)-3))
b_arrP = zeros(length(A_matP[:,1]))


iter = 1
# while resid>1E-3 || iter <1E3
# for iter2 = 1:20

#
#
#
#------------- 1) u_vel ------------#
# GENERATE COEFFICIENTS
println("u_vel")

for j = 2:length(u_val[:,1])-1
    for i = 2:length(u_val[1,:])-1

        if i==2
            #Apply Inlet Conditions
            Ae = 0.0
            Aw = 0.0
            An = 0.0
            As = 0.0
            Ap = 1.0
            b_con = U0

        elseif i==(length(u_val[1,:])-1)
            #Apply outlet conditions
            Ae = 0.0
            Aw = 1.0
            An = 0.0
            As = 0.0
            Ap = 1.0
            b_con = 0.0

        else

            De = mu/(u_x[i+1]-u_x[i])
            Fe = rho/2*(u_val_star[j,i+1]+u_val_star[j,i])
            Ae = (De + max(-Fe,0.0))*(v_y[j+1]-v_y[j])

            Dw = mu/(u_x[i]-u_x[i-1])
            Fw = rho/2*(u_val_star[j,i]+u_val_star[j,i-1])
            Aw = (Dw + max(Fw,0.0))*(v_y[j+1]-v_y[j])

            Dn = mu/(u_y[j+1]-u_y[j])
            Fn = rho/2*(v_val_star[j+1,i]+v_val_star[j+1,i-1])
            An = (Dn + max(-Fn,0.0))*(v_x[i+1]-v_x[i])

            Ds = mu/(u_y[j]-u_y[j-1])
            Fs = rho/2*(v_val_star[j,i]+v_val_star[j,i-1])
            As = (Ds + max(Fs,0.0))*(v_x[i+1]-v_x[i])


            #Apply bottom wall boundary conditions but not on inlet or outlet
            if j==2

                Ds = mu/(u_y[j]-Bot_pos)
                # println(u_y[j]-Bot_pos)
                Fs = rho/2*(v_val_star[j,i]+v_val_star[j,i-1])
                As = (Ds + max(Fs,0.0))*(v_x[i+1]-v_x[i])

                Few = (Fe-Fw)*(v_y[j+1]-v_y[j])
                Fns = (Fn-Fs)*(v_x[i+1]-v_x[i])

                Ap = Ae+Aw+An+ Few + Fns
                Ap += mu/(u_y[j]-Bot_pos)*(v_x[i]-v_x[i-1])

                As = 0.0
                Ap = Ap/u_relax
                b_con = (P_val_star[j,i-1]-P_val_star[j,i])*(v_y[j+1]-Bot_pos) + 0.0 + (1-u_relax)*Ap*u_val_star[j,i]

                #Apply top wall boundary conditions but not on inlet or outlet
            elseif j == (length(u_val[:,1])-1)

                Dn = mu/(Top_pos-u_y[j])
                # println(Top_pos-u_y[j])
                Fn = rho/2*(v_val_star[j+1,i]+v_val_star[j+1,i-1])
                An = (Dn + max(-Fn,0.0))*(v_x[i+1]-v_x[i])

                Few = (Fe-Fw)*(v_y[j+1]-v_y[j])
                Fns = (Fn-Fs)*(v_x[i+1]-v_x[i])

                Ap = Ae+Aw+As+ Few + Fns
                Ap += mu/(Top_pos-u_y[j])*(v_x[i]-v_x[i-1])

                An = 0.0
                Ap = Ap/u_relax
                b_con = (P_val_star[j,i-1]-P_val_star[j,i])*(Top_pos-v_y[j]) + 0.0 + (1-u_relax)*Ap*u_val_star[j,i]
            else
                Few = (Fe-Fw)*(v_y[j+1]-v_y[j])
                Fns = (Fn-Fs)*(v_x[i+1]-v_x[i])

                Ap = Ae+Aw+An+As+ Few + Fns
                # println(u_y[j]-u_y[j-1])
                Ap = Ap/u_relax
                b_con = (P_val_star[j,i-1]-P_val_star[j,i])*(v_y[j+1]-v_y[j]) + 0.0 + (1-u_relax)*Ap*u_val_star[j,i]
            end

        end

        Ae_matu[j,i] = Ae
        Aw_matu[j,i] = Aw
        An_matu[j,i] = An
        As_matu[j,i] = As
        Ap_matu[j,i] = Ap
        b_matu[j,i] = b_con
    end
end

#ASSEMBLE U COEFFICIENTS

k = 1
for j = 2:length(u_val[:,1])-1
    for i = 2:length(u_val[1,:])-1

        #Assemble A matrix
        if k-(length(u_val[1,:])-2)>0
            A_matu[k,k-(length(u_val[1,:])-2)] = -As_matu[j,i]
        end

        if k-1>0
            A_matu[k,k-1] = -Aw_matu[j,i]
        end

        A_matu[k,k] = Ap_matu[j,i]

        if k+1<=(length(u_val[:,1])-2)*(length(u_val[1,:])-2)
            A_matu[k,k+1] = -Ae_matu[j,i]
        end

        if k+(length(u_val[1,:])-2)<=(length(u_val[:,1])-2)*(length(u_val[1,:])-2)
            A_matu[k,k+(length(u_val[1,:])-2)] = -An_matu[j,i]
        end
        b_arru[k] = b_matu[j,i]
        k+=1
    end

end

# println("A_matu")
# for i = 1:length(A_matu[:,1])
#     println(A_matu[i,:])
# end
# println(b_arru)

# Solve for the u_velocity
u_val_column = A_matu\b_arru

# Reassemble u_val_column to correct dimensions
j = 1
for i = 1:length(u_val[:,1])-2
    u_val_new[i+1,2:end-1] = u_val_column[j:j+length(u_val[1,:])-3]
    j=j+length(u_val[1,:])-2
end

# Apply mass flowrate correction
Mout = sum(u_val_new[2:end-1,end-1])
Min = sum(u_val_new[2:end-1,2])
Min/Mout
u_val_new[2:end-1,end-1] = u_val_new[2:end-1,end-2]*Min/Mout
# u_val_new = u_val_new[end:-1:1,:]
# println("uval")
# for i = 1:length(u_val[:,1])
#     println(u_val_new[i,:])
# end


#------------ 2) v_vel, use newest u_vel ----------#
#
#
#

for j = 2:length(v_val[:,1])-1
    for i = 2:length(v_val[1,:])-0
# println("$j $i ")
        if i==2
            #Apply inlet
            Ae = 0.0
            Aw = 0.0
            An = 0.0
            As = 0.0
            Ap = 1.0
            b_con = 0.0

        elseif i==length(v_val[1,:])-0
            #Apply outlet conditions
            Ae = 0.0
            Aw = 1.0
            An = 0.0
            As = 0.0
            Ap = 1.0
            b_con = 0.0
            # println("outlet")

        else
            #Apply Top and Bottom Wall Conditions, but not on inlet or outlet
            if j==2 || j==length(v_val[:,1])-1
                Ae = 0.0
                Aw = 0.0
                An = 0.0
                As = 0.0
                Ap = 1.0
                b_con = 0.0
                # println("wall")
            else
                De = mu/(v_x[i+1]-v_x[i])
                Fe = rho/2*(u_val_new[j,i+1]+u_val_new[j-1,i+1])
                Ae = (De + max(-Fe,0.0))*(u_y[j]-u_y[j-1])

                Dw = mu/(v_x[i]-v_x[i-1])
                Fw = rho/2*(u_val_new[j,i]+u_val_new[j-1,i])
                Aw = (Dw + max(Fw,0.0))*(u_y[j]-u_y[j-1])

                Dn = mu/(v_y[j+1]-v_y[j])
                Fn = rho/2*(v_val_star[j,i]+v_val_star[j+1,i])
                An = (Dn + max(-Fn,0.0))*(u_x[i+1]-u_x[i])

                Ds = mu/(v_y[j]-v_y[j-1])
                Fs = rho/2*(v_val_star[j-1,i]+v_val_star[j,i])
                As = (Ds + max(Fs,0.0))*(u_x[i+1]-u_x[i])

                Few = (Fe-Fw)*(v_y[j+1]-v_y[j])
                Fns = (Fn-Fs)*(v_x[i+1]-v_x[i])

                Ap = Ae+Aw+An+As+ Few + Fns
                Ap = Ap/v_relax
                b_con = (P_val_star[j-1,i]-P_val_star[j,i])*(u_x[i+1]-u_x[i]) + 0.0 + (1-v_relax)*Ap*v_val_star[j,i]
                # println("$j $i")
                # println(Fw)

            end
        end
        Ae_matv[j,i] = Ae
        Aw_matv[j,i] = Aw
        An_matv[j,i] = An
        As_matv[j,i] = As
        Ap_matv[j,i] = Ap
        b_matv[j,i] = b_con
    end
end

# println("b_matv")
# for i = 1:length(b_matv[:,1])
#     println(b_matv[i,:])
# end



k = 1
for j = 2:length(v_val[:,1])-1
    for i = 2:length(v_val[1,:])-0
        if k-(length(v_val[1,:])-0)>0
            A_matv[k,k-(length(v_val[1,:])-0)] = -As_matv[j,i]
        end

        if k-1>0
            A_matv[k,k-1] = -Aw_matv[j,i]
        end

        A_matv[k,k] = Ap_matv[j,i]

        if k+1<=(length(v_val[:,1])-2)*(length(v_val[1,:])-2)
            A_matv[k,k+1] = -Ae_matv[j,i]
        end

        if k+(length(v_val[1,:])-0)<=(length(v_val[:,1])-2)*(length(v_val[1,:])-2)
            A_matv[k,k+(length(v_val[1,:])-0)] = -An_matv[j,i]
        end

        b_arrv[k] = b_matv[j,i]
        k+=1

    end
end

# println("A_matv")
for i = 1:length(A_matv[:,1])
    println(A_matv[i,:])
end
# println("b_arrv")
# println(b_arrv)

v_val_column = A_matv\b_arrv

# Reassemble v_val_column to correct dimensions
j = 1
for i = 2:length(v_val[:,1])-1
    # println(j)
    v_val_new[i,2:end] = v_val_column[j:j+length(v_val_new[1,:])-2]
    j=j+length(v_val[1,:])-1
    # println(j-1)
end

# println("vval")
# for i = 1:length(v_val[:,1])
#     println(v_val_new[i,:])
# end

# #--------- 3) Pressure --------#
# # Set last column to 0 and don't solve for it, so if 4x4 pressure nodes, solving for 4x3
# #
# #
# #
#
#
# j = 2 # x position
# i = 2 # y position + 1 since we have ghost cells of 0 value
#
# for j = 2:(length(P_val[:,1])-1)
#     for i = 2:(length(P_val[1,:])-2)
#
#         ae = v_y[j+1] - v_y[j]
#         aw = v_y[j+1] - v_y[j]
#         an = u_x[i+1] - u_x[i]
#         as = u_x[i+1] - u_x[i]
#
#         de = ae/Ap_matu[j,i+1] # u_relax already applied, don't double apply
#         dw = aw/Ap_matu[j,i]
#         dn = an/Ap_matv[j,i]
#         ds = as/Ap_matv[j+1,i]
#
#         #Apply boundary Conditions for deltas
#         if j==2
#             #Bottom
#             ds = as/Ap_matv[j+1,i]/2.0
#             # as = as/2.0
#         elseif  j==(length(P_val[:,1])-1)
#             #Top
#             dn = an/Ap_matv[j,i]/2.0
#             # an = an/2.0
#         end
#
#         if i==1
#             # Apply inlet condition
#             dw = aw/Ap_matu[j,i]/2.0
#             # aw = aw/2.0
#         elseif i==(length(P_val[1,:])-3)
#             #Apply outlet conditions
#             # de = de/2.0
#             # ae = ae/2.0
#         end
#
#         Ae = rho*de*ae
#         Aw = rho*dw*aw
#         An = rho*dn*an
#         As = rho*ds*as
#         Ap = Ae+Aw+An+As
#
#         #Apply Top and Bottom Wall Conditions
#         if j==2
#             As = 0.0
#             # as =as/2
#             # println("bottom")
#         elseif  j==(length(P_val[:,1])-1)
#             An = 0.0
#             # an = an/2
#             # println("top")
#         end
#
#         if i==1
#             # Apply inlet condition
#             Aw = 0.0
#             # aw = aw/2
#             # println("inlet")
#         elseif i==(length(P_val[1,:])-3) # We aren't calculating outlet
#             #Apply outlet conditions
#             Ae = 0.0
#             # println("outlet")
#         end
#
#         b_con = rho*u_val_new[j,i]*aw - rho*u_val_new[j,i+1]*ae + rho*v_val_new[j+1,i]*as - rho*v_val_new[j,i]*an
#
#         Ae_matP[j,i] = Ae
#         Aw_matP[j,i] = Aw
#         An_matP[j,i] = An
#         As_matP[j,i] = As
#         Ap_matP[j,i] = Ap
#         b_matP[j,i] = b_con
#
#         dw_n[j,i] = dw
#         ds_n[j,i] = ds
#
#     end
# end
#
# # println("b_matP")
# # for i = 1:length(b_matP[:,1])
# #     println(b_matP[i,:])
# # end
#
# #Assemble coefficients
# k=1
# for j = 2:(length(P_val[:,1])-1)
#     for i = 2:(length(P_val[1,:])-2)
#
#         # Assemble A matrix
#         if k-(length(P_val[1,:])-3)>0
#             A_matP[k,k-(length(P_val[1,:])-3)] = -As_matP[j,i]
#         end
#
#         if k-1>0
#             A_matP[k,k-1] = -Aw_matP[j,i]
#         end
#
#         A_matP[k,k] = Ap_matP[j,i]
#
#         if k+1<=(length(P_val[:,1])-2)*(length(P_val[1,:])-3)
#             A_matP[k,k+1] = -Ae_matP[j,i]
#         end
#
#         if k+(length(P_val[1,:])-3)<=(length(P_val[:,1])-2)*(length(P_val[1,:])-3)
#             A_matP[k,k+(length(P_val[1,:])-3)] = -An_matP[j,i]
#         end
#
#         b_arrP[k] = b_matP[j,i]
#
#         k+=1
#     end
# end
#
# # println("b_arr")
# # println(b_arr)
#
# # Solve for the P correction
# P_cor_new = A_matP\b_arrP
#
# # Reassemble P_cor_new to correct dimensions
#
#
# j = 1
# for i = 1:length(P_val[:,1])-2
#     P_cor[i+1,2:end-2] = P_cor_new[j:j+length(P_cor[1,:])-4]
#     j=j+length(P_cor[1,:])-3
# end
#
#
# P_val_new = P_val_new+p_relax*P_cor
#
# for j = 2:length(u_val[:,1])-1
#     u_val_new[j,3:end-2] = u_val_new[j,3:end-2]+dw_n[j,3:end-1].*(P_cor[j,2:end-2]-P_cor[j,3:end-1])
# end
#
# for j = 3:length(v_val[:,1])-2
#     v_val_new[j,3:end] = v_val_new[j,3:end]+ds_n[j,2:end].*(P_cor[j-1,2:end]-P_cor[j,2:end])
# end
#
# # v_val_new[:,1] = 0.0
# # v_val[:,end] = 0.0
# # P_val[:,end-1] = 0.0
#
# # println("u_val")
# # for i = 1:length(u_val[:,1])
# #     println(u_val_new[i,:])
# # end
# #
# # println("v_val")
# # for i = 1:length(v_val[:,1])
# #     println(v_val_new[i,:])
# # end
# #
# # println("P_val")
# # for i = 1:length(P_val[:,1])
# #     println(P_val_new[i,:])
# # end
#
# u_val = copy(u_val_new)
# v_val = copy(v_val_new)
# P_val = copy(P_val_new)
# # if iter2 ==2
# #     break
# # end
# # println("!!! $iter2")
# u_val_star = copy(u_val_new)
# v_val_star = copy(v_val_new)
# P_val_star = copy(P_val_new)
# #
# #
# # resid = norm(P_cor)
# # iter+=1 
# #
# # # if iter%100==0 || iter == 1
# #     # println(u_val)
# #     # figure("u_contour")
# #     # PyPlot.clf()
# #     # CS = PyPlot.contourf(Xu,Yu,u_val[2:end-1,2:end-1],interpolation = "cubic",origin = "lower",cmap = PyPlot.cm[:viridis])
# #     # PyPlot.clabel(CS)
# #     # PyPlot.colorbar(orientation = "horizontal",extend = "both")
# #     # # PyPlot.pause(0.00005)
# #     #
# #     # figure("v_contour")
# #     # PyPlot.clf()
# #     # CS = PyPlot.contourf(Xv,Yv,v_val[2:end-1,2:end-2],interpolation = "cubic",origin = "lower",cmap = PyPlot.cm[:viridis])
# #     # PyPlot.clabel(CS)
# #     # PyPlot.colorbar(orientation = "horizontal",extend = "both")
# #     # # PyPlot.pause(0.00005)
# #     #
# #     # figure("p_contour")
# #     # PyPlot.clf()
# #     # CS = PyPlot.contourf(Xp,Yp,P_val[2:end-1,2:end-1],interpolation = "cubic",origin = "lower",cmap = PyPlot.cm[:viridis])
# #     # PyPlot.clabel(CS)
# #     # PyPlot.colorbar(orientation = "horizontal",extend = "both")
# #     # PyPlot.pause(0.0000000005)
# #
# # # y_anal = linspace(Bot_pos,Top_pos,100)
# # # G = (mean(P_val[:,2]-P_val[:,end-1]))/(u_x[end-1]-u_x[2])
# # # u_anal = G/(2*mu)*y_anal.*(Top_pos-y_anal)
# # #
# # # figure("outlet")
# # # PyPlot.plot(u_val[:,end-1],u_y)
# # # PyPlot.plot(u_anal,y_anal)
# # # PyPlot.pause(0.5)
# # # end
# # # for i = 1:length(P_val[:,1])
# # #     println(P_val[i,:])
# # # end
# # # end
