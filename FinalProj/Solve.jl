using PyPlot
using Dierckx
using Distributions
close("all")
rc("figure", figsize=(5.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.9)
rc("image",interpolation = "spline16")

meshgrid(x,y) = (repmat(x',length(y),1),repmat(y,1,length(x)))

NG = 4
factor = linspace(1,NG,2)
end_error = zeros(factor)
final_iter = zeros(factor)
factor = round.(Int, factor)
total_error = 0.0
N_nodes = (4*2.^factor).^2

u_relax = 0.6
v_relax = 0.6
p_relax = 1.0


u_val = []
u_y = []
u_anal = []
y_anal = []
iter = []

Xp = []
Yp = []
u_val = []
v_val = []
P_val= []


Outlet_pos = []
Inlet_pos = []
Top_pos = []
Bot_pos = []


for grid_i = 1:NG

    refine = 2^factor[grid_i]

    N_x = 4*refine #horizontal pressure nodes
    N_y = 4*refine #vertical pressure nodes
    N_xt = N_x+2 #horizontal pressure nodes
    N_yt = N_y+2 #vertical pressure nodes
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

    delx = v_x[2]-v_x[1]
    dely = v_y[2]-v_y[1]

    Xu,Yu=meshgrid(u_x[2:end],u_y)
    Xv,Yv=meshgrid(v_x[2:end],v_y[1:end-1])
    Xp,Yp=meshgrid(P_x,P_y)

    # include ghost cells of 0 value
    P_val = ones(N_xt,N_yt)*P0
    P_val[:,1] = 0.0 # outside domain
    P_val[:,end] = 0.0
    P_val[1,:] = 0.0
    P_val[end,:] = 0.0
    P_val_star = copy(P_val)

    #
    # P_val[:,end-1] = 0.0

    u_val = ones(N_xt,N_yt)*U0
    u_val[:,1] = 0.0
    # u_val[:,end] = 0.0
    u_val[1,:] = 0.0
    u_val[end,:] = 0.0
    u_val_star = copy(u_val)
    u_val_new = copy(u_val)

    v_val = ones(N_xt,N_yt)*V0
    v_val[:,end] = 0.0
    v_val[:,1] = 0.0
    # v_val[1,:] = 0.0
    v_val[end,:] = 0.0
    v_val_star = copy(v_val)
    v_val_new = copy(v_val)
    P_cor = zeros(P_val)

    # plot(P_x,ones(P_x)*P_y[1],"x",label = "Pressure X Nodes")
    # plot(ones(P_y)*P_x[1],P_y,"x",label = "Pressure Y Nodes")
    # plot(u_x,ones(u_x)*u_y[1],".",label = "u X Nodes")
    # plot(ones(u_y)*u_x[1],u_y,".",label = "u Y Nodes")
    # plot(v_x,ones(v_x)*v_y[1],"*",label = "v X Nodes")
    # plot(ones(v_y)*v_x[1],v_y,"*",label = "v Y Nodes")
    # legend(loc = "best")
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

    A_matu = eye(N_xt*N_yt)
    b_arru = zeros(length(A_matu[:,1]))

    Ae_matv = zeros(v_val)
    Aw_matv = zeros(v_val)
    An_matv = zeros(v_val)
    As_matv = zeros(v_val)
    Ap_matv = zeros(v_val)
    b_matv = zeros(v_val)

    A_matv = eye(N_xt*N_yt)
    b_arrv = zeros(length(A_matv[:,1]))

    Ae_matP = zeros(P_val)
    Aw_matP = zeros(P_val)
    An_matP = zeros(P_val)
    As_matP = zeros(P_val)
    Ap_matP = zeros(P_val)
    b_matP = zeros(P_val)

    dw_n = zeros(P_val)
    ds_n = zeros(P_val)

    A_matP = eye(N_xt*N_yt)
    b_arrP = zeros(length(A_matP[:,1]))

    resid = 1.0
    iter = 1
    while resid>1E-6 && iter<1E3
        # for iter2 = 1:20
        # N_xt = N_xt-1
        # N_yt = N_yt-1

        #
        #
        #
        #------------- 1) u_vel ------------#
        # GENERATE COEFFICIENTS
        # println("u_vel")

        for i = 2:N_xt
            for j = 2:N_yt-1
                # println("$j $i")
                if i==2
                    #Apply Inlet Conditions
                    Ae = 0.0
                    Aw = 0.0
                    An = 0.0
                    As = 0.0
                    Ap = 1.0
                    b_con = U0

                    A_matu[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)] = Ap
                    b_arru[1+(j-1)*N_xt+(i-1)] = b_con

                elseif i==(N_xt)
                    #Apply outlet conditions
                    Ae = 0.0
                    Aw = 1.0
                    An = 0.0
                    As = 0.0
                    Ap = 1.0
                    b_con = 0.0

                    A_matu[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)-1] = -Aw
                    A_matu[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)] = Ap
                    b_arru[1+(j-1)*N_xt+(i-1)] = b_con

                else

                    De = mu/delx
                    Fe = rho/2*(u_val_star[j,i+1]+u_val_star[j,i])
                    Ae = (De + max(-Fe,0.0))*dely

                    Dw = mu/delx
                    Fw = rho/2*(u_val_star[j,i]+u_val_star[j,i-1])
                    Aw = (Dw + max(Fw,0.0))*dely

                    Dn = mu/dely
                    Fn = rho/2*(v_val_star[j-1,i]+v_val_star[j-1,i-1])
                    An = (Dn + max(-Fn,0.0))*delx

                    Ds = mu/dely
                    Fs = rho/2*(v_val_star[j,i]+v_val_star[j,i-1])
                    As = (Ds + max(Fs,0.0))*delx

                    Few = (Fe-Fw)*dely
                    Fns = (Fn-Fs)*delx

                    Ap = Ae+Aw+An+As+ Few + Fns
                    Ap = Ap/u_relax

                    #Apply top wall boundary conditions but not on inlet or outlet
                    if j==2

                        An = 0.0
                        Ap = Ae+Aw+As+An+ Few + Fns
                        Ap += mu/(dely/2)*delx

                        Ap = Ap/u_relax

                        #Apply top wall boundary conditions but not on inlet or outlet
                    elseif j == N_yt-1

                        As = 0.0
                        Ap = Ae+Aw+As+An+ Few + Fns
                        Ap += mu/(dely/2)*delx

                        Ap = Ap/u_relax

                    end

                    b_con = (P_val_star[j,i-1]-P_val_star[j,i])*(dely) + 0.0 + (1-u_relax)*Ap*u_val_star[j,i]

                    A_matu[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)+1] = -Ae
                    A_matu[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)-1] = -Aw
                    A_matu[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)-N_xt] = -An
                    A_matu[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)+N_xt] = -As
                    A_matu[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)] = Ap
                    b_arru[1+(j-1)*N_xt+(i-1)] = b_con

                end

                Ae_matu[j,i] = Ae
                Aw_matu[j,i] = Aw
                An_matu[j,i] = An
                As_matu[j,i] = As
                Ap_matu[j,i] = Ap
                b_matu[j,i] = b_con

                # println("$j $i")

            end
        end

        # println("A_matu")
        # for i = 1:length(A_matu[:,1])
        #     println(A_matu[i,:])
        # end
        # println(b_arru)

        # Solve for the u_velocity
        u_val_column = A_matu\b_arru
        # println(u_val_column)
        # Reassemble u_val_column to correct dimensions
        i = 1
        for j = 1:N_yt
            u_val_new[j,:] = u_val_column[i:i+N_xt-1]
            i = i+N_xt
        end

        # Apply mass flowrate correction
        Mout = sum(u_val_new[2:end-1,end-1])
        Min = sum(u_val_new[2:end-1,2])

        u_val_new[2:end-1,end] = u_val_new[2:end-1,end-1]*Min/Mout

        # println("uval")
        # for i = 1:length(u_val[:,1])
        #     println(u_val_new[i,:])
        # end


        #------------ 2) v_vel, use newest u_vel ----------#
        #
        #
        #

        for i = 1:N_xt
            for j = 1:N_yt-1

                if i==1
                    #Apply inlet
                    Ae = 0.0
                    Aw = 0.0
                    An = 0.0
                    As = 0.0
                    Ap = 1.0
                    b_con = 0.0

                    A_matv[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)] = Ap
                    b_arrv[1+(j-1)*N_xt+(i-1)] = b_con

                elseif i==N_xt
                    #Apply outlet conditions
                    Ae = 0.0
                    Aw = 1.0
                    An = 0.0
                    As = 0.0
                    Ap = 1.0
                    b_con = 0.0
                    # println("outlet")

                    A_matv[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)-1] = -Aw
                    A_matv[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)] = Ap
                    b_arrv[1+(j-1)*N_xt+(i-1)] = b_con

                else
                    # Apply Top and Bottom Wall Conditions, but not on inlet or outlet
                    if j==1 || j==N_yt-1
                        Ae = 0.0
                        Aw = 0.0
                        An = 0.0
                        As = 0.0
                        Ap = 1.0
                        b_con = 0.0

                        A_matv[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)] = Ap
                        b_arrv[1+(j-1)*N_xt+(i-1)] = b_con
                    else
                        De = mu/delx
                        Fe = rho/2*(u_val_new[j,i+1]+u_val_new[j+1,i+1])
                        Ae = (De + max(-Fe,0.0))*dely

                        Dw = mu/delx
                        Fw = rho/2*(u_val_new[j,i]+u_val_new[j+1,i])
                        Aw = (Dw + max(Fw,0.0))*dely

                        Dn = mu/dely
                        Fn = rho/2*(v_val_star[j-1,i]+v_val_star[j,i])
                        An = (Dn + max(-Fn,0.0))*delx

                        Ds = mu/dely
                        Fs = rho/2*(v_val_star[j+1,i]+v_val_star[j,i])
                        As = (Ds + max(Fs,0.0))*delx

                        Few = (Fe-Fw)*dely
                        Fns = (Fn-Fs)*delx

                        Ap = Ae+Aw+An+As+ Few + Fns
                        Ap = Ap/v_relax
                        b_con = (P_val_star[j+1,i]-P_val_star[j,i])*delx + 0.0 + (1-v_relax)*Ap*v_val_star[j,i]

                        A_matv[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)+1] = -Ae
                        A_matv[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)-1] = -Aw
                        A_matv[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)-N_xt] = -An
                        A_matv[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)+N_xt] = -As
                        A_matv[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)] = Ap
                        b_arrv[1+(j-1)*N_xt+(i-1)] = b_con

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


        # println("A_matv")
        # for i = 1:length(A_matv[:,1])
        #     println(A_matv[i,:])
        # end
        # println(b_arrv)

        v_val_column = A_matv\b_arrv

        # Reassemble v_val_column to correct dimensions
        i = 1
        for j = 1:N_yt
            v_val_new[j,:] = v_val_column[i:i+N_xt-1]
            i = i+N_xt
        end

        # println("uval")
        # for i = 1:length(u_val[:,1])
        #     println(u_val_new[i,:])
        # end
        #
        # println("vval")
        # for i = 1:length(v_val_new[:,1])
        #     println(v_val_new[i,:])
        # end

        #
        #
        #
        #
        #--------- 3) Pressure --------#
        # Set last column to 0 and don't solve for it, so if 4x4 pressure nodes, solving for 4x3
        #
        #
        #

        for i = 2:N_xt-1
            for j = 2:N_yt-1

                ae = dely
                aw = dely
                an = delx
                as = delx

                de = ae/Ap_matu[j,i+1] # u_relax already applied, don't double apply
                dw = aw/Ap_matu[j,i]
                dn = an/Ap_matv[j-1,i]
                ds = as/Ap_matv[j,i]

                # println("$de $dw $ds $dn")

                #Apply boundary Conditions for deltas
                if j==2 #top
                    dn = dn/2.0
                elseif  j==N_yt-1 #bottom
                    ds = ds/2.0
                end

                if i==2 # Apply inlet condition
                    dw = dw/2.0
                end #outlet not calculated or included due to the x-1 nature

                Ae = rho*de*ae
                Aw = rho*dw*aw
                An = rho*dn*an
                As = rho*ds*as
                Ap = Ae+Aw+An+As

                b_con = rho*u_val_new[j,i]*aw - rho*u_val_new[j,i+1]*ae + rho*v_val_new[j,i]*as - rho*v_val_new[j-1,i]*an
                # println("$j $i $(v_val_new[j-1,i]) ")
                #Apply Top and Bottom Wall Conditions
                if j==2 #top
                    An = 0.0
                elseif  j==N_yt-1 #bot
                    As = 0.0
                end

                if i==2 # Apply inlet condition
                    Aw = 0.0
                elseif i==N_xt-2 #Apply outlet conditions
                    Ae = 0.0
                end

                if i<N_xt-1
                    A_matP[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)+1] = -Ae
                    A_matP[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)-1] = -Aw
                    A_matP[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)-N_xt] = -An
                    A_matP[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)+N_xt] = -As
                    A_matP[1+(j-1)*N_xt+(i-1),1+(j-1)*N_xt+(i-1)] = Ap
                    b_arrP[1+(j-1)*N_xt+(i-1)] = b_con
                end


                Ae_matP[j,i] = Ae
                Aw_matP[j,i] = Aw
                An_matP[j,i] = An
                As_matP[j,i] = As
                Ap_matP[j,i] = Ap
                b_matP[j,i] = b_con

                dw_n[j,i] = dw
                ds_n[j,i] = ds

            end
        end


        # println("A_matP")
        # for i = 1:length(A_matP[:,1])
        #     println(A_matP[i,:])
        # end
        # println(b_arrP)

        # Solve for the P correction
        P_cor_new_col = A_matP\b_arrP
        # println("P_cor")
        # println(P_cor_new_col)
        # Reassemble P_cor_new to correct dimensions


        i = 1
        for j = 1:N_yt
            P_cor[j,:] = P_cor_new_col[i:i+N_xt-1]
            i = i+N_xt
        end

        # println("P_cor")
        # for i = 1:length(P_cor[:,1])
        #     println(P_cor[i,:])
        # end


        P_val = P_val+p_relax*P_cor
        # P_val[:,end] = 0.0
        # println("P_val")
        # for i = 1:length(P_val[:,1])
        #     println(P_val[i,:])
        # end


        for j = 1:N_yt
            u_val[j,3:end] = u_val_new[j,3:end]+dw_n[j,3:end].*(P_cor[j,2:end-1]-P_cor[j,3:end])
        end

        for j = 1:N_yt-2
            v_val[j,2:end-1] = v_val_new[j,2:end-1]+ds_n[j,2:end-1].*(P_cor[j+1,2:end-1]-P_cor[j,2:end-1])
        end
        v_val[end-1,:] = v_val_new[end-1,:]
        # println("uval")
        # for i = 1:length(u_val[:,1])
        #     println(u_val[i,:])
        # end

        # println("vval")
        # for i = 1:length(v_val[:,1])
        #     println(v_val[i,:])
        # end



        u_val_star = (u_val_new)
        v_val_star = (v_val_new)
        P_val_star = (P_val)


        resid = norm(P_cor)

        println(resid)
        if iter%10==0 || iter == 1

            y_anal = u_y
            # ell = linspace(Inlet_pos,Outlet_pos,4)
            # for i = 1:length(ell)
            ell = Outlet_pos-Inlet_pos
            D_H = (Top_pos-Bot_pos)*2
            ReD = rho/mu*U0*D_H
            eff = 24/ReD
            DP = eff*4*ell/D_H*(0.5*rho*U0^2)
            DPDX = DP/ell
            u_anal = DPDX/(2*mu)*y_anal.*(Top_pos-y_anal)
            for i = 1:length(u_anal)
                if u_anal[i]<0
                    u_anal[i] = 0.0
                end
            end

            # G = (mean(P_val[:,2]-P_val[:,end]))/(u_x[end-1]-u_x[2])
            figname = "outlet_$grid_i"
            figure("outlet_$grid_i")
            PyPlot.plot(u_val[:,end],u_y,"k.")
            PyPlot.plot(u_anal,y_anal,"bx-")
            PyPlot.xlabel("U Velocity (m/s)")
            PyPlot.ylabel("y-position (m)")
            # PyPlot.savefig("./figures/$figname.pdf",transparent = true)
            PyPlot.pause(0.0000005)

            spl_U = Dierckx.Spline1D(u_y,u_val[:,end])
            y_anal2 = linspace(Bot_pos,Top_pos,100)

            ell = Outlet_pos-Inlet_pos
            D_H = (Top_pos-Bot_pos)*2
            ReD = rho/mu*U0*D_H
            eff = 24/ReD
            DP = eff*4*ell/D_H*(0.5*rho*U0^2)
            DPDX = DP/ell
            u_anal2 = DPDX/(2*mu)*y_anal2.*(Top_pos-y_anal2)
            for i = 1:length(u_anal)
                if u_anal[i]<0
                    u_anal[i] = 0.0
                end
            end

            u_val_for_error = spl_U(y_anal2)
            total_error = sum(abs.((u_anal2-u_val_for_error)))/sum(u_anal2)


            figname = "error_$grid_i"
            println(total_error)
            figure("error_$grid_i")
            PyPlot.semilogy(iter,total_error*100,"k.")
            PyPlot.xlabel("Iteration (#)")
            PyPlot.ylabel("Integral Error (%)")
            # PyPlot.savefig("./figures/$figname.pdf",transparent = true)

        end
        iter+=1
    end

    y_anal = u_y
    # ell = linspace(Inlet_pos,Outlet_pos,4)
    # for i = 1:length(ell)
    ell = Outlet_pos-Inlet_pos
    D_H = (Top_pos-Bot_pos)*2
    ReD = rho/mu*U0*D_H
    eff = 24/ReD
    DP = eff*4*ell/D_H*(0.5*rho*U0^2)
    DPDX = DP/ell
    u_anal = DPDX/(2*mu)*y_anal.*(Top_pos-y_anal)
    for i = 1:length(u_anal)
        if u_anal[i]<0
            u_anal[i] = 0.0
        end
    end

    rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.9)

    # G = (mean(P_val[:,2]-P_val[:,end]))/(u_x[end-1]-u_x[2])
    figname = "outlet_$grid_i"
    figure("outlet_$grid_i")
    PyPlot.plot(u_val[:,end]*1E4,u_y,"k.",label = "Numerical U-Velocity")
    PyPlot.plot(u_anal*1E4,y_anal,"bx-",label = "Analytical U-Velocity")
    PyPlot.xlabel("U-Velocity *1E4 (m/s)")
    PyPlot.ylabel("y-position *1E4 (m)")
    PyPlot.legend(loc = "best")
    PyPlot.savefig("./figures/$figname.pdf",transparent = true)


    spl_U = Dierckx.Spline1D(u_y,u_val[:,end])
    y_anal2 = linspace(Bot_pos,Top_pos,100)

    ell = Outlet_pos-Inlet_pos
    D_H = (Top_pos-Bot_pos)*2
    ReD = rho/mu*U0*D_H
    eff = 24/ReD
    DP = eff*4*ell/D_H*(0.5*rho*U0^2)
    DPDX = DP/ell
    u_anal2 = DPDX/(2*mu)*y_anal2.*(Top_pos-y_anal2)
    for i = 1:length(u_anal)
        if u_anal[i]<0
            u_anal[i] = 0.0
        end
    end

    u_val_for_error = spl_U(y_anal2)
    total_error = sum(abs.((u_anal2-u_val_for_error)))/sum(u_anal2)

    figname = "error_$grid_i"
    println(total_error)
    figure("error_$grid_i")
    PyPlot.semilogy(iter,total_error*100,"k.")
    PyPlot.xlabel("Iteration (#)")
    PyPlot.ylabel("Integral Error (%)")
    PyPlot.savefig("./figures/$figname.pdf",transparent = true)

    rc("figure.subplot", left=0.17, bottom=0.3, top=0.97, right=.9)
    # CONTOUR PLOTS
    figname = "u_contour$grid_i"
    figure(figname)
    CS = PyPlot.contourf(Xu,Yu,u_val*1E4,interpolation = "spline16",origin = "lower",cmap = PyPlot.cm[:viridis],levels = linspace(minimum(u_val*1E4),maximum(u_val*1E4),1000))
    PyPlot.colorbar(orientation = "horizontal",extend = "both",label = "U Velocity x 1E4 (m/s)",pad = 0.3,aspect = 30,ticks = round.(linspace(minimum(u_val*1E4),maximum(u_val*1E4),8),1))
    PyPlot.xlabel("X-Position (m)")
    PyPlot.ylabel("Y-Position (m)")
    PyPlot.savefig("./figures/$figname.pdf",transparent = true)

    figname = "v_contour$grid_i"
    figure(figname)
    CS = PyPlot.contourf(Xv,Yv,v_val*1E4,interpolation = "spline16",origin = "lower",cmap = PyPlot.cm[:viridis],levels = linspace(minimum(v_val*1E4),maximum(v_val*1E4),1000))
    PyPlot.colorbar(orientation = "horizontal",extend = "both",label = "V Velocity x 1E4 (m/s)",pad = 0.3,aspect = 30,ticks = round.(linspace(minimum(v_val*1E4),maximum(v_val*1E4),8),3))
    PyPlot.xlabel("X-Position (m)")
    PyPlot.ylabel("Y-Position (m)")
    PyPlot.savefig("./figures/$figname.pdf",transparent = true)


    figname = "p_contour$grid_i"
    figure(figname)
    CS = PyPlot.contourf(Xp,Yp,P_val,interpolation = "spline16",origin = "lower",cmap = PyPlot.cm[:viridis],levels = linspace(minimum(P_val),maximum(P_val),1000))
    PyPlot.colorbar(orientation = "horizontal",extend = "both",label = "Pressure (Pa)",pad = 0.3,aspect = 30,ticks = round.(linspace(minimum(P_val),maximum(P_val),8),3))
    PyPlot.xlabel("X-Position (m)")
    PyPlot.ylabel("Y-Position (m)")
    PyPlot.savefig("./figures/$figname.pdf",transparent = true)


    figname = "stream$grid_i"
    speed = sqrt.(u_val.^2 + v_val.^2)
    lw = 1 * speed / maximum(speed) # Line Widths
    figure(figname)
    CS = PyPlot.streamplot(Xp,Yp,u_val,-v_val,color = P_val,linewidth=lw,cmap = PyPlot.cm[:viridis],arrowsize = 0.000001)
    PyPlot.plot([Inlet_pos,Outlet_pos,Outlet_pos,Inlet_pos,Inlet_pos],[Bot_pos,Bot_pos,Top_pos,Top_pos,Bot_pos],"k")
    PyPlot.xlabel("X-Position (m)")
    PyPlot.ylabel("Y-Position (m)")
    PyPlot.colorbar(orientation = "horizontal",extend = "both",label = "Pressure (Pa)",pad = 0.3,aspect = 30)
    PyPlot.savefig("./figures/$figname.pdf",transparent = true)

    PyPlot.pause(0.0000000005)

    final_iter[grid_i] = iter
    end_error[grid_i] = total_error
end

# rc("figure.subplot", left=0.17, bottom=0.18, top=0.97, right=.9)
# figname = "GridCon_Error"
# PyPlot.figure(figname)
# PyPlot.plot(N_nodes,end_error*100,"k.-")
# PyPlot.xlabel("Nodes (#)")
# PyPlot.ylabel("Converged Integral Error (%)")
# PyPlot.savefig("./figures/$figname.pdf",transparent = true)
#
# figname = "GridCon_Iters"
# PyPlot.figure(figname)
# PyPlot.plot(N_nodes,final_iter,"k.-")
# PyPlot.xlabel("Nodes (#)")
# PyPlot.ylabel("Requred Iterations(#)")
# PyPlot.savefig("./figures/$figname.pdf",transparent = true)
