A_mat = [5.96283
-0.1746
0
-5.707
0
0
0
0
0
0
0
0
-0.1746
6.15476
-0.1746
0
-5.7275
0
0
0
0
0
0
0
0
-0.1746
6.15708
0
0
-5.7298
0
0
0
0
0
0
-5.707
0
0
11.5226
-0.2208
0
-5.5916
0
0
0
0
0
0
-5.7275
0
-0.2208
11.775
-0.2208
0
-5.6059
0
0
0
0
0
0
-5.7298
0
-0.2208
11.7796
0
0
-5.6081
0
0
0
0
0
0
-5.5916
0
0
11.5073
-0.2208
0
-5.6917
0
0
0
0
0
0
-5.6059
0
-0.2208
11.7574
-0.2208
0
-5.7099
0
0
0
0
0
0
-5.6081
0
-0.2208
11.7617
0
0
-5.7119
0
0
0
0
0
0
-5.6917
0
0
5.93618
-0.1632
0
0
0
0
0
0
0
0
-5.7099
0
-0.1632
6.11437
-0.1632
0
0
0
0
0
0
0
0
-5.7119
0
-0.1632
6.11644]

A2 = zeros(12,12)
j = 1
for i = 1:12
    A2[i,:] = A_mat[j:j+11]
    j = j+12
end


b_arr = [-7.84E-05
-0.000837921
-0.000901208
3.38E-05
-0.000156806
-0.000194101
0.000278606
0.000173685
0.000149205
0.001535294
0.001006104
0.000967065]


p_cor_new = A2\b_arr

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
p_relax = 0.5


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

p_cor = zeros(P_val)
dw_n = zeros(P_val)
ds_n = zeros(P_val)
j = 1
for i = 1:length(P_val[:,1])-2
    p_cor[i+1,2:end-2] = p_cor_new[j:j+length(p_cor[1,:])-4]
    j=j+length(p_cor[1,:])-3
end
