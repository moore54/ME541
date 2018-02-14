# Kevin Moore 1/18/18 HW1 prob 7

using PyPlot

rc("figure", figsize=(3.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.95)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])


function fun(x)
    return sin(x)
end


function curvature(a,b,c,dx)
    d2dx2 = (a-2*b+c)/(dx^2)
    return d2dx2
end

dx_array = [.2,.1,0.5,0.025]
x = 1.0

#Case 1 Central Diff
d2dx2_central = zeros(dx_array)
for i = 1:length(dx_array)
    # i=1
    dx = dx_array[i]
    a = fun(x-dx)
    b = fun(x)
    c = fun(x+dx)

    d2dx2_central[i] = curvature(a,b,c,dx)
end

#Case 2 Backwards Diff
d2dx2_backward = zeros(dx_array)
for i = 1:length(dx_array)
    dx = dx_array[i]
    a = fun(x-2*dx)
    b = fun(x-dx)
    c = fun(x)

    d2dx2_backward[i] = curvature(a,b,c,dx)
end

d2dx2_real = ones(dx_array)*-fun(x)

PyPlot.close("all")
PyPlot.figure("Central_Reverse_DiffComparison")
PyPlot.plot(dx_array,d2dx2_real,"o",label = "Exact")
PyPlot.plot(dx_array,d2dx2_central,"o",label = "Central")
PyPlot.plot(dx_array,d2dx2_backward,"o",label = "Backward")
PyPlot.xlabel("dx")
PyPlot.ylabel("d2dx2")
PyPlot.legend(loc = "best")


d2dx2_central_error = (d2dx2_real-d2dx2_central)./d2dx2_real
d2dx2_backward_error = (d2dx2_real-d2dx2_backward)./d2dx2_real

PyPlot.figure("Central_Reverse_DiffError")
PyPlot.semilogy(dx_array,d2dx2_central_error*100.0,"o",label = "Central")
PyPlot.semilogy(dx_array,d2dx2_backward_error*100.0,"o",label = "Backward")
PyPlot.xlabel("dx")
PyPlot.ylabel("Error (%)")
PyPlot.legend(loc = "best")
