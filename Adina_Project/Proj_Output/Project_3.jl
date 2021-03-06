using PyPlot
using CSV
PyPlot.close("all")
my_mesh = [0.5,0.75,1.0,1.25,1.5,1.75,2.25,3.0,4.0,5.0]
time = [1,8,15,25,39,59,119.0,274.0,672,1391]
elements = [240,675,1440,2625,4320,6615,13365.0,30240,69120,132000]
nodes = [287,760,1573,2816,4579,6952,13888.0,31117,70609,134261]

rc("figure", figsize=(6.0, 3.0))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.1, bottom=0.18, top=0.97, right=0.82)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7", "000486","700000","006907","4C0099"])


filename = []
fig1 = "fig3b.pdf"
PyPlot.figure(fig1)
idx_MaxVor = 0
for i = [1,2,3,4,5,7,8,9,10]
# i = 9
    filename = "Mesh$(round.(my_mesh[i],2))_vorticy.txt"
    vorticity = CSV.read(filename,header = false,delim = ",")
    # println(filename)
    vor = vorticity[:,3].*0.01/0.05

    PyPlot.plot(vorticity[:,2],vor,"-",label = "$(round.(Int,elements[i]))")
end
PyPlot.xlabel("Y-Location (m)")
PyPlot.ylabel("Vorticity")
legend(loc="center left", bbox_to_anchor=(1, 0.5),title = "# Elements")
PyPlot.savefig(fig1,transparent = true)


rc("figure", figsize=(4.0, 3.0))
rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.92)
fig2 = "fig3b_time_cells.pdf"
PyPlot.figure(fig2)
PyPlot.semilogx(elements,time,"o",label = "Elements")
# PyPlot.semilogx(nodes,time,"o-",label = "Nodes")
PyPlot.xlabel("# Elements")
PyPlot.ylabel("Time to Solve (s)")
# PyPlot.legend(loc = "best")
PyPlot.savefig(fig2,transparent = true)
# legend(loc="center left", bbox_to_anchor=(1, 0.5))

filename = []
fig3 = "fig3b_max.pdf"
rc("figure", figsize=(6.0, 3.0))
rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.92)
PyPlot.figure(fig3)
vor_max = zeros(length(elements))
for i = [1,2,3,4,5,7,8,9,10]
# i = 9
    filename = "Mesh$(round.(my_mesh[i],2))_vorticy.txt"
    vorticity = CSV.read(filename,header = false,delim = ",")
    # println(filename)
    vor_max[i] = maximum(vorticity[:,3].*0.01/0.05)
end
PyPlot.semilogx(elements,vor_max,"o")
PyPlot.xlabel("# Elements")
PyPlot.ylabel("Maximum Vorticity")
# legend(loc="center left", bbox_to_anchor=(1, 0.5),title = "# Elements")
PyPlot.savefig(fig3,transparent = true)


#
# ##########
# #  Plot Time to Solve and Elements #
# ##########
# fig = figure("pyplot_multiaxis",figsize=(10,10))
# p = plot(x,y1,linestyle="-",marker="o",label="First") # Plot a basic line
# ax = gca()
# title("Multi-axis Plot")
#
# xlabel("X Axis")
# font1 = Dict("color"=>"blue")
# ylabel("Left Axis",fontdict=font1)
# setp(ax[:get_yticklabels](),color="blue") # Y Axis font formatting
#
# ################
# #  Other Axes  #
# ################
# new_position = [0.06;0.06;0.77;0.91] # Position Method 2
# ax[:set_position](new_position) # Position Method 2: Change the size and position of the axis
# #fig[:subplots_adjust](right=0.85) # Position Method 1
#
# ax2 = ax[:twinx]() # Create another axis on top of the current axis
# font2 = Dict("color"=>"purple")
# ylabel("Right Axis",fontdict=font2)
# p = plot_date(x,y2,color="purple",linestyle="-",marker="o",label="Second") # Plot a basic line
# ax2[:set_position](new_position) # Position Method 2
# setp(ax2[:get_yticklabels](),color="purple") # Y Axis font formatting
#
# ax3 = ax[:twinx]() # Create another axis on top of the current axis
# ax3[:spines]["right"][:set_position](("axes",1.12)) # Offset the y-axis label from the axis itself so it doesn't overlap the second axis
# font3 = Dict("color"=>"green")
# ylabel("Far Right Axis",fontdict=font3)
# p = plot_date(x,y3,color="green",linestyle="-",marker="o",label="Third") # Plot a basic line
# ax3[:set_position](new_position) # Position Method 2
# setp(ax3[:get_yticklabels](),color="green") # Y Axis font formatting
#
# axis("tight")
#
# # Enable just the right part of the frame
# ax3[:set_frame_on](true) # Make the entire frame visible
# ax3[:patch][:set_visible](false) # Make the patch (background) invisible so it doesn't cover up the other axes' plots
# ax3[:spines]["top"][:set_visible](false) # Hide the top edge of the axis
# ax3[:spines]["bottom"][:set_visible](false) # Hide the bottom edge of the axis
#
# fig[:canvas][:draw]() # Update the figure
