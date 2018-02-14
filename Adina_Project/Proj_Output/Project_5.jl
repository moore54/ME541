using PyPlot
using CSV
PyPlot.close("all")
my_B = [1.0,25.0,100.0,129.0,130.0,131.0,132.0,140.0,150.0]
my_B2 = ["1.0","25.0","100.0","129.0","130.0","131.0","132.0","140.0","150.0"]




rc("figure", figsize=(6.0, 3.0))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.15, bottom=0.18, top=0.97, right=0.82)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7", "000486","700000","006907","4C0099"])
color_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7", "#000486","#700000","#006907","#4C0099"]

filename = []
# #,Y-pos,Vorticity, Displacement, Pressure
defl=zeros(841,length(my_B))
press=zeros(841,length(my_B))
vor=zeros(841,length(my_B))
file_data = []
for i = 1:length(my_B)
# i = 9
    #Load Adina Results
    filename = "B$((my_B2[i]))_mesh3.0_TOL0.0001.txt"
    file_data = CSV.read(filename,header = false,delim = ",")
    # println(filename)
    defl[:,i] = file_data[:,4]
    press[:,i] = file_data[:,5]
    vor[:,i] = file_data[:,3]#.*0.01/0.05
end


fig1 = "fig5_vor.pdf"
PyPlot.figure(fig1)
for i = 1:length(my_B)
    PyPlot.plot(file_data[:,2],vor[:,i],"-",color = color_cycle[i],label = "$(round.(my_B[i],6))")

    if my_B[i]<=129
        filename_paper= "B$(convert(Int,my_B[i]))"
        file_data_paper = CSV.read("$(filename_paper)C.txt",header = false,delim = "\t")
        PyPlot.plot(file_data_paper[:,1]/100,file_data_paper[:,2],".",color = color_cycle[i])
    end
end
PyPlot.xlabel("Y-Location (m)")
PyPlot.ylabel("Vorticity")
legend(loc="center left", bbox_to_anchor=(1, 0.5),title = L"\beta")
PyPlot.savefig(fig1,transparent = true)

fig1 = "fig5_press.pdf"
PyPlot.figure(fig1)
for i = 1:length(my_B)
    PyPlot.plot(file_data[:,2],press[:,i]-.93,"-",color = color_cycle[i],label = "$(round.(my_B[i],6))")

    if my_B[i]<=129
        filename_paper= "B$(convert(Int,my_B[i]))"
        file_data_paper = CSV.read("$(filename_paper)B.txt",header = false,delim = "\t")
        PyPlot.plot((file_data_paper[:,1]/100+1500)/10000,file_data_paper[:,2]/100000,".",color = color_cycle[i])
    end
end
PyPlot.xlabel("Y-Location (m)")
PyPlot.ylabel("Pressure")
legend(loc="center left", bbox_to_anchor=(1, 0.5),title = L"\beta")
PyPlot.savefig(fig1,transparent = true)


fig1 = "fig5_def.pdf"
PyPlot.figure(fig1)
for i = 1:length(my_B)
    PyPlot.plot(file_data[:,2],defl[:,i],"-",color = color_cycle[i],label = "$(round.(my_B[i],6))")

    if my_B[i]<=129
        filename_paper= "B$(convert(Int,my_B[i]))"
        file_data_paper = CSV.read("$(filename_paper)A.txt",header = false,delim = "\t")
        PyPlot.plot(file_data_paper[:,1]/100,file_data_paper[:,2]/100,".",color = color_cycle[i])
    end
end
PyPlot.xlabel("Y-Location (m)")
PyPlot.ylabel("Z-Position (m)")
legend(loc="center left", bbox_to_anchor=(1, 0.5),title = L"\beta")
PyPlot.savefig(fig1,transparent = true)
