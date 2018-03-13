using CSV
using PyPlot
using Dierckx

PyPlot.close("all")

rc("figure", figsize=(6.0, 3.0))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.15, bottom=0.18, top=0.97, right=0.92)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7", "000486","700000","006907","4C0099"])
color_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7", "#000486","#700000","#006907","#4C0099"]

#-------- Validation 1 --------#
# filename = "Cylinder_graph.txt"
# file_data = CSV.read(filename,header = false,delim = ",")

# filename2 = "cylinder2.txt"
# file_data2 = CSV.read(filename2,header = false,delim = ",")

filename_other = "cylinder_other.txt"
file_data_other = CSV.read(filename_other,header = false,delim = ",")

cyl_spl = Dierckx.Spline1D(10.0.^(file_data_other[:,1]),10.0.^(file_data_other[:,2]);k=1)
Re_splined = logspace(-2,6,100)
splined_cd = cyl_spl(Re_splined)

figname = "validation"
PyPlot.figure(figname)
# PyPlot.loglog(10.0.^(file_data[:,1]),10.0.^(file_data[:,2]),".-")
# PyPlot.loglog(10.0.^(file_data2[:,1]),10.0.^(file_data2[:,2]),".-")
PyPlot.loglog(10.0.^(file_data_other[:,1]),10.0.^(file_data_other[:,2]),".",label="Cheng")
# PyPlot.loglog(Re_splined,splined_cd,"-")
PyPlot.plot(20,2.2027,"x",label = "CFD Re 20")
PyPlot.plot(150,1.32,"x",label = "CFD Re 150")
PyPlot.xlabel("Reynolds Number")
PyPlot.ylabel("Drag Coefficient")
legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)



Error_20 = (2.2027-cyl_spl(20))/cyl_spl(20)
Error_150 = (1.32-cyl_spl(150))/cyl_spl(150)


#-------- Grid Convergence --------#
V20 = 0.30286
V150 = 2.2714
rho = 1.225
mu = 1.855E-5
diam = 0.0001

cells = [3.366E3, 6.163E3,1.508E4,4.5395E4,1.526E5]
D20 = [1.2283E-5,1.23045E-5,1.2333E-5,1.2332E-5,1.23345E-5]
D150 = [0.000407748,0.000412198,0.00042197,0.00042997,0.0004323]
VC = [3.846,4.08,4.17,4.26,4.237]
D20 = D20/(0.5*rho*V20^2*diam)
D150 = D150/(0.5*rho*V150^2*diam)

rc("figure", figsize=(3.0, 3.0))
rc("figure.subplot", left=0.25, bottom=0.18, top=0.97, right=0.92)

figname = "grid20"
PyPlot.figure(figname)
PyPlot.semilogx(cells,round.(D20,3),".",label="Re 20")
PyPlot.xlabel("Cell Number")
PyPlot.ylabel("Drag Coefficient")
PyPlot.savefig("$figname.pdf",transparent = true)

figname = "grid150"
PyPlot.figure(figname)
PyPlot.semilogx(cells,round.(D150,3),".",label="Re 150")
PyPlot.xlabel("Cell Number")
PyPlot.ylabel("Drag Coefficient")
PyPlot.savefig("$figname.pdf",transparent = true)


#--------- BoundaryConversion ---------#
Diameter=[18,16,14,12,10]
Cells=[4.598E4,4.7921E4,4.9904E4,4.951E4,4.5395E4]
D20=[1.176E-5,1.1857E-5,1.1947E-5,1.2102E-5,1.2332E-5]
D150=[.0004196,.0004212,.0004231,.0004262,.00042997]
Frequency=[4.22,4.18,4.199,4.23,4.26]

D20 = D20/(0.5*rho*V20^2*diam)
D150 = D150/(0.5*rho*V150^2*diam)

rc("figure", figsize=(3.5, 3.0))
rc("figure.subplot", left=0.25, bottom=0.18, top=0.97, right=0.92)

figname = "boundary20"
PyPlot.figure(figname)
PyPlot.plot(Diameter,round.(D20,3),".",label="Re 20")
PyPlot.xlabel("Inlet Radius (diameters)")
PyPlot.ylabel("Drag Coefficient")
PyPlot.savefig("$figname.pdf",transparent = true)

figname = "boundary150"
PyPlot.figure(figname)
PyPlot.plot(Diameter,round.(D150,3),".",label="CD")
# PyPlot.plot(Diameter,Frequency,".",label = "Frequency (Hz)")
PyPlot.xlabel("Inlet Radius (diameters)")
PyPlot.ylabel("Drag Coefficient")
# PyPlot.ticklabel_format(style="sci", axis="x")
# legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)

#TimeStep
TimeStep=[.01,.005,.0025,.00125]
D150=[.00041909,.0004212,.0004217,.000421698]
Frequency=[4.098,4.18,4.42,4.24]

D150 = D150/(0.5*rho*V150^2*diam)

figname = "timestep150"
PyPlot.figure(figname)
PyPlot.plot(TimeStep,round.(D150,3),".",label="CD")
# PyPlot.plot(Diameter,Frequency,".",label = "Frequency (Hz)")
PyPlot.xlabel("Time Step (s)")
PyPlot.ylabel("Drag Coefficient")
# PyPlot.ticklabel_format(style="sci", axis="x")
# legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)

figname = "frequencyTS150"
PyPlot.figure(figname)
PyPlot.plot(TimeStep,round.(Frequency,3),".",label="CD")
# PyPlot.plot(Diameter,Frequency,".",label = "Frequency (Hz)")
PyPlot.xlabel("Time Step (s)")
PyPlot.ylabel("Von Karman Frequency (Hz)")
# PyPlot.ticklabel_format(style="sci", axis="x")
# legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)

#InnerIterations
InnerIterations=[8,6,4,2]
D150=[.00042176,.00042176,.000421707,.0004205]
Frequency=[4.195,4.43,4.21,4.18]

D150 = D150/(0.5*rho*V150^2*diam)


figname = "II150"
PyPlot.figure(figname)
PyPlot.plot(InnerIterations,round.(D150,3),".",label="CD")
# PyPlot.plot(Diameter,Frequency,".",label = "Frequency (Hz)")
PyPlot.xlabel("Inner Iterations (#)")
PyPlot.ylabel("Drag Coefficient")
# PyPlot.ticklabel_format(style="sci", axis="x")
# legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)

figname = "frequencyII150"
PyPlot.figure(figname)
PyPlot.plot(TimeStep,round.(Frequency,3),".",label="CD")
# PyPlot.plot(Diameter,Frequency,".",label = "Frequency (Hz)")
PyPlot.xlabel("Inner Iterations (#)")
PyPlot.ylabel("Von Karman Frequency (Hz)")
# PyPlot.ticklabel_format(style="sci", axis="x")
# legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)

#------ Time Histories ------#

F20 = CSV.read("./force_history_Re20.csv",header = true,delim = ",")
F150 = CSV.read("./force_history.csv",header = true,delim = ",")
Vor20 = CSV.read("./average_vorticity_Re20.csv",header = true,delim = ",")
Vor150 = CSV.read("./average_vorticity.csv",header = true,delim = ",")

D20 = F20[:,2]/(0.5*rho*V20^2*diam)
D150 = F150[:,2]/(0.5*rho*V150^2*diam)

rc("figure", figsize=(3.9, 3.0))
rc("figure.subplot", left=0.20, bottom=0.18, top=0.97, right=0.98)

figname = "CD_time_20"
PyPlot.figure(figname)
PyPlot.plot(F20[:,1],round.(D20,3),"-",label="Re 20")
# PyPlot.plot(F150[:,1],round.(D150,3),"-",label = "Re 150")
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Drag Coefficient")
PyPlot.ylim([2.1,3.0])
# PyPlot.ticklabel_format(style="sci", axis="x")
# legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)

figname = "CD_time_150"
PyPlot.figure(figname)
# PyPlot.plot(F20[:,1],round.(D20,3),"-",label="Re 20")
PyPlot.plot(F150[:,1],round.(D150,3),"-",label = "Re 150")
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Drag Coefficient")
PyPlot.ylim([1.0,1.4])
# PyPlot.ticklabel_format(style="sci", axis="x")
# legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)

figname = "Ave_20"
PyPlot.figure(figname)
PyPlot.plot(Vor20[:,1],round.(Vor20[:,2],3),"-",label="Re 20")
# PyPlot.plot(Vor150[:,1],round.(Vor150[:,2],3),"-",label = "Re 150")
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Average Domain Vorticity (/s)")
# PyPlot.ylim([2.1,3.0])
# PyPlot.ticklabel_format(style="sci", axis="x")
# legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)

figname = "Ave_150"
PyPlot.figure(figname)
# PyPlot.plot(Vor20[:,1],round.(Vor20[:,2],3),"-",label="Re 20")
PyPlot.plot(Vor150[:,1],round.(Vor150[:,2],3),"-",label = "Re 150")
PyPlot.xlabel("Time (s)")
PyPlot.ylabel("Average Domain Vorticity (/s)")
# PyPlot.ylim([1.0,1.4])
# PyPlot.ticklabel_format(style="sci", axis="x")
# legend(loc="best")
PyPlot.savefig("$figname.pdf",transparent = true)
