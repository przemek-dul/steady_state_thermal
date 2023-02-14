from Analysis import *

l1 = Line(p1=[0, 0], p2=[60, 0])
l2 = Line(p1=[60, 0], p2=[60, 40])
l3 = Line(p1=[60, 40], p2=[50, 40])
l4 = Line(p1=[50, 40], p2=[50, 20])
l5 = Line(p1=[50, 20], p2=[10, 20])
l6 = Line(p1=[10, 20], p2=[10, 40])
l7 = Line(p1=[10, 40], p2=[0, 40])
l8 = Line(p1=[0, 40], p2=[0, 0])

lines = [l1, l2, l3, l4, l5, l6, l7, l8]
g = Geometry(lines)

bc1 = Temperature_bc(temp=273 + 150, line=l1)
bc2 = Convection_bc(h=100, temp=273 + 20, line=Line(p1=[60, 0], p2=[60, 10]))
bc3 = Convection_bc(h=85, temp=273 + 20, line=Line(p1=[60, 10], p2=[60, 40]))
bc4 = Convection_bc(h=85, temp=273 + 20, line=l3)
bc5 = Convection_bc(h=85, temp=273 + 20, line=l4)
bc6 = Convection_bc(h=85, temp=273 + 20, line=l5)
bc7 = Convection_bc(h=85, temp=273 + 20, line=l6)
bc8 = Convection_bc(h=85, temp=273 + 20, line=l7)
bc9 = Convection_bc(h=100, temp=273 + 20, line=Line(p1=[0, 0], p2=[0, 10]))
bc10 = Convection_bc(h=100, temp=273 + 20, line=Line(p1=[0, 10], p2=[0, 40]))

bc = [bc1, bc2, bc3, bc4, bc5, bc6, bc7, bc8, bc9, bc10]

a = Analysis(element_size=0.5, k=55, geometry=g, boundary_conditions=bc)
a.mesh_info()
a.plot_results()
