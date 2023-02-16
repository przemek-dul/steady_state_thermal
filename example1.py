from Analysis import *

l1 = Line(p1=[0, 50], p2=[0, 0])
l2 = Line(p1=[0, 0], p2=[10, 0])
l3 = Line(p1=[10, 0], p2=[10, 20])
l4 = Line(p1=[10, 20], p2=[60, 20])
l5 = Line(p1=[60, 20], p2=[60, 30])
l6 = Line(p1=[60, 30], p2=[10, 30])
l7 = Line(p1=[10, 30], p2=[10, 50])
l8 = Line(p1=[10, 50], p2=[0, 50])

lines = [l1, l2, l3, l4, l5, l6, l7, l8]
g = Geometry(lines)

bc1 = Temperature_bc(temp=273+300, line=l1)
bc2 = Convection_bc(h=85, temp=273+28, line=l3)
bc3 = Convection_bc(h=85, temp=273+28, line=l4)
bc4 = Convection_bc(h=85, temp=273+28, line=l5)
bc5 = Convection_bc(h=85, temp=273+28, line=l6)
bc6 = Convection_bc(h=85, temp=273+28, line=l7)


bc = [bc1, bc2, bc3, bc4, bc5, bc6]

a = Analysis(element_size=10, k=55, geometry=g, boundary_conditions=bc)
a.mesh_info()
a.plot_results()
