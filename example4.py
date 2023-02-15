from Analysis import *

l1 = Line(p1=[0, 0], p2=[200, 0])
l2 = Line(p1=[200, 0], p2=[200, 150])
l3 = Line(p1=[200, 150], p2=[400, 150])
l4 = Line(p1=[400, 150], p2=[400, 300])
l5 = Line(p1=[400, 300], p2=[0, 300])
l6 = Line(p1=[0, 300], p2=[0, 0])

lines = [l1, l2, l3, l4, l5, l6]
g = Geometry(lines)

bc1 = Temperature_bc(temp=273 + 100, line=l1)
bc2 = Temperature_bc(temp=273 + 30, line=l4)
bc3 = Convection_bc(h=85, temp=273 + 28, line=l5)
bc4 = Convection_bc(h=85, temp=273 + 28, line=l6)

bc = [bc1, bc2, bc3, bc4]

a = Analysis(element_size=10, k=55, geometry=g, boundary_conditions=bc)
a.mesh_info()
a.plot_results()