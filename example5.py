from Analysis import *

l1 = Line(p1=[0, 0], p2=[60, 0])
l2 = Line(p1=[60, 0], p2=[60, 10])
l3 = Line(p1=[60, 10], p2=[10, 10])
l4 = Line(p1=[10, 10], p2=[10, 70])
l5 = Line(p1=[10, 70], p2=[60, 70])
l6 = Line(p1=[60, 70], p2=[60, 80])
l7 = Line(p1=[60, 80], p2=[0, 80])
l8 = Line(p1=[0, 80], p2=[0, 0])

lines = [l1, l2, l3, l4, l5, l6, l7, l8]
g = Geometry(lines)

bc1 = Temperature_bc(temp=273 + 100, line=Line(p1=[0, 50], p2=[0, 30]))
bc2 = Convection_bc(h=100, temp=273 + 20, line=Line(p1=[0, 30], p2=[0, 0]))
bc3 = Convection_bc(h=100, temp=273 + 20, line=l1)
bc4 = Convection_bc(h=100, temp=273 + 20, line=l2)
bc5 = Convection_bc(h=100, temp=273 + 20, line=l3)
bc6 = Convection_bc(h=100, temp=273 + 20, line=Line(p1=[10, 10], p2=[10, 30]))
bc7 = Convection_bc(h=85, temp=273 + 20, line=Line(p1=[10, 30], p2=[10, 50]))
bc8 = Convection_bc(h=100, temp=273 + 20, line=Line(p1=[10, 50], p2=[10, 70]))
bc9 = Convection_bc(h=100, temp=273 + 20, line=l5)
bc10 = Convection_bc(h=100, temp=273 + 20, line=l6)
bc11 = Convection_bc(h=100, temp=273 + 20, line=l7)
bc12 = Convection_bc(h=100, temp=273 + 20, line=Line(p1=[0, 80], p2=[0, 50]))

bc = [bc1, bc2, bc3, bc4, bc5, bc6, bc7, bc8, bc9, bc10, bc11, bc12, ]

a = Analysis(element_size=5, k=55, geometry=g, boundary_conditions=bc)
a.mesh_info()
a.plot_results()