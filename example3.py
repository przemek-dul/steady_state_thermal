from Analysis import *

l1 = Line(p1=[0, 0], p2=[350, 0])
l2 = Line(p1=[350, 0], p2=[350, 200])
l3 = Line(p1=[350, 200], p2=[300, 200])
l4 = Line(p1=[300, 200], p2=[300, 50])
l5 = Line(p1=[300, 50], p2=[200, 50])
l6 = Line(p1=[200, 50], p2=[200, 150])
l7 = Line(p1=[200, 150], p2=[150, 150])
l8 = Line(p1=[150, 150], p2=[150, 50])
l9 = Line(p1=[150, 50], p2=[50, 50])
l10 = Line(p1=[50, 50], p2=[50, 200])
l11 = Line(p1=[50, 200], p2=[0, 200])
l12 = Line(p1=[0, 200], p2=[0, 0])


lines = [l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12]
g = Geometry(lines)

bc1 = Heat_flow_bc(q=60000, line=l1)
bc2 = Convection_bc(h=85, temp=273 + 30, line=l2)
bc3 = Convection_bc(h=85, temp=273 + 30, line=l3)
bc4 = Convection_bc(h=85, temp=273 + 30, line=l4)
bc5 = Convection_bc(h=85, temp=273 + 30, line=l5)
bc6 = Convection_bc(h=85, temp=273 + 30, line=l6)
bc7 = Convection_bc(h=85, temp=273 + 30, line=l7)
bc8 = Convection_bc(h=85, temp=273 + 30, line=l8)
bc9 = Convection_bc(h=85, temp=273 + 30, line=l9)
bc10 = Convection_bc(h=85, temp=273 + 30, line=l10)
bc11 = Convection_bc(h=85, temp=273 + 30, line=l11)
bc12 = Convection_bc(h=85, temp=273 + 30, line=l12)


bc = [bc1, bc2, bc3, bc4, bc5, bc6, bc7, bc8, bc9, bc10, bc11, bc12]

a = Analysis(element_size=10, k=55, geometry=g, boundary_conditions=bc)
a.mesh_info()
a.plot_results()