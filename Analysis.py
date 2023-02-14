import numpy as np
from numpy.linalg import inv
import sympy as sp
import matplotlib.pyplot as plt
from symbolic import Kk, k_alpha_1_2, \
    k_alpha_2_3, k_alpha_3_4, k_alpha_4_1, r_beta_1_2, r_beta_2_3, r_beta_3_4, r_beta_4_1, N, rq, u
from matplotlib.cm import ScalarMappable

fast_Kk = np.vectorize(sp.lambdify(('a', 'b', 'k1'), Kk, 'numpy'))

k_alpha_1_2 = np.vectorize(sp.lambdify(('a', 'b', 'h'), k_alpha_1_2, 'numpy'))
k_alpha_2_3 = np.vectorize(sp.lambdify(('a', 'b', 'h'), k_alpha_2_3, 'numpy'))
k_alpha_3_4 = np.vectorize(sp.lambdify(('a', 'b', 'h'), k_alpha_3_4, 'numpy'))
k_alpha_4_1 = np.vectorize(sp.lambdify(('a', 'b', 'h'), k_alpha_4_1, 'numpy'))

r_beta_1_2 = np.vectorize(sp.lambdify(('a', 'b', 'Tot', 'h'), r_beta_1_2, 'numpy'))
r_beta_2_3 = np.vectorize(sp.lambdify(('a', 'b', 'Tot', 'h'), r_beta_2_3, 'numpy'))
r_beta_3_4 = np.vectorize(sp.lambdify(('a', 'b', 'Tot', 'h'), r_beta_3_4, 'numpy'))
r_beta_4_1 = np.vectorize(sp.lambdify(('a', 'b', 'Tot', 'h'), r_beta_4_1, 'numpy'))

rq = np.vectorize(sp.lambdify(('a', 'b', 'q'), rq, 'numpy'))


class Line:
    def __init__(self, p1, p2):
        self.p1 = p1
        self.p2 = p2


class Temperature_bc:
    def __init__(self, temp, line):
        self.temp = temp
        self.line = line


class Convection_bc:
    def __init__(self, h, temp, line):
        self.h = h
        self.temp = temp
        self.line = line


class Heat_flow_bc:
    def __init__(self, q, line):
        self.q = q
        self.line = line


class Geometry:
    def __init__(self, lines):
        self.lines = lines

    def is_point_in_geometry(self, point):
        cnt = 0
        for line in self.lines:
            x1, y1 = line.p1
            x2, y2 = line.p2

            if np.max([x1, x2]) >= point[0] >= np.min([x1, x2]):
                if np.max([y1, y2]) >= point[1] >= np.min([y1, y2]):
                    if x1 == x2:
                        if point[0] == x1:
                            return True
                    else:
                        k = (y2 - y1) / (x2 - x1)
                        b = y1 - k * x1
                        if point[1] == k * point[0] + b:
                            return True

            if min(y1, y2) < point[1] <= max(y1, y2) and point[0] <= max(x1, x2):
                if y1 != y2:
                    x_intersection = (point[1] - y1) * (x2 - x1) / (y2 - y1) + x1
                else:
                    x_intersection = x1
                if x1 == x2 or point[0] <= x_intersection:
                    cnt += 1

        return cnt % 2 == 1

    def mesh(self, size):
        minX = np.min(
            [np.min(np.array([r.p1[0] for r in self.lines])), np.min(np.array([r.p2[0] for r in self.lines]))])
        maxX = np.max(
            [np.max(np.array([r.p1[0] for r in self.lines])), np.max(np.array([r.p2[0] for r in self.lines]))])

        minY = np.min(
            [np.min(np.array([r.p1[1] for r in self.lines])), np.min(np.array([r.p2[1] for r in self.lines]))])
        maxY = np.max(
            [np.max(np.array([r.p1[1] for r in self.lines])), np.max(np.array([r.p2[1] for r in self.lines]))])

        mesh_array = np.zeros(
            (int(abs(minY) / size) + int(abs(maxY) / size), int(abs(minX) / size) + int(abs(maxX) / size), 8))
        elementsData = []
        print('Tworzenie siatki elementów skończonych...')
        for i in np.arange(int(minY / size), int(maxY / size), 1):
            for j in np.arange(int(minX / size), int(maxX / size), 1):
                x = j * size + size
                y = i * size + size
                if self.is_point_in_geometry((x, y)) and self.is_point_in_geometry((x, y - size)) and \
                        self.is_point_in_geometry((x - size, y - size)) and self.is_point_in_geometry((x - size, y)):
                    if mesh_array[i - 1, j, -2] != 0:
                        p1 = mesh_array[i - 1, j, -2]
                        p2 = mesh_array[i - 1, j, -3]
                        p3 = mesh_array[i - 1, j, -4]
                    elif mesh_array[i, j - 1, 2] != 0:
                        p1 = mesh_array[i, j - 1, 2]
                        p2 = np.max(mesh_array) + 1
                        p3 = p2 + 1
                    else:
                        p1 = np.max(mesh_array) + 1
                        p2 = p1 + 1
                        p3 = p2 + 1

                    mesh_array[i, j, 0] = p1
                    mesh_array[i, j, 1] = p2
                    mesh_array[i, j, 2] = p3

                    p4 = np.max(mesh_array) + 1
                    p5 = p4 + 1
                    p6 = p5 + 1

                    mesh_array[i, j, 3] = p4
                    mesh_array[i, j, 4] = p5
                    mesh_array[i, j, 5] = p6

                    if mesh_array[i, j - 1, -2] != 0:
                        p7 = mesh_array[i, j - 1, 4]
                        p8 = mesh_array[i, j - 1, 3]
                    else:
                        p7 = np.max(mesh_array) + 1
                        p8 = p7 + 1

                    mesh_array[i, j, 6] = p7
                    mesh_array[i, j, 7] = p8

                    structure = [p1, p2, p3, p4, p5, p6, p7, p8]
                    structure = [int(x) for x in structure]
                    key = {'structure': structure, 'x': j * size, 'y': i * size}
                    elementsData.append(key)
        print('Dyskretyzacja zakończona!!')
        return elementsData, int(np.max(mesh_array))


class Element:
    def __init__(self, k, a, structure, x, y):
        self.k = k
        self.a = 0.5 * a
        self.x = x
        self.y = y
        self.structure = structure
        self.K_local = fast_Kk(self.a, self.a, self.k)
        self.rbeta = 0
        self.rq = 0


class Analysis:
    def __init__(self, element_size, k, geometry, boundary_conditions):
        self.k = k
        self.element_size = element_size
        self.geometry = geometry
        self.boundary_conditions = boundary_conditions
        self.nodes_number = 0

        self.elements = []

        self.K_matrix = None
        self.C_matrix = None

        self.nodes_temps = None
        self.bc_temps = []
        self.bc_heat = np.array([])

        self.meshing()
        self.set_bc()
        self.solution()

    # Nakładanie siatki
    def meshing(self):
        data, self.nodes_number = self.geometry.mesh(self.element_size)
        print('Nakładanie warunków brzegowych...')
        for element in data:
            e = Element(self.k, self.element_size * 0.001, element['structure'], element['x'], element['y'])

            p1 = element['x'], element['y']
            p2 = element['x'] + 0.5 * self.element_size, element['y']
            p3 = element['x'] + self.element_size, element['y']
            p4 = element['x'] + self.element_size, element['y'] + 0.5 * self.element_size
            p5 = element['x'] + self.element_size, element['y'] + self.element_size
            p6 = element['x'] + 0.5 * self.element_size, element['y'] + self.element_size
            p7 = element['x'], element['y'] + self.element_size
            p8 = element['x'], element['y'] + 0.5 * self.element_size

            points = [p1, p2, p3, p4, p5, p6, p7, p8]
            for bc in self.boundary_conditions:
                added = False
                nodes = []
                x1, y1 = bc.line.p1
                x2, y2 = bc.line.p2
                for p in range(0, len(points)):
                    if np.max([x1, x2]) >= points[p][0] >= np.min([x1, x2]):
                        if np.max([y1, y2]) >= points[p][1] >= np.min([y1, y2]):
                            if x1 == x2:
                                if points[p][0] == x1:
                                    nodes.append(element['structure'][p])
                                    added = True
                            else:
                                k = (y2 - y1) / (x2 - x1)
                                b = y1 - k * x1
                                if points[p][1] == k * points[p][0] + b:
                                    nodes.append(element['structure'][p])
                                    added = True
                if added:
                    if type(bc) == Temperature_bc:
                        key = {'nodes': nodes, 'temp': bc.temp}
                        self.bc_temps.append(key)
                    if type(bc) == Convection_bc:
                        if len(nodes) == 3:
                            if nodes[0] == element['structure'][0] and nodes[1] == element['structure'][1] \
                                    and nodes[2] == element['structure'][2]:
                                e.K_local = e.K_local + k_alpha_1_2(e.a, e.a, bc.h)
                                e.rbeta = e.rbeta + r_beta_1_2(e.a, e.a, bc.temp, bc.h)

                            elif nodes[0] == element['structure'][2] and nodes[1] == element['structure'][3] \
                                    and nodes[2] == element['structure'][4]:
                                e.K_local = e.K_local + k_alpha_2_3(e.a, e.a, bc.h)
                                e.rbeta = e.rbeta + r_beta_2_3(e.a, e.a, bc.temp, bc.h)

                            elif nodes[0] == element['structure'][4] and nodes[1] == element['structure'][5] \
                                    and nodes[2] == element['structure'][6]:
                                e.K_local = e.K_local + k_alpha_3_4(e.a, e.a, bc.h)
                                e.rbeta = e.rbeta + r_beta_3_4(e.a, e.a, bc.temp, bc.h)

                            elif nodes[0] == element['structure'][0] and nodes[1] == element['structure'][-2] \
                                    and nodes[2] == element['structure'][-1]:
                                e.K_local = e.K_local + k_alpha_4_1(e.a, e.a, bc.h)
                                e.rbeta = e.rbeta + r_beta_4_1(e.a, e.a, bc.temp, bc.h)
                    if type(bc) == Heat_flow_bc:
                        e.rq = e.rq + rq(e.a, e.a, bc.q)
            self.elements.append(e)
        print('Nakładanie warunków brzegowych zakończone!!!')

    def set_bc(self):
        # Utworzenie pustych macierzy sztywności i wymuszeń
        self.K_matrix = np.zeros((self.nodes_number, self.nodes_number))
        self.C_matrix = np.zeros((self.nodes_number, 1))

        # Obliczenie macierzy globalnej
        print('Agregacja macierzy...')
        for e in self.elements:

            for k in range(0, 8):
                for i in range(0, 8):
                    q = e.K_local[k, i]
                    x = int(e.structure[k])
                    y = int(e.structure[i])
                    self.K_matrix[x - 1, y - 1] = self.K_matrix[x - 1, y - 1] + q

            if type(e.rq) != int:
                for i in range(0, 8):
                    self.C_matrix[e.structure[i] - 1, 0] = self.C_matrix[e.structure[i] - 1, 0] + e.rq[0, i]

            # Agragacja macierzy Rbeta
            if type(e.rbeta) != int:
                for i in range(0, 8):
                    self.C_matrix[e.structure[i] - 1, 0] = self.C_matrix[e.structure[i] - 1, 0] + e.rbeta[0, i]

        # Nadanie warunku brzegowego temperatury
        for bc in self.bc_temps:
            for p in bc['nodes']:
                self.K_matrix[p - 1, :] = 0
                self.C_matrix[p - 1, 0] = bc['temp']

        # Zerowanie głównej przekątnej macierzy sztywności
        for i in range(0, self.nodes_number):
            if self.K_matrix[i, i] == 0:
                self.K_matrix[i, i] = 1

        print('Agregacja zakończona!!')

    def solution(self):
        # Obliczenie temperatur w węzłach
        print('Rozwiązywanie równania...')
        self.nodes_temps = np.linalg.solve(self.K_matrix, self.C_matrix)
        #self.nodes_temps = inv(self.K_matrix).dot(self.C_matrix)
        print('Rozwiązywanie gotowe!!')

    # Informacje o siatce
    def mesh_info(self):
        print(f"Liczba elementów: {len(self.elements)}")
        print(f"Liczba węzłów: {self.nodes_number}")
        print(f"Rozmiar elementu: {self.element_size} [mm]")

    def plot_results(self):
        fun1 = N.dot(u.transpose())
        fun1 = np.vectorize(sp.lambdify((
            'u_1', 'u_2', 'u_3', 'u_4', 'u_5', 'u_6', 'u_7', 'u_8', 'a', 'b', 's', 't'), fun1, 'numpy'))
        cmap = plt.get_cmap('jet')
        norm = plt.Normalize(vmin=np.min(self.nodes_temps), vmax=np.max(self.nodes_temps))
        sm = ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        for e in self.elements:
            s = np.linspace(e.x * 0.001, e.x * 0.001 + e.a * 2, int(self.element_size / 0.1))
            t = np.linspace(e.y * 0.001, e.y * 0.001 + e.a * 2, int(self.element_size / 0.1))

            xs = e.x * 0.001 + e.a
            ys = e.y * 0.001 + e.a

            ss, tt = np.meshgrid(s, t)

            nodes = np.zeros((8, 1))
            for k in range(0, 8):
                nodes[k, 0] = self.nodes_temps[e.structure[k] - 1, 0]

            temps = fun1(
                nodes[0, 0], nodes[1, 0], nodes[2, 0], nodes[3, 0], nodes[4, 0], nodes[5, 0], nodes[6, 0], nodes[7, 0],
                e.a, e.a, ss - xs, tt - ys)

            plt.pcolormesh(ss * 1000, tt * 1000, temps, cmap=cmap, norm=norm)

            plt.plot([e.x, e.x + self.element_size], [e.y, e.y], c='black', linewidth=0.5)
            plt.plot([e.x + self.element_size, e.x + self.element_size], [e.y, e.y + self.element_size], c='black',
                     linewidth=0.5)
            plt.plot([e.x, e.x], [e.y, e.y + self.element_size], c='black', linewidth=0.5)
            plt.plot([e.x, e.x + self.element_size], [e.y + self.element_size, e.y + self.element_size], c='black',
                     linewidth=0.5)

        cbar = plt.colorbar(sm)
        cbar.set_label('Temperatura [K]', fontsize=14)
        cbar.ax.text(1.05, 1.0, f"Max: {round(np.max(self.nodes_temps), 2)}", transform=cbar.ax.transAxes, va='bottom',
                     ha='left')
        cbar.ax.set_xlabel(f"Min: {round(np.min(self.nodes_temps), 2)}", rotation=0, ha='left')

        plt.title('Rozkład temperatury w domenie')
        plt.xlabel('Odległość [mm]')
        plt.ylabel('Odległość [mm]')
        plt.show()
