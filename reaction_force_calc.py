import numpy as np
import math
import pyvista as pv
import time


class Force_analysis:
    def __init__(self, mesh_size, mesh_lx, mesh_ly, mesh_lz, point5, point6, point7,
                 point8, point9, point10, force_point, force, E, v):

        ### mesh info
        self.mesh_size =float(mesh_size)
        self.mesh_lx = float(mesh_lx)
        self.mesh_ly = float(mesh_ly)
        self.mesh_lz = float(mesh_lz)
        
        #### Dirichelet B.C
        self.point1 = [0.0,0.0,0.0]
        self.point2 = [self.mesh_lx,0.0,0.0]
        self.point3 = [0.0,0.0,self.mesh_lz]
        self.point4 = [self.mesh_lx,0.0,self.mesh_lz]
        self.point5 = point5
        self.point6 = point6
        self.point7 = point7
        self.point8 = point8
        self.point9 = point9
        self.point10 = point10
        
        
        ### material
        self.E = float(E)
        self.v = float(v)

        ### Neuman B.C
        self.force_point = force_point
        self.force = force
        
        
        ### derived
        self.mesh_ex = round((self.mesh_lx/self.mesh_size)*1000)  
        self.mesh_ey = round((self.mesh_ly/self.mesh_size)*1000)    
        self.mesh_ez = round((self.mesh_lz/self.mesh_size)*1000)  

        self.mesh_nx = self.mesh_ex + 1
        self.mesh_ny = self.mesh_ey + 1
        self.mesh_nz = self.mesh_ez + 1
        self.num_nodes = self.mesh_nx * self.mesh_ny * self.mesh_nz
        self.num_elements = self.mesh_ex * self.mesh_ey * self.mesh_ez
        
        ### mesh creation

        self.mesh_point_x, self.mesh_point_y, self.mesh_point_z, self.mesh_point_x_ratio, self.mesh_point_y_ratio, self.mesh_point_z_ratio = self.mesh_ratio(self.point1, self.point2, self.point3, self.point4, 
                                                                                                                                                             self.point5, self.point6, self.point7, self.point8, self.point9, 
                                                                                                                                                             self.point10, self.force_point, self.mesh_nx, self.mesh_ny, 
                                                                                                                                                             self.mesh_nz, self.mesh_lx, self.mesh_ly, self.mesh_lz)


        self.nodes = self.mesh_creation(self.mesh_point_x, self.mesh_point_y, self.mesh_point_z, self.mesh_point_x_ratio, 
                              self.mesh_point_y_ratio, self.mesh_point_z_ratio)

        

        self.mesh_plot(self.nodes, self.point1, self.point2, self.point3, self.point4, self.point5, self.point6, self.point7, self.point8, self.point9, self.point10, self.force_point)


    def run_simulation(self):
        ### Global matrix K creatoin
        self.conn = self.create_connectivity()
        self.C = self.c_matrix(self.E, self.v)
        self.K = self.global_stiffness_matrix()

        
        ### Assing boundary condition
        print('assign nodal forces and boundary conditions')
        self.f = np.zeros((3*self.num_nodes))
        self.u = np.zeros((3*self.num_nodes))
        self.node_filter, self.node_filter_reverse, self.point_dic, self.point_order = self.apply_boundary_conditions(self.nodes, self.num_nodes, self.point5, self.point6, self.point7, self.point8, 
                                       self.point9, self.point10, self.force_point, self.u, self.f, self.force)
        

        ### KPP, KPu, Up, Fp
        self.Kpp = self.Kpp_indexing(self.K, self.node_filter_reverse, self.node_filter)
        self.Kpu = self.Kpu_indexing(self.K, self.node_filter_reverse)
        self.Up = self.u[self.node_filter]

        self.Fp = self.f[self.node_filter_reverse]
        
        
        
        ###############################
        print('solving linear static system')
        self.f_ = self.Fp - np.matmul(self.Kpp, self.Up)
        self.Uu = np.linalg.solve(self.Kpu, self.f_)
        self.u[self.node_filter_reverse] = self.Uu
        self.f = np.matmul(self.K,self.u)
        

        if self.point9[0] == '_':
            self.num = 4
        elif self.point9[0] != '_' and self.point10[0] == '_':
            self.num = 5
        elif self.point10[0] != '_':
            self.num = 6
            

        self.print_(self.f, self.node_filter, self.num)

      
        ###############################
        print('plotting displacement')
        print('max u_x=', max(self.u[0::2]))
        print('max u_y=', max(self.u[1::2]))
        print('max u_z=', max(self.u[2::2]))

        self.x = self.nodes[:,0]
        self.y = self.nodes[:,1]
        self.z = self.nodes[:,2]

        
        
    def plot(self, force_dir):
        self.results_plot(self.u, self.nodes, self.f[self.node_filter], force_dir, self.force_point, self.force,
                          self.mesh_nx, self.mesh_ny, self.mesh_nz, self.point5, self.point6, self.point7, self.point8,
                          self.point9, self.point10)
        
 
        
    def mesh_ratio(self, point1, point2, point3, point4, point5, point6, point7, point8, point9, point10, force_point, mesh_nx, 
                   mesh_ny, mesh_nz, mesh_lx, mesh_ly, mesh_lz):
        mesh_point = [point1, point2, point3, point4, point5, point6, point7, point8, point9, point10, force_point]
        mesh_point_x = []
        mesh_point_y = []
        mesh_point_z = []
        
        for val in mesh_point:
            if val[0] == "_":
                pass
            else:
                mesh_point_x.append(val[0])          ### boundary condition x coordinate
                mesh_point_y.append(val[1])
                mesh_point_z.append(val[2])
        mesh_point_x.append(0)
        mesh_point_y.append(0)
        mesh_point_z.append(0)       
        mesh_point_x.append(mesh_lx)
        mesh_point_y.append(mesh_ly)
        mesh_point_z.append(mesh_lz)
        x = set(mesh_point_x)
        y = set(mesh_point_y)
        z = set(mesh_point_z)
        mesh_point_x = sorted(list(x))
        mesh_point_y = sorted(list(y))
        mesh_point_z = sorted(list(z))
        
        mesh_point_x_ratio = (np.array(mesh_point_x)/mesh_point_x[-1])*mesh_nx
        mesh_point_y_ratio = (np.array(mesh_point_y)/mesh_point_y[-1])*mesh_ny
        mesh_point_z_ratio = (np.array(mesh_point_z)/mesh_point_z[-1])*mesh_nz
        
        return mesh_point_x, mesh_point_y, mesh_point_z, mesh_point_x_ratio, mesh_point_y_ratio, mesh_point_z_ratio

    
    def mesh_creation(self, mesh_point_x, mesh_point_y, mesh_point_z, mesh_point_x_ratio, 
                      mesh_point_y_ratio, mesh_point_z_ratio):
        
        x_interval = mesh_point_x_ratio[1:]-mesh_point_x_ratio[:-1]  ### mesh_interval
        y_interval = mesh_point_y_ratio[1:]-mesh_point_y_ratio[:-1]
        z_interval = mesh_point_z_ratio[1:]-mesh_point_z_ratio[:-1]
        x_interval = [round(x) for x in x_interval]
        x_sum = sum(x_interval[:-1])
        x_interval[-1] = int(mesh_point_x_ratio[-1]-x_sum)
        y_interval = [round(y) for y in y_interval]
        z_interval = [round(z) for z in z_interval]
        z_sum = sum(z_interval[:-1])
        z_interval[-1] = int(mesh_point_z_ratio[-1]-z_sum)
    
        
        mesh_point_x = mesh_point_x[1:] #remove 0
        mesh_point_y = mesh_point_y[1:] 
        mesh_point_z = mesh_point_z[1:]

        
        nodes = []
        li_x = []
        li_y = []
        li_z = []
        i = 0
        j = 0
        k = 0
        for x,y in zip(mesh_point_x, x_interval):
            if len(mesh_point_x) == 1:
                li_x.append(np.linspace(0.0,x,y,endpoint=True))
            else:
                if i == 0:
                    li_x.append(np.linspace(0.0,x,y,endpoint=False))
                    
                elif x == mesh_point_x[-1]:
                    li_x.append(np.linspace(mesh_point_x[i-1],x,y,endpoint=True))
        
                    
                else:
                    li_x.append(np.linspace(mesh_point_x[i-1],x,y,endpoint=False))
            i += 1
            
                
        for l,m in zip(mesh_point_y, y_interval):
            if len(mesh_point_y) == 1:
                li_y.append(np.linspace(0.0,l,m,endpoint=True))
                
            else:    
                if j == 0:
                    li_y.append(np.linspace(0.0,l,m,endpoint=False))
                    
                elif l == mesh_point_y[-1]:
                    li_y.append(np.linspace(mesh_point_y[j-1],l,m,endpoint=True))
    
                    
                else:
                    li_y.append(np.linspace(mesh_point_y[j-1],l,m,endpoint=False))
            j += 1
                
    
            
        for t,w in zip(mesh_point_z, z_interval):
            if len(mesh_point_z) == 1:
                li_z.append(np.linspace(0.0,t,w,endpoint=True))
            else:
                if k == 0:
                    li_z.append(np.linspace(0.0,t,w,endpoint=False))
                    
                elif t == mesh_point_z[-1]:
                    li_z.append(np.linspace(mesh_point_z[k-1],t,w,endpoint=True))
        
                else:
                    li_z.append(np.linspace(mesh_point_z[k-1],t,w,endpoint=False))
            k += 1
             
    
        interval_x = np.concatenate(li_x, axis = None)
        interval_y = np.concatenate(li_y, axis = None)
        interval_z = np.concatenate(li_z, axis = None)

    
        
        for y in interval_y:
            for z in interval_z:
                for x in interval_x:
                    nodes.append([x,y,z])
        nodes = np.array(nodes)
        return nodes


    def mesh_plot(self, nodes, point1, point2, point3, point4, point5, point6, point7, point8, point9, point10, force_point):
        x = np.unique(nodes[:,0:1])
        y = np.unique(nodes[:,1:2])
        z = np.unique(nodes[:,2:3])
        X, Y, Z = np.meshgrid(x, y, z)
        grid = pv.StructuredGrid(X, Y, Z)
        plotter = pv.Plotter()
        plotter.add_mesh(grid, show_edges=True, cmap="viridis")
        plotter.add_point_labels(
            point1,  # Point coordinates
            [f'{point1}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"label_1")
        plotter.add_point_labels(
            point2,  # Point coordinates
            [f'{point2}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"label_2")
        plotter.add_point_labels(
            point3,  # Point coordinates
            [f'{point3}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"label_3")
        plotter.add_point_labels(
            point4,  # Point coordinates
            [f'{point4}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"label_4")
        plotter.add_point_labels(
            point5,  # Point coordinates
            [f'{point5}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"label_5")
        plotter.add_point_labels(
            point6,  # Point coordinates
            [f'{point6}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"label_6")
        plotter.add_point_labels(
            point7,  # Point coordinates
            [f'{point7}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"label_7")
        plotter.add_point_labels(
            point8,  # Point coordinates
            [f'{point8}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"label_8")
    

        if point9[0] != '_':
            plotter.add_point_labels(
                point9,  # Point coordinates
                [f'{point9}'],  # Scalar values as labels
                point_size=10,
                font_size=12,
                text_color="black",
                name=f"label_9")
        if point10[0] != '_':
            plotter.add_point_labels(
                point10,  # Point coordinates
                [f'{point10}'],  # Scalar values as labels
                point_size=10,
                font_size=12,
                text_color="black",
                name=f"label_10")
        plotter.add_point_labels(
            force_point,  # Point coordinates
            [f'{force_point}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"force_point")
        plotter.show()
    
    
    def create_connectivity(self):
        conn = []
        for k in range(self.mesh_ey):
            for j in range(self.mesh_ez):
            	for i in range(self.mesh_ex):
            		n0 = i + j*self.mesh_nx +k*self.mesh_nx*self.mesh_nz
            		conn.append([n0, n0 + self.mesh_nx, n0 + self.mesh_nx + self.mesh_nx*self.mesh_nz, 
                           n0 + self.mesh_nx*self.mesh_nz, n0 + 1, n0 + 1 + self.mesh_nx, n0 + 1 + self.mesh_nx + self.mesh_nx*self.mesh_nz, 
                           n0 + 1 + self.mesh_nx*self.mesh_nz])
        return conn
    
    
    def c_matrix(self, E, v):
        print('material model - plane strain')
        C = E/(1.0+v)/(1.0-2.0*v) * np.array([[1.0-v,        v,        v,        0.0,        0.0,        0.0],
        								      [    v,    1.0-v,        v,        0.0,        0.0,        0.0],
        								      [    v,        v,    1.0-v,        0.0,        0.0,        0.0],
                                              [  0.0,      0.0,      0.0,    (1.0-2.0*v)/2,    0.0,        0.0],
                                              [  0.0,      0.0,      0.0,        0.0,    (1.0-2.0*v)/2,    0.0],
                                              [  0.0,      0.0,      0.0,        0.0,        0.0, (1.0-2.0*v)/2]])
        return C
    
    

    
    def shape(self, xi):
        x,y,z = tuple(xi)
        N = [(1.0-x)*(1.0-y)*(1.0+z), (1.0-x)*(1.0-y)*(1.0-z), (1.0-x)*(1.0+y)*(1.0-z), (1.0-x)*(1.0+y)*(1.0+z),
          (1.0+x)*(1.0-y)*(1.0+z), (1.0+x)*(1.0-y)*(1.0-z), (1.0+x)*(1.0+y)*(1.0-z), (1.0+x)*(1.0+y)*(1.0+z)]
        return 0.125*np.array(N)
    
    def gradshape(self, xi):
    	s,t,z = tuple(xi)
    	dN = [[-(1.0-t)*(1.0+z),  -(1.0-t)*(1.0-z), -(1.0+t)*(1.0-z), -(1.0+t)*(1.0+z), (1.0-t)*(1.0+z), (1.0-t)*(1.0-z), (1.0+t)*(1.0-z), (1.0+t)*(1.0+z)],
    		  [-(1.0-s)*(1.0+z), -(1.0-s)*(1.0-z), (1.0-s)*(1.0-z), (1.0-s)*(1.0+z), -(1.0+s)*(1.0+z), -(1.0+s)*(1.0-z), (1.0+s)*(1.0-z), (1.0+s)*(1.0+z)],
              [(1.0-s)*(1.0-t), -(1.0-s)*(1.0-t), -(1.0-s)*(1.0+t), (1.0-s)*(1.0+t), (1.0+s)*(1.0-t), -(1.0+s)*(1.0-t), -(1.0+s)*(1.0+t), (1.0+s)*(1.0+t)]]
    	return 0.125*np.array(dN)
    
    
    
    def global_stiffness_matrix(self):
        print('create global stiffness matrix')
        K = np.zeros((3*self.num_nodes, 3*self.num_nodes))
        q4 = [[x/math.sqrt(3.0),y/math.sqrt(3.0),z/math.sqrt(3.0)] for z in [-1.0,1.0] for y in [-1.0,1.0] for x in [-1.0,1.0]]
        B = np.zeros((6,24))
        for c in self.conn:
            xIe = self.nodes[c,:]
            Ke = np.zeros((24,24))
            for q in q4:
                dN = self.gradshape(q)
                J  = np.dot(dN, xIe)
                J_inv = np.linalg.inv(J)
                dN = np.dot(J_inv, dN)
           
                B[0,0::3] = dN[0,:]
                B[1,1::3] = dN[1,:]
                B[2,2::3] = dN[2,:]
                B[3,0::3] = dN[1,:]
                B[3,1::3] = dN[0,:]
                B[4,0::3] = dN[2,:]
                B[4,2::3] = dN[0,:]
                B[5,1::3] = dN[2,:]
                B[5,2::3] = dN[1,:]
                Ke += np.dot(np.dot(B.T,self.C),B) * np.linalg.det(J)
                
            for i,I in enumerate(c):
                for j,J in enumerate(c):
                    K[3*I,3*J]     += Ke[3*i,3*j]
                    K[3*I,3*J+1]   += Ke[3*i,3*j+1]
                    K[3*I+1,3*J]   += Ke[3*i+1,3*j]
                    K[3*I+1,3*J+1]   += Ke[3*i+1,3*j+1]
                    K[3*I,3*J+2]   += Ke[3*i,3*j+2]
                    K[3*I+2,3*J]   += Ke[3*i+2,3*j]
                    K[3*I+1,3*J+2]   += Ke[3*i+1,3*j+2]
                    K[3*I+2,3*J+1]   += Ke[3*i+2,3*j+1]
                    K[3*I+2,3*J+2]   += Ke[3*i+2,3*j+2]
        return K
    
    
    def apply_boundary_conditions(self, nodes, num_nodes, point5, point6, point7, point8, point9, point10, 
                                  force_point, u, f, force):

        point_order = []
        point_dic = {}
        node_filter = []
        for i in range(num_nodes):
            if nodes[i,0] == point5[0] and nodes[i,1] == point5[1]  and nodes[i,2]== point5[2]:
                u[3*i] = 0.0     ### X
                u[3*i+1] = 0.0   ### y
                u[3*i+2] = 0.0   ### z
                node_filter.append(3*i)
                node_filter.append(3*i+1)
                node_filter.append(3*i+2)
                point_dic.update({'point5':[3*i,3*i+1,3*i+2]})
                point_order.append(0)
                
            if nodes[i,0] == point6[0] and nodes[i,1] == point6[1] and nodes[i,2] == point6[2]:
                u[3*i] = 0.0     ### X
                u[3*i+1] = 0.0   ### y
                u[3*i+2] = 0.0   ### z
                node_filter.append(3*i)
                node_filter.append(3*i+1)
                node_filter.append(3*i+2)
                point_dic.update({'point6':[3*i,3*i+1,3*i+2]})
                point_order.append(1)
            
            if nodes[i,0] == point7[0] and nodes[i,1] == point7[1] and nodes[i,2] == point7[2]:
                u[3*i] = 0.0     ### i,2X
                u[3*i+1] = 0.0   ### y
                u[3*i+2] = 0.0   ### z
                node_filter.append(3*i)
                node_filter.append(3*i+1)
                node_filter.append(3*i+2)
                point_dic.update({'point7':[3*i,3*i+1,3*i+2]})
                point_order.append(2)

            if nodes[i,0] == point8[0] and nodes[i,1] == point8[1] and nodes[i,2] == point8[2]:
                u[3*i] = 0.0     ### X
                u[3*i+1] = 0.0   ### y
                u[3*i+2] = 0.0   ### z
                node_filter.append(3*i)
                node_filter.append(3*i+1)
                node_filter.append(3*i+2)
                point_dic.update({'point8':[3*i,3*i+1,3*i+2]})
                point_order.append(3)
            
            if nodes[i,0] == point9[0] and nodes[i,1] == point9[1] and nodes[i,2] == point9[2]:
                u[3*i] = 0.0     ### X
                u[3*i+1] = 0.0   ### y
                u[3*i+2] = 0.0   ### z
                node_filter.append(3*i)
                node_filter.append(3*i+1)
                node_filter.append(3*i+2)
                point_dic.update({'point9':[3*i,3*i+1,3*i+2]})
                point_order.append(4)
            
            if nodes[i,0] == point10[0] and nodes[i,1] == point10[1] and nodes[i,2] == point10[2]:
                u[3*i] = 0.0     ### X
                u[3*i+1] = 0.0   ### y
                u[3*i+2] = 0.0   ### z
                node_filter.append(3*i)
                node_filter.append(3*i+1)
                node_filter.append(3*i+2)
                point_dic.update({'point10':[3*i,3*i+1,3*i+2]})
                point_order.append(5)
                
                
                
                
            if nodes[i,0] == force_point[0] and nodes[i,1] == force_point[1] and nodes[i,2] == force_point[2]:
                f[3*i] = force[0]     ### X
                f[3*i+1] = force[1]   ### y
                f[3*i+2] = force[2]   ### z

           

        node_filter_reverse = [x for x in range(3*num_nodes) if x not in node_filter]
        return node_filter, node_filter_reverse, point_dic, point_order

    def Kpp_indexing(self, K,node_filter_reverse,node_filter):
        Kpp = K[node_filter_reverse,:]
        Kpp = Kpp[:,node_filter]
        return Kpp
    
    def Kpu_indexing(self, K,node_filter_reverse):
        Kpu = K[node_filter_reverse,:]
        Kpu = Kpu[:,node_filter_reverse]
        return Kpu
    
    def print_(self, f, node_filter, num):
        if num >= 4:
            print("force R1x", f[self.point_dic["point5"][0]])
            print("force R1y", f[self.point_dic["point5"][1]])
            print("force R1z", f[self.point_dic["point5"][2]])
            print("force R2x", f[self.point_dic["point6"][0]])
            print("force R2y", f[self.point_dic["point6"][1]])
            print("force R2z", f[self.point_dic["point6"][2]] )
            print("force R3x", f[self.point_dic["point7"][0]])
            print("force R3y", f[self.point_dic["point7"][1]])
            print("force R3z", f[self.point_dic["point7"][2]])
            print("force R4x", f[self.point_dic["point8"][0]])
            print("force R4y", f[self.point_dic["point8"][1]])
            print("force R4z", f[self.point_dic["point8"][2]])
        if num == 4:
            print("R1x+R2x+R3x+R4x", f[self.point_dic["point5"][0]]+f[self.point_dic["point6"][0]]+f[self.point_dic["point7"][0]]+f[self.point_dic["point8"][0]])
            print("R1y+R2y+R3y+R4y", f[self.point_dic["point5"][1]]+f[self.point_dic["point6"][1]]+f[self.point_dic["point7"][1]]+f[self.point_dic["point8"][1]])
            print("R1z+R2z+R3z+R4z", f[self.point_dic["point5"][2]]+f[self.point_dic["point6"][2]]+f[self.point_dic["point7"][2]]+f[self.point_dic["point8"][2]])
        if num >= 5:
            print("force R5x", f[self.point_dic["point9"][0]])
            print("force R5y", f[self.point_dic["point9"][1]])
            print("force R5z", f[self.point_dic["point9"][2]])
        if num == 5:
            print("R1x+R2x+R3x+R4x+R5x", f[self.point_dic["point5"][0]]+f[self.point_dic["point6"][0]]+f[self.point_dic["point7"][0]]+f[self.point_dic["point8"][0]]+f[self.point_dic["point9"][0]])
            print("R1y+R2y+R3y+R4y+R5y", f[self.point_dic["point5"][1]]+f[self.point_dic["point6"][1]]+f[self.point_dic["point7"][1]]+f[self.point_dic["point8"][1]]+f[self.point_dic["point9"][1]])
            print("R1z+R2z+R3z+R4z+R5z", f[self.point_dic["point5"][2]]+f[self.point_dic["point6"][2]]+f[self.point_dic["point7"][2]]+f[self.point_dic["point8"][2]]+f[self.point_dic["point9"][2]])
            
        if num == 6:
            print("force R6x", f[self.point_dic["point10"][0]])
            print("force R6y", f[self.point_dic["point10"][1]])
            print("force R6z", f[self.point_dic["point10"][2]])
            print("R1x+R2x+R3x+R4x+R5x+R6x", f[self.point_dic["point5"][0]]+f[self.point_dic["point6"][0]]+f[self.point_dic["point7"][0]]+f[self.point_dic["point8"][0]]+f[self.point_dic["point9"][0]]+f[self.point_dic["point10"][0]])
            print("R1y+R2y+R3y+R4y+R5y+R6y", f[self.point_dic["point5"][1]]+f[self.point_dic["point6"][1]]+f[self.point_dic["point7"][1]]+f[self.point_dic["point8"][1]]+f[self.point_dic["point9"][1]]+f[self.point_dic["point10"][1]])
            print("R1z+R2z+R3z+R4z+R5z+R6z", f[self.point_dic["point5"][2]]+f[self.point_dic["point6"][2]]+f[self.point_dic["point7"][2]]+f[self.point_dic["point8"][2]]+f[self.point_dic["point9"][2]]+f[self.point_dic["point10"][2]])
    
    
    def results_plot(self, u, nodes, f_filtered, force_dir,force_point, force, mesh_nx, mesh_ny, mesh_nz, point5,
                     point6, point7, point8, point9, point10):
        
        if force_dir == 'X':
            num = 0
            arrow_dir = [1,0,0]
        elif force_dir == 'Y':
            num = 1
            arrow_dir = [0,1,0]
        else:
            num = 2
            arrow_dir = [0,0,1]
        
        u = u[1::3] #y축 변위만 추출
        u_matrix = u.reshape(mesh_ny, mesh_nz, mesh_nx)
        u_matrix = np.transpose(u_matrix, (1, 2, 0))
        u_1d = u_matrix.flatten()
    
        x = np.unique(nodes[:,0:1])
        y = np.unique(nodes[:,1:2])
        z = np.unique(nodes[:,2:3])
        X, Y, Z = np.meshgrid(x, y, z)
        grid = pv.StructuredGrid(X, Y, Z)
        grid.point_data["Displacement"] =u_1d
        plotter = pv.Plotter()
        plotter.add_mesh(grid, show_edges=True, cmap="viridis")
        points = [point5, point6, point7, point8, point9, point10]
        points = [points[x] for x in self.point_order]
        
        for i, value in enumerate(f_filtered[num::3]):
            plotter.add_point_labels(
                points[i],  # Point coordinates
                [round(value,3)],  # Scalar values as labels
                point_size=10,
                font_size=15,
                text_color="black",
                name=f'label_{i}')
            
            # Add arrows to the corners
            points[i][num] = points[i][num]
            arrows = pv.Arrow(points[i], arrow_dir, scale=((self.mesh_lx+self.mesh_ly+self.mesh_lz)/12))
            plotter.add_mesh(arrows, color='red')
        
        # Add force values and arrows
        plotter.add_point_labels(
            force_point,  # Point coordinates
            [f'{force_point}'],  # Scalar values as labels
            point_size=10,
            font_size=12,
            text_color="black",
            name=f"label_8")
        force_point[1] = force_point[1]
        arrows = pv.Arrow(force_point, force, scale=((self.mesh_lx+self.mesh_ly+self.mesh_lz)/12))
        plotter.add_mesh(arrows, color='red')
        
        plotter.show()
    




    





