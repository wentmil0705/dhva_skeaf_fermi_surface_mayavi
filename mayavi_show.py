# from traits.etsconfig.api import ETSConfig
# ETSConfig.toolkit = 'qt4'
# from math import gamma, ceil
import os
# from turtle import bgcolor
# from unittest.case import doModuleCleanups
os.environ['ETS_TOOLKIT'] = 'qt'
os.environ['QT_API'] = 'pyqt'
from pyface.qt import QtGui, QtCore
from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor
import numpy as np
import os
import re
from scipy.interpolate import interpn
import matplotlib.pyplot as plt
import pandas as pd
from skimage.measure import marching_cubes as marching_cubes
from scipy.spatial import cKDTree
from scipy.spatial import Voronoi
import hdbscan


class Visualization(HasTraits):
    CONVAU2ANG = 0.529177209
    scene = Instance(MlabSceneModel, (), editor=SceneEditor())

    def init_plot(self):
        
        self.scene.scene_editor.background = (1, 1, 1)
        self.scene.mlab.plot3d([0,1], [0,1], [0,1], color=(1, 1, 1))
        self.orbit_data = []
        self.bz = np.array([])
        pass


    def update_plot(self, bcell, ebands, efermi, para_dict):

        self.scene.mlab.clf()
        self.scene.scene_editor.background = para_dict['background_color']
        all_l = []
        tree = []
        self.orbit_data = para_dict['orbit_data']
        self.bz = np.asarray(para_dict['bz_number'])
        if para_dict['mode'] == 0:
            all_l, tree = self.plot_simple_fermi_surface(bcell, ebands, efermi, para_dict)
        else:
            all_l, tree = self.plot_seperate_fermi_surface(bcell, ebands, efermi, para_dict)
        if para_dict['show_slice']:
            self.plot_slice(bcell, para_dict, all_l, tree)
        if para_dict['trait']:
            self.plot_trait(bcell, para_dict, tree)
        if para_dict['axes']:
            self.scene.mlab.orientation_axes()
        if para_dict['magnetic_field']:
            if (para_dict['trait'] == True) & (para_dict['is_mag_h'] == True):
                mag_h = para_dict['mag_h']
            else:
                mag_h = para_dict['section_v']
                if para_dict['ang_flag']:
                    mag_h = mag_h * Visualization.CONVAU2ANG / (2 * np.pi)
                else:
                    mag_h = mag_h / (2 * np.pi)
            mag_h = mag_h / np.sqrt(np.sum(mag_h ** 2))
            self.scene.mlab.quiver3d(0.06, 0.05, 0.06, mag_h[0], mag_h[1], mag_h[2], color=(0, 0, 0), scale_factor=.022, mode='arrow') 
        pass


    def get_brillouin_zone_3d(self, cell):

        cell = np.asarray(cell, dtype=float)
        assert cell.shape == (3, 3)

        px, py, pz = np.tensordot(cell, np.mgrid[-1:2, -1:2, -1:2], axes=[0, 0])
        points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

        vor = Voronoi(points)

        bz_facets = []
        bz_ridges = []
        bz_vertices = []

        for pid, rid in zip(vor.ridge_points, vor.ridge_vertices):
            if(pid[0] == 13 or pid[1] == 13):
                bz_ridges.append(vor.vertices[np.r_[rid, [rid[0]]]])
                bz_facets.append(vor.vertices[rid])
                bz_vertices += rid

        bz_vertices = list(set(bz_vertices))

        return vor.vertices[bz_vertices], bz_ridges, bz_facets
    

    def get_primitive_cell_3d(self, cell):

        cell = np.asarray(cell, dtype=float)
        assert cell.shape == (3, 3)

        dx, dy, dz = np.mgrid[0:2, 0:2, 0:2]
        dxyz = np.c_[dx.ravel(), dy.ravel(), dz.ravel()]
        px, py, pz = np.tensordot(cell, [dx, dy, dz], axes=[0, 0])
        points = np.c_[px.ravel(), py.ravel(), pz.ravel()]

        lines = []
        faces = None

        for ii in range(len(points)):
            for jj in range(ii):
                if np.abs(dxyz[ii] - dxyz[jj]).sum() == 1:
                    lines.append(np.vstack([points[ii], points[jj]]))

        return points, lines, faces

    def plot_primitive(self, bcell, ebands, efermi, para_dict):

        p, l, f = self.get_primitive_cell_3d(bcell)
        bz = self.bz.astype(int)
        b1, b2, b3 = np.linalg.norm(bcell, axis=1)
        b3d = np.tile(ebands, tuple(bz))
        nx, ny, nz = ebands.shape
        b3d = np.pad(b3d, (0,1), mode='wrap')
        verts, faces_in_fs, normals, values = marching_cubes(b3d,
                                                        level=efermi,
                                                        spacing=(
                                                            b1/nx, b2/ny, b3/nz)
                                                        )
        
        bz_list = []
        for i in range(bz[0]):
            for j in range(bz[1]):
                for k in range(bz[2]):
                    bz_list.append([i, j, k])
        bz_list = bz_list[1:]

        all_l = [l]
        for bz_array in bz_list:
            bz_array = np.asarray(bz_array)
            new_l = []
            for part in l:
                new_part = np.dot(part, np.linalg.inv(bcell)) + bz_array
                new_part = np.dot(new_part, bcell)
                new_l.append(new_part)
            all_l.append(new_l)

        if para_dict['line']:
            for p_l in all_l:
                for xx in p_l:
                    self.scene.mlab.plot3d(xx[:, 0], xx[:, 1], xx[:, 2],
                                tube_radius=para_dict['line_width'],
                                color=para_dict['line_color'],
                                transparent=True,
                                opacity=para_dict['line_opacity'])
        

        verts_cart = np.dot(
                verts / np.array([b1, b2, b3]),
                bcell
            )
        
        return verts_cart, faces_in_fs, all_l
    

    def plot_bz(self, bcell, ebands, efermi, para_dict):

        p, l, f = self.get_brillouin_zone_3d(bcell)
        b1, b2, b3 = np.linalg.norm(bcell, axis=1)
        # bz = np.array(list(map(lambda x: int(x), para_dict['bz_number'])))
        bz = self.bz.astype(int)
        bz_max = int(np.max(bz)) + 1
        px, py, pz = np.tensordot(
            bcell,
            np.mgrid[-1:bz_max, -1:bz_max, -1:bz_max],
            axes=[0, 0]
        )
        points = np.c_[px.ravel(), py.ravel(), pz.ravel()]
        tree = cKDTree(points)
        gamma_region_id = tree.query([0, 0, 0])[1]

        b3d_2uc = np.tile(ebands, (bz_max, bz_max, bz_max))
        nx, ny, nz = b3d_2uc.shape
        verts, faces, normals, values = marching_cubes(b3d_2uc,
                                        level=efermi,
                                        spacing=(
                                            bz_max*b1/nx, bz_max*b2/ny, bz_max*b3/nz)
                                        )
        region_list = [gamma_region_id]

        bz_list = []
        for i in range(bz[0]):
            for j in range(bz[1]):
                for k in range(bz[2]):
                    bz_list.append([i, j, k])
        bz_list = bz_list[1:]

        dm = int(np.cbrt(np.max(tree.indices) + 1))
        test = np.reshape(np.arange(int(dm**3)), (dm, dm, dm))
        all_l = [l]
        for bz_array in bz_list:
            new_l = []
            for part in l:
                new_part = np.dot(part, np.linalg.inv(bcell)) + bz_array
                new_part = np.dot(new_part, bcell)
                new_l.append(new_part)
            all_l.append(new_l)
            new_id_loc = np.asarray(np.where(test == gamma_region_id))[:, 0] + bz_array
            region_list.append(test[new_id_loc[0]][new_id_loc[1]][new_id_loc[2]])

        if para_dict['line']:
            for p_l in all_l:
                for xx in p_l:
                    self.scene.mlab.plot3d(xx[:, 0], xx[:, 1], xx[:, 2],
                                tube_radius=para_dict['line_width'] / 1.5,
                                color=para_dict['line_color'],
                                transparent=True,
                                opacity=para_dict['line_opacity'])

        verts_cart = np.dot(
                                    verts / np.array([b1, b2, b3]) - np.ones(3),
                                    bcell
                                )
                                
        verts_region_id = tree.query(verts_cart)[1]
        verts_in_bz = np.array(list(map(lambda x: x in region_list, verts_region_id)))        
        faces_in_fs = faces[np.all(verts_in_bz[faces], axis=1)]
        vertices_old_id = np.unique(faces_in_fs)
        vertices_new_id = range(vertices_old_id.size)
        old_new_map = dict(np.c_[vertices_old_id, vertices_new_id])

        verts_cart = verts_cart[vertices_old_id]
        faces_in_fs = [[old_new_map[v] for v in f] for f in faces_in_fs]

        
        return verts_cart, faces_in_fs, all_l, tree
 
        
    def plot_simple_fermi_surface(self, bcell, ebands, efermi, para_dict):

        tree = []
        if para_dict['bz_mode'] == 2:            
            if (para_dict['trait']) & (len(self.orbit_data) != 0):
                self.move_bz_orbit_data(bcell, para_dict)
            verts_cart, faces_in_fs, all_l, tree = self.plot_bz(bcell, ebands, efermi, para_dict)
            if para_dict['inner_outer']:
                verts_cart_inner, faces_in_fs_inner, all_l, tree = self.plot_bz(bcell, ebands, efermi+0.0005, para_dict)

        else:
            if (para_dict['trait']) & (len(self.orbit_data) != 0):
                self.move_pm_orbit_data(bcell, para_dict)
            verts_cart, faces_in_fs, all_l = self.plot_primitive(bcell, ebands, efermi, para_dict)
            if para_dict['inner_outer']:
                verts_cart_inner, faces_in_fs_inner, all_l = self.plot_primitive(bcell, ebands, efermi+0.0005, para_dict)

        self.scene.mlab.triangular_mesh(verts_cart[:,0], verts_cart[:,1], verts_cart[:,2],
                                faces_in_fs,
                                color=para_dict['outer_color'],
                                transparent=True,
                                opacity=para_dict['fs_opacity'],
                                scalars=np.linalg.norm(verts_cart, axis=1)
                                )
        
        if para_dict['inner_outer']:
            self.scene.mlab.triangular_mesh(verts_cart_inner[:,0], verts_cart_inner[:,1], verts_cart_inner[:,2],
                                    faces_in_fs_inner,
                                    color=para_dict['inner_color'],
                                    transparent=True,
                                    opacity=para_dict['fs_opacity'],
                                    scalars=np.linalg.norm(verts_cart_inner, axis=1)
            )
        
        return all_l, tree


    def plot_seperate_fermi_surface(self, bcell, ebands, efermi, para_dict):

        tree = []

        if para_dict['bz_mode'] == 2:
            if (para_dict['trait']) & (len(self.orbit_data) != 0):
                self.move_bz_orbit_data(bcell, para_dict)
            verts_cart, faces_in_fs, all_l, tree = self.plot_bz(bcell, ebands, efermi, para_dict)
        else:
            if (para_dict['trait']) & (len(self.orbit_data) != 0):
                self.move_pm_orbit_data(bcell, para_dict)
            verts_cart, faces_in_fs, all_l = self.plot_primitive(bcell, ebands, efermi, para_dict)

        knn = hdbscan.HDBSCAN(prediction_data=True)

        verts_cart_sorts = knn.fit_predict(verts_cart)
        sort_count = np.max(verts_cart_sorts) + 1
        verts_cart_copy = verts_cart.copy()
        colors = plt.get_cmap('rainbow', sort_count)
        self.plot_trait(bcell, para_dict, tree)
        select_zones = []

        for i in para_dict['select_row']:
            orbit_data = self.orbit_data[i]
            x_new = np.mean(orbit_data, axis=0)
            predict_zone, s = hdbscan.approximate_predict(knn, [x_new])
            predict_zone = predict_zone[0]
            select_zones.append(predict_zone)
            # self.scene.mlab.text3d(x_new[0], x_new[1], x_new[2], str(predict_zone), color=(1, 0, 0), scale=.005)
        ## grey模式下的整个布里渊区
        self.scene.mlab.triangular_mesh(verts_cart[:,0], verts_cart[:,1], verts_cart[:,2],
                                         faces_in_fs,
                                         colormap='Greys',
                                         opacity=.3,
                                         transparent=True,
                                         scalars=np.linalg.norm(verts_cart, axis=1)
                                         )
        
        # print(select_zones)

        ## 特殊模块选取绘图                               
        verts_cart_copy = verts_cart.copy()
        colors = []
        for i in select_zones:        
            color_map = plt.get_cmap('rainbow', sort_count)
            colors.append(tuple(list(color_map(i))[:-1]))
        df_copy = pd.DataFrame(np.c_[verts_cart_copy, verts_cart_sorts], columns=['x', 'y', 'z', 's'])
        df = df_copy.copy()
        # for i in range(sort_count):
        #     df = pd.DataFrame(np.c_[verts_cart_copy, verts_cart_sorts], columns=['x', 'y', 'z', 's'])
        #     df = df[df['s'] == i]
        #     mid = np.mean(df[['x', 'y', 'z']].values, axis=0)
        #     self.scene.mlab.text3d(mid[0], mid[1], mid[2], str(i), color=(0, 0, 0), scale=.005)
        for count, i in enumerate(select_zones):
            plot_sort = i
            df = pd.DataFrame(np.c_[verts_cart_copy, verts_cart_sorts], columns=['x', 'y', 'z', 's'])
            df.loc[df['s'] != plot_sort, 'x'] = np.nan
            df.loc[df['s'] != plot_sort, 'y'] = np.nan
            df.loc[df['s'] != plot_sort, 'z'] = np.nan
            verts_cart = df[['x', 'y', 'z']].values

            self.scene.mlab.triangular_mesh(verts_cart[:,0], verts_cart[:,1], verts_cart[:,2],
                                                    faces_in_fs,
                                                    color=colors[count],
                                                    opacity=para_dict['fs_opacity'],
                                                    transparent=True,
                                                    scalars=np.linalg.norm(verts_cart, axis=1)
                                                    )
        
        return all_l, tree


    def plot_slice(self, bcell, para_dict, all_l, tree):
        b1, b2, b3 = np.linalg.norm(bcell, axis=1)
        if (para_dict['trait'] == True) & (para_dict['is_mag_h'] == True):
            mag_h = para_dict['mag_h']
        else:
            mag_h = para_dict['section_v']
            mag_h = np.dot(mag_h, bcell)
            # if para_dict['ang_flag']:
            #     mag_h = mag_h * Visualization.CONVAU2ANG / (2 * np.pi)
            # else:
            #     mag_h = mag_h / (2 * np.pi)
        mag_h = mag_h / np.sqrt(np.sum(mag_h ** 2))
        self.scene.mlab.quiver3d(0.06, 0.05, 0.06, mag_h[0], mag_h[1], mag_h[2], color=(0, 0, 0), scale_factor=.022, mode='arrow')

        ## 寻找截面用到的函数
        def inner_line(P, P1, P2):
            P_xmax = P1[0] if P1[0] > P2[0] else P2[0]
            P_xmin = P1[0] if P1[0] < P2[0] else P2[0]
            P_ymax = P1[1] if P1[1] > P2[1] else P2[1]
            P_ymin = P1[1] if P1[1] < P2[1] else P2[1]
            P_zmax = P1[2] if P1[2] > P2[2] else P2[2]
            P_zmin = P1[2] if P1[2] < P2[2] else P2[2]
            if (P[0] > P_xmax) | (P[0] < P_xmin):
                return False
            if (P[1] > P_ymax) | (P[1] < P_ymin):
                return False
            if (P[2] > P_zmax) | (P[2] < P_zmin):
                return False
            return True

        def Find_intersection(P1, P2, a, b, c, d):
            m = P1[0] - P2[0]
            n = P1[1] - P2[1]
            p = P1[2] - P2[2]
            t = (-a * P1[0] - b * P1[1] - c * P1[2] - d) / (a * m + b * n + c * p)

            x = m * t + P1[0]
            y = n * t + P1[1]
            z = p * t + P1[2]
            X = np.array([x, y, z])

            return X
        
        slices = []
        if para_dict['trait'] is False:
            slices = [para_dict['section_v']]
        else:
            for i in para_dict['select_row']:
                orbit_data = self.orbit_data[i]
                slices.append(np.mean(orbit_data, axis=0))

        for a in slices:
            a_new = a

            h_slice = np.r_[mag_h, -np.dot(a_new, mag_h)]


            slice_all_points = []
            all_triangles = []
            for p_l in all_l:
                slice_points = []
                plus_num = len(slice_all_points)
                for bzplotx in p_l:
                    for i in range(len(bzplotx)-1):
                        vertices_1 = bzplotx[i]
                        vertices_2 = bzplotx[i+1]
                        # mlab.text3d(vertices_mid[0], vertices_mid[1], vertices_mid[2], str(i), color=(0, 1, 0), scale=.005)
                        cross_point = Find_intersection(vertices_1, vertices_2, h_slice[0], h_slice[1], h_slice[2], h_slice[3])           
                        if (np.isnan(cross_point).any() == 0) & (inner_line(cross_point, vertices_1, vertices_2)):
                            slice_points.append(cross_point)

                slice_points = np.array(slice_points)
                slice_points = np.around(slice_points, 6)
                slice_points = np.unique(slice_points, axis=0)
                if len(slice_points) != 0:
                    slice_points = slice_points[np.lexsort((slice_points[:, 0], slice_points[:, 1], slice_points[:, 2]))]
                    # for i, p in enumerate(slice_points):
                    #     mlab.text3d(p[0], p[1], p[2], str(i), color=(0, 1, 0), scale=.005)
                    triangles = []
                    target_index_list = []
                    source_index_list = []
                    src_id = 0
                    tar_id = 1
                    all_mid = np.mean(slice_points, axis=0)
                    avg_dist = np.mean(np.sqrt(np.sum((slice_points - all_mid)**2, axis=1)))
                    while len(target_index_list) < len(slice_points)-1:
                        source_index_list.append(src_id)
                        near_dist = 100.
                        for i in [x for x in range(len(slice_points)) if (x not in target_index_list) & (x != src_id) & (x not in source_index_list)]:
                            dist = np.sqrt(np.sum((slice_points[i] - slice_points[src_id]) ** 2))
                            mid = 1/2 * (slice_points[i] + slice_points[src_id])
                            mm_dist = np.sqrt(np.sum((mid - all_mid) ** 2)) 
                            if (dist < near_dist) & (mm_dist > avg_dist/10):
                                near_dist = dist
                                tar_id = i
                        target_index_list.append(tar_id)
                        triangles.append((len(slice_points)+plus_num, src_id+plus_num, tar_id+plus_num))
                        src_id = tar_id
                    triangles.append((len(slice_points)+plus_num, src_id+plus_num, 0+plus_num))
                    slice_points = np.r_[slice_points, [np.mean(slice_points, axis=0)]]
                    all_triangles = all_triangles + triangles
                    slice_all_points = slice_all_points + slice_points.tolist()
            
            if len(slice_all_points) != 0:
                slice_all_points = np.asarray(slice_all_points)
                slice_x = slice_all_points[:, 0]
                slice_y = slice_all_points[:, 1]
                slice_z = slice_all_points[:, 2]

                # print(all_triangles)
                # for i, p in enumerate(slice_all_points):
                #     self.scene.mlab.text3d(p[0], p[1], p[2], str(i), scale=.002, color=(0, 0, 0))
                self.scene.mlab.triangular_mesh(slice_x, slice_y, slice_z, all_triangles, color=para_dict['slice_color'], transparent=True, opacity=para_dict['slice_opacity'])
                
    
    def plot_trait(self, bcell, para_dict, tree):

        if para_dict['bz_mode'] == 2:
            gamma_region_id = tree.query([0, 0, 0])[1]

            for row in para_dict['select_row']:
                orbit_data = self.orbit_data[row]

                if len(orbit_data) != 0:
                    if para_dict['full_trait'] == False:
                        orbit_region_id = tree.query(orbit_data)[1]
                        orbit_in_bz = (orbit_region_id == gamma_region_id)
                        orbit_data_in_bz = orbit_data[orbit_in_bz]
                        orbit_data = orbit_data_in_bz
                        self.scene.mlab.points3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2],
                                scale_factor=para_dict['trait_width']*4,
                                color=para_dict['trait_color'],
                                transparent=True,
                                opacity=para_dict['trait_opacity'])
                    else:
                        self.scene.mlab.plot3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2],
                                tube_radius=para_dict['trait_width'],
                                color=para_dict['trait_color'],
                                transparent=True,
                                opacity=para_dict['trait_opacity'])

                else:
                    print("The orbit of frequency "+str(para_dict['freqs'][row])+" kT not in first bz")
        
        else:
            for row in para_dict['select_row']:
                orbit_data = self.orbit_data[row]
                if len(orbit_data) != 0:
                    if para_dict['full_trait'] == False:
                        orig_orbit_data = np.dot(orbit_data, np.linalg.inv(bcell))
                        more_in_l = np.c_[orig_orbit_data[:, 0] < 0, orig_orbit_data[:, 1] < 0, orig_orbit_data[:, 2] < 0]
                        more_in_l = np.sum(more_in_l, axis=1)
                        orig_in_pm = orig_orbit_data[more_in_l < 1]
                        orig_orbit_data = orig_in_pm
                        more_in_l = np.c_[orig_orbit_data[:, 0] > 1, orig_orbit_data[:, 1] > 1, orig_orbit_data[:, 2] > 1]
                        more_in_l = np.sum(more_in_l, axis=1)
                        orig_in_pm = orig_orbit_data[more_in_l < 1]
                        orbit_data = np.dot(orig_in_pm, bcell)

                        self.scene.mlab.points3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2],
                                scale_factor=para_dict['trait_width']*4,
                                color=para_dict['trait_color'],
                                transparent=True,
                                opacity=para_dict['trait_opacity'])
                    else:
                        self.scene.mlab.plot3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2],
                                tube_radius=para_dict['trait_width'],
                                color=para_dict['trait_color'],
                                transparent=True,
                                opacity=para_dict['trait_opacity'])
                else:
                    print("The orbit of frequency "+str(para_dict['freqs'][row])+" kT not in first bz")
    
    def move_pm_orbit_data(self, bcell, para_dict):

        new_orbit_data = []
        bz_max_list = []
        for orbit_data in self.orbit_data:
            orbit_data = np.asarray(orbit_data)
            if para_dict['ang_flag']:
                orbit_data = orbit_data * Visualization.CONVAU2ANG
            orbit_data = orbit_data / (2 * np.pi)

            def move(a, a_ref):
                moveArray = []
                for i in a:
                    if i > 0:
                        moveArray.append(-int(i/a_ref))
                    else:
                        moveArray.append(-int(i/a_ref)+1)
                return np.array(moveArray)

            orig_data = np.dot(orbit_data, np.linalg.inv(bcell))
            if para_dict['full_fs_full_trait']:
                moveArray = move(np.min(orig_data, axis=0), 1)
                orig_data = orig_data + moveArray
                bz_max_list.append(np.ceil(np.max(orig_data, axis=0)))
                new_data = np.dot(orig_data, bcell)
                new_orbit_data.append(new_data)
            else:
                moveArray = move(np.mean(orig_data, axis=0), 1)
                orig_data = orig_data + moveArray
                new_data = np.dot(orig_data, bcell)
                new_orbit_data.append(new_data)
        self.orbit_data = new_orbit_data
        if para_dict['full_fs_full_trait']:
            bz_max_list = np.asarray(bz_max_list)
            self.bz = np.max(bz_max_list, axis=0)            

    
    def move_bz_orbit_data(self, bcell, para_dict):
        new_orbit_data = []
        bz_max_list = []
        bz_max = np.max(self.bz)
        px, py, pz = np.tensordot(
            bcell,
            np.mgrid[-1:bz_max+1, -1:bz_max+1, -1:bz_max+1],
            axes=[0, 0]
        )
        points = np.c_[px.ravel(), py.ravel(), pz.ravel()]
        tree = cKDTree(points)
        for orbit_data in self.orbit_data:
            orbit_data = np.asarray(orbit_data)
            if para_dict['ang_flag']:
                orbit_data = orbit_data * Visualization.CONVAU2ANG
            orbit_data = orbit_data / (2 * np.pi)

            orig_data = np.dot(orbit_data, np.linalg.inv(bcell))
            if para_dict['full_fs_full_trait']:
                dm = int(np.cbrt(np.max(tree.indices) + 1))
                test = np.reshape(np.arange(int(dm**3)), (dm, dm, dm))
                gamma_region_id = tree.query([0, 0, 0])[1]
                zone_id = np.unique(tree.query(orbit_data)[1])
                zone_id_loc = []
                for i in zone_id:
                    zone_id_loc.append(np.asarray(np.where(test == i))[:, 0])
                zone_id_dif = -np.asarray(zone_id_loc) + np.asarray(np.where(test == gamma_region_id))[:, 0]
                moveArray = np.max(zone_id_dif, axis=0)
                new_data = np.dot(np.dot(orbit_data, np.linalg.inv(bcell)) + moveArray, bcell)
                new_orbit_data.append(new_data)
                zone_id = np.unique(tree.query(new_data)[1])
                zone_id_loc = []
                for i in zone_id:
                    zone_id_loc.append(np.asarray(np.where(test == i))[:, 0])
                zone_id_dif = np.asarray(zone_id_loc) - np.asarray(np.where(test == gamma_region_id))[:, 0] + 1
                bz_max_list.append(np.ceil(np.max(zone_id_dif, axis=0))) 
            else:
                orbit_id = tree.query(np.dot(np.mean(orig_data, axis=0), bcell))[1]
                dm = int(np.cbrt(np.max(tree.indices) + 1))
                test = np.reshape(np.arange(int(dm**3)), (dm, dm, dm))
                gamma_region_id = tree.query([0, 0, 0])[1]
                moveArray = np.asarray(np.where(test == gamma_region_id))[:, 0] - np.asarray(np.where(test == orbit_id))[:, 0]
                orig_data = orig_data + moveArray
                new_data = np.dot(orig_data, bcell)
                new_orbit_data.append(new_data)
        if para_dict['full_fs_full_trait']:
            bz_max_list = np.asarray(bz_max_list)
            self.bz = np.max(bz_max_list, axis=0)
        self.orbit_data = new_orbit_data 


    # the layout of the dialog screated
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True  # We need this to resize with the parent widget
                )



class MayaviQWidget(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self.visualization = Visualization()

        self.ui = self.visualization.edit_traits(parent=self, kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)
        self.visualization.init_plot()
    
    def bxsfMayavi(self, bcell, ebands, efermi, para_dict):
        self.visualization.update_plot(bcell, ebands, efermi, para_dict)

