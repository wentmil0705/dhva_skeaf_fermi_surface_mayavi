
__author__ = "Guanzhang Liu"
__maintainer__ = "Guanzhang Liu"
__email__ = "MG21340029@smail.nju.edu.cn"
__date__ = "October 23, 2022"


'''
mayavi_show.py为mayavi_widget的设计文件文件
'''

import os

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
from mayavi.tools.pipeline import glyph, surface, tube
from mayavi import mlab
import gc
from share import *
from scipy.spatial import ConvexHull

class Visualization(HasTraits):
    # 创建全局绘图界面（scene仅为一个界面，若多创建则会产生多个界面）
    # CONVAU2ANG = 0.529177209
    scene = Instance(MlabSceneModel, (), editor=SceneEditor())

    '''
    初始化界面
    '''
    def init_plot(self):
        
        self.scene.scene_editor.background = (1, 1, 1)
        # 初始化图形，可通过设计实现，目前仅是简单的图案
        self.initial = self.scene.mlab.plot3d([0,1], [0,1], [0,1], color=(1, 1, 1))
        self.AXES = self.scene.mlab.orientation_axes()
        self.AXES.stop()
        self.MAG = None
        self.orbit_data_list_full_bz = []
        self.orbit_data_list_not_full_bz = []
        self.orbit_data_list_full_pm = []
        self.orbit_data_list_not_full_pm = []
        self.orbit_data_list_full_bz_tube = []
        self.orbit_data_list_full_pm_tube = []
        self.fermi_bz_inner = None
        self.fermi_bz_outer = None
        self.fermi_pm_inner = None
        self.fermi_pm_outer = None
        self.seperate_bg_bz = None
        self.seperate_bg_pm = None
        self.bz_edge = []
        self.pm_edge = []
        self.bz_edge_tube = []
        self.pm_edge_tube = []
        self.seperate_orbit_area_bz = []
        self.seperate_orbit_area_pm = []
        self.all_l = []
        self.tree = None
        self.bz_size_full = []
        self.pm_size_full = []
        self.bz = []
        self.tree_pm = None
        self.all_l_pm = []
        self.tree_bz = None
        self.all_l_bz = []
        self.knn_bz = None
        self.knn_pm = None
        self.verts_cart_sorts_bz = []
        self.verts_cart_sorts_pm = []
        self.sort_count_bz = None
        self.sort_count_pm = None
        self.verts_cart_bz = []
        self.verts_cart_pm = []
        self.faces_in_fs_bz = []
        self.faces_in_fs_pm = []
        self.orbit_data_list_bz = []
        self.orbit_data_list_pm = []
        self.select_mag = None
        self.gv_vector_mlab_bz = None
        self.gv_mag_mlab_bz = None
        self.gv_magText_mlab_bz = None
        self.gv_vector_mlab_pm = None
        self.gv_mag_mlab_pm = None
        self.gv_magText_mlab_pm = None
        self.region_list = []
        self.orbit_slice_list = None
        
        pass

    # B_direction窗口关闭触发取消磁场设定方向显示
    def update_select_mag_end(self):
        if self.select_mag is not None:
            self.select_mag.stop()
            del self.select_mag
        gc.collect()
    
    # slice窗口sync后触发的同步截面位置
    def update_synSlice(self, para_dict, mid_points):
        if self.orbit_slice_list is not None:
            for p in self.orbit_slice_list:
                p.stop()
            del self.orbit_slice_list
            gc.collect()
        self.orbit_slice_list = []
        try:
            self.sub_plot_slice(mid_points, para_dict)
        except:
            pass 

    # 计算群速度
    def update_calc_gv(self, bcell, ebands):
        if self.gv_vector_mlab_bz is not None:
                self.gv_vector_mlab_bz.stop()
                self.gv_mag_mlab_bz.stop()
                self.gv_magText_mlab_bz.stop()
                self.gv_vector_mlab_pm.stop()
                self.gv_mag_mlab_pm.stop()
                self.gv_magText_mlab_pm.stop()
                del self.gv_vector_mlab_bz, self.gv_mag_mlab_bz, self.gv_magText_mlab_bz, self.gv_vector_mlab_pm, self.gv_mag_mlab_pm, self.gv_magText_mlab_pm
                gc.collect()
        b1, b2, b3 = np.linalg.norm(bcell, axis=1)
        new_verts_cart_bz = self.verts_cart_bz / [b1, b2, b3]
        new_verts_cart_bz = np.dot(self.verts_cart_bz, np.linalg.inv(bcell)) + np.ones(3)
        new_verts_cart_pm = self.verts_cart_pm / [b1, b2, b3]
        new_verts_cart_pm = np.dot(self.verts_cart_pm, np.linalg.inv(bcell))
        c_1 = np.fromiter(map(lambda x: x-(x // 1) if ((x > 1) | (x < 0)) else x, new_verts_cart_bz[:, 0]), dtype=float)
        c_2 = np.fromiter(map(lambda x: x-(x // 1) if ((x > 1) | (x < 0)) else x, new_verts_cart_bz[:, 1]), dtype=float)
        c_3 = np.fromiter(map(lambda x: x-(x // 1) if ((x > 1) | (x < 0)) else x, new_verts_cart_bz[:, 2]), dtype=float)
        new_verts_cart_bz = np.c_[c_1, c_2, c_3]
        c_1 = np.fromiter(map(lambda x: x-(x // 1) if ((x > 1) | (x < 0)) else x, new_verts_cart_pm[:, 0]), dtype=float)
        c_2 = np.fromiter(map(lambda x: x-(x // 1) if ((x > 1) | (x < 0)) else x, new_verts_cart_pm[:, 1]), dtype=float)
        c_3 = np.fromiter(map(lambda x: x-(x // 1) if ((x > 1) | (x < 0)) else x, new_verts_cart_pm[:, 2]), dtype=float)
        new_verts_cart_pm = np.c_[c_1, c_2, c_3]

        ## 梯度求值
        def sub_calc_gv(values, length, new_edge, new_verts_cart):
            points = (np.linspace(0, new_edge[0], length[0]), np.linspace(0, new_edge[1], length[1]), np.linspace(0, new_edge[2], length[2]))
            nx = np.linspace(0, new_edge[0], length[0])
            ny = np.linspace(0, new_edge[1], length[1]) 
            nz = np.linspace(0, new_edge[2], length[2])
            X,Y,Z = np.meshgrid(nx,ny,nz)
            XYZ = np.concatenate((X[:,:,:,None], Y[:,:,:,None], Z[:,:,:,None]),axis=-1)
            delta = new_edge / length / 10
            x_plus_delta = new_verts_cart + [delta[0], 0, 0]
            y_plus_delta = new_verts_cart + [0, delta[1], 0]
            z_plus_delta = new_verts_cart + [0, 0, delta[2]]
            x_plus_delta_e = interpn(points, values, x_plus_delta, method='linear')
            y_plus_delta_e = interpn(points, values, y_plus_delta, method='linear')
            z_plus_delta_e = interpn(points, values, z_plus_delta, method='linear')
            orig_e = interpn(points, values, new_verts_cart, method='linear')
            grad_v_x = (x_plus_delta_e - orig_e) / delta[0]
            grad_v_y = (y_plus_delta_e - orig_e) / delta[1]
            grad_v_z = (z_plus_delta_e - orig_e) / delta[2]
            grad_v_list = np.c_[grad_v_x, grad_v_y, grad_v_z]

            ## 转换单位
            import math
            METER_ANGSTROM = 10**(-10) #m /A
            lattice = np.linalg.inv(bcell.T).T     
            grad_v_list = grad_v_list/(2*math.pi)
            grad_v_list = np.multiply(grad_v_list, np.array([np.linalg.norm(lattice[:,0])*METER_ANGSTROM * len(np.unique(XYZ[:,0])),
                                                                    np.linalg.norm(lattice[:,1])*METER_ANGSTROM * len(np.unique(XYZ[:,1])),
                                                                    np.linalg.norm(lattice[:,2])*METER_ANGSTROM * len(np.unique(XYZ[:,2]))])
                                )

            # grad_v_list_cart = np.array([np.matmul(lattice, gradient) for gradient in grad_v_list])

            ## 群速度
            HBAR_EV = 6.582119 *10**(-16) #eV*s
            gv_x = grad_v_list[:,0]/HBAR_EV
            gv_y = grad_v_list[:,1]/HBAR_EV
            gv_z = grad_v_list[:,2]/HBAR_EV

            gv_vector = [gv_x,gv_y,gv_z]

            gv_mag = np.c_[gv_vector[0], gv_vector[1], gv_vector[2]]
            gv_mag = np.sum(gv_mag**2, axis=1)

            return gv_vector, gv_mag
        
        values = np.pad(ebands, (0,1), mode='wrap')
        length = np.array(values.shape, dtype=int)
        new_edge = 1 + 1 / (length-1)
        gv_vector_pm, gv_mag_pm = sub_calc_gv(values, length, new_edge, new_verts_cart_pm)
        gv_vector_bz, gv_mag_bz = sub_calc_gv(values, length, new_edge, new_verts_cart_bz)


        self.gv_vector_mlab_bz = self.scene.mlab.quiver3d(self.verts_cart_bz[:, 0], self.verts_cart_bz[:, 1], self.verts_cart_bz[:, 2], gv_vector_bz[0], gv_vector_bz[1], gv_vector_bz[2])
        self.gv_mag_mlab_bz = self.scene.mlab.triangular_mesh(self.verts_cart_bz[:,0], self.verts_cart_bz[:,1], self.verts_cart_bz[:,2],
                                        self.faces_in_fs_bz,
                                        scalars=np.sqrt(gv_mag_bz)
                                    )
        self.gv_magText_mlab_bz = self.scene.mlab.title(str(np.mean(np.sqrt(gv_mag_bz))  / 13.61)+' m/s')
        self.gv_vector_mlab_pm = self.scene.mlab.quiver3d(self.verts_cart_pm[:, 0], self.verts_cart_pm[:, 1], self.verts_cart_pm[:, 2], gv_vector_pm[0], gv_vector_pm[1], gv_vector_pm[2])
        self.gv_mag_mlab_pm = self.scene.mlab.triangular_mesh(self.verts_cart_pm[:,0], self.verts_cart_pm[:,1], self.verts_cart_pm[:,2],
                                        self.faces_in_fs_pm,
                                        scalars=np.sqrt(gv_mag_pm)
                                    )
        self.gv_magText_mlab_pm = self.scene.mlab.title(str(np.mean(np.sqrt(gv_mag_pm))  / 13.61)+' m/s')
        self.gv_vector_mlab_bz.stop()
        self.gv_mag_mlab_bz.stop()
        self.gv_magText_mlab_bz.stop()
        self.gv_vector_mlab_pm.stop()
        self.gv_mag_mlab_pm.stop()
        self.gv_magText_mlab_pm.stop() 
    
    # 更新群速度显示
    def update_gv(self, para_dict):
        self.fermi_bz_inner.stop()
        self.fermi_bz_outer.stop()
        self.fermi_pm_inner.stop()
        self.fermi_pm_outer.stop()
        self.shut_every_trait()
        self.seperate_bg_bz.stop()
        self.seperate_bg_pm.stop()
        self.gv_vector_mlab_bz.stop()
        self.gv_mag_mlab_bz.stop()
        self.gv_magText_mlab_bz.stop()
        self.gv_vector_mlab_pm.stop()
        self.gv_mag_mlab_pm.stop()
        self.gv_magText_mlab_pm.stop()           
        self.bz_edge[-1].start()
        if para_dict['bz_mode'] == 2:
            if para_dict['vector']:
                self.gv_vector_mlab_bz.start()
            if para_dict['mag']:
                self.gv_mag_mlab_bz.start()
            if para_dict['avg_mag']:
                self.gv_magText_mlab_bz.start()
        else:
            if para_dict['vector']:
                self.gv_vector_mlab_pm.start()
            if para_dict['mag']:
                self.gv_mag_mlab_pm.start()
            if para_dict['avg_mag']:
                self.gv_magText_mlab_pm.start()
        self.update_edge(para_dict)

    # 更新B_direction窗口所选磁场方向
    def update_select_mag(self, bcell, mag):
        if self.select_mag is not None:
            self.select_mag.stop()
            del self.select_mag
        pos = np.sum(bcell, axis=0)
        mag = mag / np.sqrt(np.sum(mag**2))
        self.select_mag = self.scene.mlab.quiver3d(pos[0], pos[1], pos[2], mag[0], mag[1], mag[2], color=(0, 0, 1), scale_factor=.022, mode='arrow')
        self.select_mag.start()

    # 更新背景颜色
    def update_bg(self, para_dict):
        self.scene.scene_editor.background = para_dict['background_color']

    # 更新图像显示位置
    def update_camera(self):
        return self.scene.mlab.view()[0:3]

    # 更新切片位置
    def update_show_slice(self, bcell, para_dict):
        self.plot_slice(bcell, para_dict)
    
    # 更新轨道计算
    def update_calc_trait(self, bcell, para_dict):
        self.shut_every_trait()
        self.calculate_trait(bcell, para_dict)
        gc.collect()
    
    # 更新轴显示
    def update_axes(self, para_dict):
        if para_dict['axes']:
            self.AXES.start()
        else:
            self.AXES.stop()
    
    # 更新磁场方向显示（目前与section-v无关）
    def update_mag(self, para_dict):
        ## 绘制磁场方向
        if self.MAG is not None:
            self.MAG.stop()
            del self.MAG
        # if para_dict['is_mag_h'] == True:
        #     mag_h = para_dict['mag_h']
        # else:
        #     mag_h = para_dict['mag_h']
        #     if para_dict['ang_flag']:
        #         mag_h = mag_h * Visualization.CONVAU2ANG / (2 * np.pi)
        #     else:
        #         mag_h = mag_h / (2 * np.pi)
        mag_h = para_dict['mag_h']
        mag_h = mag_h / np.sqrt(np.sum(mag_h ** 2))
        self.MAG = self.scene.mlab.quiver3d(0.06, 0.05, 0.06, mag_h[0], mag_h[1], mag_h[2], color=(0, 0, 0), scale_factor=.022, mode='arrow')
        if para_dict['magnetic_field']:
            self.MAG.start()
        else:
            self.MAG.stop()
    
    # 更新轨道显示
    def update_circular_trait(self, orbit_list, para_dict, mode, orbit_list2=None):
        if self.orbit_data_list_full_bz != []:
            for row in para_dict['select_row']:
                orbit_list[row].start()
                if mode == 'p':
                    tmp = glyph(orbit_list[row], 
                                scale_factor=para_dict['trait_width']*4,
                                color=para_dict['trait_color'],
                                transparent=True,
                                opacity=para_dict['trait_opacity'])
                else:
                    orbit_list2[row].start()
                    tmp_tube = tube(orbit_list2[row], tube_radius=para_dict['trait_width'])
                    tmp_tube.start()
                    tmp = surface(tmp_tube,
                                  color=para_dict['trait_color'],
                                  transparent=True,
                                  opacity=para_dict['trait_opacity'])
                    tmp_tube.stop()
                    orbit_list2[row].stop()
                orbit_list[row].stop()
                # del orbit_list[row]                
                # gc.collect()
                orbit_list[row] = tmp
                orbit_list[row].start()
    def shut_every_trait(self):
        len_freqs = len(self.orbit_data_list_full_bz)
        if self.orbit_data_list_full_bz != []:
            for i in range(len_freqs):
                self.orbit_data_list_full_bz[i].stop()
                self.orbit_data_list_not_full_bz[i].stop()
                self.orbit_data_list_full_pm[i].stop()
                self.orbit_data_list_not_full_pm[i].stop()
                self.seperate_orbit_area_bz[i].stop()
                self.seperate_orbit_area_pm[i].stop()
    def update_trait(self, para_dict):
        self.shut_every_trait()
        if para_dict['bz_mode'] == 2:
            if para_dict['full_trait']:
                self.update_circular_trait(self.orbit_data_list_full_bz, para_dict, 'l', self.orbit_data_list_full_bz_tube)
            else:               
                self.update_circular_trait(self.orbit_data_list_not_full_bz, para_dict, 'p')
        else:
            if para_dict['full_trait']:
                self.update_circular_trait(self.orbit_data_list_full_pm, para_dict, 'l', self.orbit_data_list_full_pm_tube)
            else:
                self.update_circular_trait(self.orbit_data_list_not_full_pm, para_dict, 'p')
        if para_dict['mode'] == 1:
            self.update_seperate(para_dict)
        self.plot_slice(para_dict)
        
    # 更新边界显示
    def update_edge(self, para_dict):
        for p in self.bz_edge:
            p.stop()
        for p in self.pm_edge:
            p.stop()
        if para_dict['bz_mode'] == 2:
            if para_dict['line']:
                del self.bz_edge
                self.bz_edge = []
                for p in range(len(self.bz_edge_tube)):
                    self.bz_edge_tube[p].start()
                    tmp = tube(self.bz_edge_tube[p], tube_radius=para_dict['line_width'])
                    tmp_new = surface(tmp,                                
                                color=para_dict['line_color'],
                                transparent=True,
                                opacity=para_dict['line_opacity'])
                    self.bz_edge_tube[p].stop()
                    self.bz_edge.append(tmp_new)
                    tmp.stop()
                    del tmp                    
                    self.bz_edge[-1].start()
        else:
            if para_dict['line']:
                del self.pm_edge
                self.pm_edge = []
                for p in range(len(self.pm_edge_tube)):
                    self.pm_edge_tube[p].start()
                    tmp = tube(self.pm_edge_tube[p], tube_radius=para_dict['line_width'])
                    tmp_new = surface(tmp,                                
                                color=para_dict['line_color'],
                                transparent=True,
                                opacity=para_dict['line_opacity']) 
                    self.pm_edge_tube[p].stop()
                    self.pm_edge.append(tmp_new)
                    tmp.stop()
                    del tmp                    
                    self.pm_edge[-1].start()
        
    # 更新费米面显示
    def update_fs(self, para_dict):
        self.fermi_bz_inner.stop()
        self.fermi_bz_outer.stop()
        self.fermi_pm_inner.stop()
        self.fermi_pm_outer.stop()
        self.seperate_bg_pm.stop()
        self.seperate_bg_bz.stop()
        if self.gv_mag_mlab_bz is not None:
            self.gv_vector_mlab_bz.stop()
            self.gv_mag_mlab_bz.stop()
            self.gv_magText_mlab_bz.stop()
            self.gv_vector_mlab_pm.stop()
            self.gv_mag_mlab_pm.stop()
            self.gv_magText_mlab_pm.stop() 
        self.update_trait(para_dict)
        self.update_edge(para_dict)
        if para_dict['bz_mode'] == 2:
            if para_dict['mode'] == 0:
                self.fermi_bz_outer.start()
                tmp = surface(self.fermi_bz_outer,
                            color=para_dict['outer_color'],
                            transparent=True,
                            opacity=para_dict['fs_opacity'],
                            )
                self.fermi_bz_outer.stop()
                del self.fermi_bz_outer
                self.fermi_bz_outer = tmp
                self.fermi_bz_outer.start()
                if para_dict['inner_outer']:
                    self.fermi_bz_inner.start()
                    tmp = surface(self.fermi_bz_inner,
                                color=para_dict['inner_color'],
                                transparent=True,
                                opacity=para_dict['fs_opacity'],
                                )
                    self.fermi_bz_inner.stop()
                    del self.fermi_bz_inner                    
                    self.fermi_bz_inner = tmp
                    self.fermi_bz_inner.start()
            else:
                self.seperate_bg_bz.start()
        else:
            if para_dict['mode'] == 0:
                self.fermi_pm_outer.start()
                tmp = surface(self.fermi_pm_outer,
                            color=para_dict['outer_color'],
                            transparent=True,
                            opacity=para_dict['fs_opacity'],
                            )
                self.fermi_pm_outer.stop()
                del self.fermi_pm_outer                
                self.fermi_pm_outer = tmp
                self.fermi_pm_outer.start()
                if para_dict['inner_outer']:
                    self.fermi_pm_inner.start()
                    tmp = surface(self.fermi_pm_inner,
                                color=para_dict['inner_color'],
                                transparent=True,
                                opacity=para_dict['fs_opacity'],
                                )
                    self.fermi_pm_inner.stop()
                    del self.fermi_pm_inner                    
                    self.fermi_pm_inner = tmp
                    self.fermi_pm_inner.start()
            else:
                self.seperate_bg_pm.start()
        gc.collect()
        
    # 更新seperate mode显示
    def update_seperate(self, para_dict):
        if para_dict['bz_mode'] == 2:
            self.seperate_bg_bz.start()
        else:
            self.seperate_bg_pm.start()
        for i in para_dict['select_row']:
            if para_dict['bz_mode'] == 2:
                self.seperate_orbit_area_bz[i].start()
            else:
                self.seperate_orbit_area_pm[i].start()
    
    # 更新费米面计算
    def update_calc_fs(self, bcell, ebands, efermi, para_dict):
        # 暂时的解决方案
        if self.bz_edge is not None:
            for p in self.bz_edge:
                p.stop()
        if self.pm_edge is not None:
            for p in self.pm_edge:
                p.stop()
        mlab.clf(figure=self.scene.mayavi_scene)
        self.scene.scene_editor.background = para_dict['background_color']
        del self.all_l
        del self.tree
        del self.bz        
        self.all_l = []  # 记录所需显示的bz区域
        self.tree = []   # 记录bz树编号
        if para_dict['full_fs_full_trait']:
            self.bz = self.bz_size_full if (para_dict['bz_mode'] == 2) else self.pm_size_full
        else:
            self.bz = np.asarray(para_dict['bz_number'])
        self.calc_simple_fermi_surface(bcell, ebands, efermi, para_dict)
        self.all_l = self.all_l_bz
        self.tree = self.tree_bz
        gc.collect()

    '''
    更新参数清除图片后的绘图
    '''
    def update_plot(self, bcell, ebands, efermi, para_dict):

        ## 清清除目前的图像
        self.scene.mlab.clf()
        self.scene.scene_editor.background = para_dict['background_color']
        self.all_l = []  # 记录所需显示的bz区域
        self.tree = []   # 记录bz树编号
        # self.orbit_data = para_dict['orbit_data']
        self.bz = np.asarray(para_dict['bz_number'])

        ## mode绘图
        if para_dict['mode'] == 0:
            self.calc_simple_fermi_surface(bcell, ebands, efermi, para_dict)
            self.update_fs(para_dict)
            self.all_l = self.all_l_bz
            self.tree = self.tree_bz
        else:
            self.calc_seperate_fermi_surface(bcell, ebands, efermi, para_dict)
            self.update_fs(para_dict)
            self.all_l = self.all_l_pm
            self.tree = self.tree_pm
        
        ## slice绘图
        if para_dict['show_slice']:
            self.update_show_slice(bcell, para_dict)
        
        ## 轨道绘图
        if para_dict['trait']:
            self.update_trait(bcell, para_dict)
        
        ## 绘制轴显示
        if para_dict['axes']:
            self.scene.mlab.orientation_axes()
        
        ## 绘制磁场方向
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

    # 布里渊区
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
    
    # 初基元胞
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

    # 绘制初基元胞费米面
    def plot_primitive(self, bcell, ebands, efermi, para_dict):

        del self.pm_edge_tube, self.pm_edge
        self.pm_edge_tube = []
        self.pm_edge = []
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

        for p_l in all_l:
            for xx in p_l:
                tmp = self.scene.mlab.plot3d(xx[:, 0], xx[:, 1], xx[:, 2])
                tmp.stop()
                self.pm_edge_tube.append(self.scene.mlab.points3d(xx[:, 0], xx[:, 1], xx[:, 2], opacity=0))
                self.pm_edge_tube[-1].mlab_source.dataset.lines = tmp.mlab_source.dataset.lines
                tmp = tube(self.pm_edge_tube[-1], tube_radius=para_dict['line_width'] / 1.5)
                self.pm_edge.append(surface(tmp,
                            color=para_dict['line_color'],
                            transparent=True,
                            opacity=para_dict['line_opacity']))
                tmp.stop()
                self.pm_edge_tube[-1].stop()
                self.pm_edge[-1].stop()
        

        verts_cart = np.dot(
                verts / np.array([b1, b2, b3]),
                bcell
            )
        
        return verts_cart, faces_in_fs, all_l
    
    # 绘制布里渊区费米面
    def plot_bz(self, bcell, ebands, efermi, para_dict):

        del self.bz_edge_tube, self.bz_edge
        self.bz_edge_tube = []
        self.bz_edge = []
        p, l, f = self.get_brillouin_zone_3d(bcell)
        b1, b2, b3 = np.linalg.norm(bcell, axis=1)
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
        self.region_list = [gamma_region_id]

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
            self.region_list.append(test[new_id_loc[0]][new_id_loc[1]][new_id_loc[2]])

        for p_l in all_l:
            for xx in p_l:
                tmp = self.scene.mlab.plot3d(xx[:, 0], xx[:, 1], xx[:, 2])
                tmp.stop()
                self.bz_edge_tube.append(self.scene.mlab.points3d(xx[:, 0], xx[:, 1], xx[:, 2], opacity=0))
                self.bz_edge_tube[-1].mlab_source.dataset.lines = tmp.mlab_source.dataset.lines
                tmp = tube(self.bz_edge_tube[-1], tube_radius=para_dict['line_width'] / 1.5)
                self.bz_edge.append(surface(tmp,
                            color=para_dict['line_color'],
                            transparent=True,
                            opacity=para_dict['line_opacity']))
                tmp.stop()
                self.bz_edge_tube[-1].stop()
                self.bz_edge[-1].stop()

        verts_cart = np.dot(
                                    verts / np.array([b1, b2, b3]) - np.ones(3),
                                    bcell
                                )
                                
        verts_region_id = tree.query(verts_cart)[1]
        verts_in_bz = np.array(list(map(lambda x: x in self.region_list, verts_region_id)))        
        faces_in_fs = faces[np.all(verts_in_bz[faces], axis=1)]
        vertices_old_id = np.unique(faces_in_fs)
        vertices_new_id = range(vertices_old_id.size)
        old_new_map = dict(np.c_[vertices_old_id, vertices_new_id])

        verts_cart = verts_cart[vertices_old_id]
        faces_in_fs = [[old_new_map[v] for v in f] for f in faces_in_fs]

        
        return verts_cart, faces_in_fs, all_l, tree
 
    # 计算简单费米面
    def calc_simple_fermi_surface(self, bcell, ebands, efermi, para_dict):

        del self.tree_pm
        self.tree_pm = []
        del self.verts_cart_bz, self.faces_in_fs_bz, self.all_l_bz, self.tree_bz
        del self.fermi_bz_outer
        self.verts_cart_bz, self.faces_in_fs_bz, self.all_l_bz, self.tree_bz = self.plot_bz(bcell, ebands, efermi, para_dict)        
        self.fermi_bz_outer = self.scene.mlab.triangular_mesh(self.verts_cart_bz[:,0], self.verts_cart_bz[:,1], self.verts_cart_bz[:,2],
                            self.faces_in_fs_bz,
                            color=para_dict['outer_color'],
                            transparent=True,
                            opacity=para_dict['fs_opacity'],
                            scalars=np.linalg.norm(self.verts_cart_bz, axis=1)
                            )
        del self.tree_bz, self.all_l_bz
        del self.fermi_bz_inner
        verts_cart_inner, faces_in_fs_inner, self.all_l_bz, self.tree_bz = self.plot_bz(bcell, ebands, efermi+0.0005, para_dict)
        self.fermi_bz_inner = self.scene.mlab.triangular_mesh(verts_cart_inner[:,0], verts_cart_inner[:,1], verts_cart_inner[:,2],
                            faces_in_fs_inner,
                            color=para_dict['inner_color'],
                            transparent=True,
                            opacity=para_dict['fs_opacity'],
                            scalars=np.linalg.norm(verts_cart_inner, axis=1)
                        )                 
        self.fermi_bz_inner.stop()
        self.fermi_bz_outer.stop()

        del self.verts_cart_pm, self.faces_in_fs_pm, self.all_l_pm
        del self.fermi_pm_outer
        self.verts_cart_pm, self.faces_in_fs_pm, self.all_l_pm = self.plot_primitive(bcell, ebands, efermi, para_dict)        
        self.fermi_pm_outer = self.scene.mlab.triangular_mesh(self.verts_cart_pm[:,0], self.verts_cart_pm[:,1], self.verts_cart_pm[:,2],
                                self.faces_in_fs_pm,
                                color=para_dict['outer_color'],
                                transparent=True,
                                opacity=para_dict['fs_opacity'],
                                scalars=np.linalg.norm(self.verts_cart_pm, axis=1)
                                )
        del self.all_l_pm
        del self.fermi_pm_inner
        verts_cart_inner, faces_in_fs_inner, self.all_l_pm = self.plot_primitive(bcell, ebands, efermi+0.0005, para_dict)
        self.fermi_pm_inner = self.scene.mlab.triangular_mesh(verts_cart_inner[:,0], verts_cart_inner[:,1], verts_cart_inner[:,2],
                                faces_in_fs_inner,
                                color=para_dict['inner_color'],
                                transparent=True,
                                opacity=para_dict['fs_opacity'],
                                scalars=np.linalg.norm(verts_cart_inner, axis=1)
                                )
        self.fermi_pm_outer.stop()
        self.fermi_pm_inner.stop()
    
    # 计算seperate费米面
    def calc_seperate_fermi_surface(self):
        del self.knn_bz, self.knn_pm
        self.knn_bz = hdbscan.HDBSCAN(prediction_data=True)
        self.knn_pm = hdbscan.HDBSCAN(prediction_data=True)

        del self.verts_cart_sorts_bz, self.verts_cart_sorts_pm, self.sort_count_bz, self.sort_count_pm, self.seperate_bg_bz, self.seperate_bg_pm
        self.verts_cart_sorts_bz = self.knn_bz.fit_predict(self.verts_cart_bz)
        self.verts_cart_sorts_pm = self.knn_pm.fit_predict(self.verts_cart_pm)
        self.sort_count_bz = np.max(self.verts_cart_sorts_bz) + 1
        self.sort_count_pm = np.max(self.verts_cart_sorts_pm) + 1
        self.seperate_bg_bz = self.scene.mlab.triangular_mesh(self.verts_cart_bz[:,0], self.verts_cart_bz[:,1], self.verts_cart_bz[:,2],
                                         self.faces_in_fs_bz,
                                         colormap='Greys',
                                         opacity=.3,
                                         transparent=True,
                                         scalars=np.linalg.norm(self.verts_cart_bz, axis=1)
                                         )
        self.seperate_bg_pm = self.scene.mlab.triangular_mesh(self.verts_cart_pm[:,0], self.verts_cart_pm[:,1], self.verts_cart_pm[:,2],
                                         self.faces_in_fs_pm,
                                         colormap='Greys',
                                         opacity=.3,
                                         transparent=True,
                                         scalars=np.linalg.norm(self.verts_cart_pm, axis=1)
                                         )
        self.seperate_bg_bz.stop()
        self.seperate_bg_pm.stop()
    
    def calc_seperate_orbit(self, para_dict):
        select_zones_bz = []
        select_zones_pm = []
        select_rows = range(len(para_dict['freqs']))

        for i in select_rows:
            x_new_bz = np.mean(self.orbit_data_list_bz[i], axis=0)
            x_new_pm = np.mean(self.orbit_data_list_pm[i], axis=0)
            predict_zone_bz, s = hdbscan.approximate_predict(self.knn_bz, [x_new_bz])
            predict_zone_pm, s = hdbscan.approximate_predict(self.knn_pm, [x_new_pm])
            select_zones_bz.append(predict_zone_bz[0])
            select_zones_pm.append(predict_zone_pm[0])
        
        def show_every_trait_area(verts_cart, faces_in_fs, select_zones, sort_count, verts_cart_sorts, area_list):
            verts_cart_copy = verts_cart.copy()
            colors = []            
            color_map = plt.get_cmap('rainbow', sort_count)
            for i in select_zones:        
                colors.append(tuple(list(color_map(i))[:-1]))
            df_copy = pd.DataFrame(np.c_[verts_cart_copy, verts_cart_sorts], columns=['x', 'y', 'z', 's'])
            df = df_copy.copy()
            for count, i in enumerate(select_zones):
                plot_sort = i
                df = pd.DataFrame(np.c_[verts_cart_copy, verts_cart_sorts], columns=['x', 'y', 'z', 's'])
                df.loc[df['s'] != plot_sort, 'x'] = np.nan
                df.loc[df['s'] != plot_sort, 'y'] = np.nan
                df.loc[df['s'] != plot_sort, 'z'] = np.nan
                verts_cart = df[['x', 'y', 'z']].values

                area_list.append(self.scene.mlab.triangular_mesh(verts_cart[:,0], verts_cart[:,1], verts_cart[:,2],
                                                        faces_in_fs,
                                                        color=colors[count],
                                                        opacity=para_dict['fs_opacity'],
                                                        transparent=True,
                                                        scalars=np.linalg.norm(verts_cart, axis=1)
                                                        ))
                area_list[-1].stop()
        
        del self.seperate_orbit_area_bz, self.seperate_orbit_area_pm
        self.seperate_orbit_area_bz = []
        self.seperate_orbit_area_pm = []
        show_every_trait_area(self.verts_cart_bz, self.faces_in_fs_bz, select_zones_bz, self.sort_count_bz, self.verts_cart_sorts_bz, self.seperate_orbit_area_bz)
        show_every_trait_area(self.verts_cart_pm, self.faces_in_fs_pm, select_zones_pm, self.sort_count_pm, self.verts_cart_sorts_pm, self.seperate_orbit_area_pm)
        
    # 绘制所选轨道切片
    def plot_slice(self, para_dict):
        if self.orbit_slice_list is not None:
            for p in self.orbit_slice_list:
                p.stop()
            del self.orbit_slice_list
            gc.collect()
        self.orbit_slice_list = []        
        if para_dict['show_slice'] & para_dict['trait']:
            for row in para_dict['select_row']:
                if para_dict['bz_mode'] == 2:
                    mid_points = np.mean(self.orbit_data_list_bz[row],axis=0)
                else:
                    mid_points = np.mean(self.orbit_data_list_pm[row],axis=0) 
                self.sub_plot_slice(mid_points, para_dict) 

    # 绘制切片
    def sub_plot_slice(self, mid_points, para_dict):
        h_slice = np.r_[para_dict['mag_h'], -np.dot(mid_points, para_dict['mag_h'])]
        if para_dict['bz_mode'] == 2:
            all_l = self.all_l_bz
        else:
            all_l = self.all_l_pm
        for p_l in all_l:
            slice_points = []
            for bzplotx in p_l:
                for i in range(len(bzplotx)-1):
                    vertices_1 = bzplotx[i]
                    vertices_2 = bzplotx[i+1]
                    cross_point = Find_intersection(vertices_1, vertices_2, h_slice[0], h_slice[1], h_slice[2], h_slice[3])           
                    if (np.isnan(cross_point).any() == 0) & (inner_line(cross_point, vertices_1, vertices_2)):
                        slice_points.append(cross_point)
            if len(slice_points) != 0:
                slice_points = np.array(slice_points)
                slice_points = np.around(slice_points, 6)
                slice_points = np.unique(slice_points, axis=0)
                points = np.array(slice_points)
                points_3d = np.dot(points, np.linalg.inv(para_dict['sc'] * np.sqrt(np.sum(para_dict['bcell']**2, axis=1))))    
                points_2d = points_3d[:,0:2]
                hull = ConvexHull(points_2d)
                slice_ps_xy = points_3d[hull.vertices]
                trans_dict = dict(zip(hull.vertices, np.arange(len(hull.vertices))))
                trans_angular = np.c_[np.fromiter(map(lambda x: trans_dict[x], hull.simplices[:,0]), dtype=int), np.fromiter(map(lambda x: trans_dict[x], hull.simplices[:,1]), dtype=int)]
                trans_angular = list(map(lambda x: tuple(list(x)+[len(trans_angular)]), trans_angular))
                slice_ps_xy = np.r_[slice_ps_xy, [np.mean(slice_ps_xy, axis=0)]]
                slice_ps_xy = np.dot(slice_ps_xy, para_dict['sc'] * np.sqrt(np.sum(para_dict['bcell']**2, axis=1))) 
                self.orbit_slice_list.append(self.scene.mlab.triangular_mesh(slice_ps_xy[:,0], slice_ps_xy[:,1], slice_ps_xy[:,2], trans_angular, color=para_dict['slice_color'], transparent=True, opacity=para_dict['slice_opacity']))
    
    # 计算轨道
    def calculate_trait(self, bcell, para_dict):

        del self.orbit_data_list_full_bz, self.orbit_data_list_full_pm, self.orbit_data_list_full_bz_tube, self.orbit_data_list_full_pm_tube
        del self.orbit_data_list_not_full_bz, self.orbit_data_list_not_full_pm
        select_row = range(len(para_dict['freqs']))
        self.orbit_data_list_not_full_bz = []
        self.orbit_data_list_full_bz = []
        self.orbit_data_list_full_bz_tube = []
        self.orbit_data_list_not_full_pm = []
        self.orbit_data_list_full_pm = []
        self.orbit_data_list_full_pm_tube = []
        gamma_region_id = self.tree.query([0, 0, 0])[1]
        
        del self.orbit_data_list_bz
        self.orbit_data_list_bz = self.move_bz_orbit_data(bcell, para_dict)
        for row in select_row:
            orbit_data = self.orbit_data_list_bz[row]

            if len(orbit_data) != 0:              
                tmp = self.scene.mlab.plot3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2])
                tmp.stop()
                self.orbit_data_list_full_bz_tube.append(self.scene.mlab.points3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2], opacity=0))
                self.orbit_data_list_full_bz_tube[-1].mlab_source.dataset.lines = tmp.mlab_source.dataset.lines
                tmp = tube(self.orbit_data_list_full_bz_tube[-1], tube_radius=para_dict['trait_width'])
                self.orbit_data_list_full_bz.append(surface(tmp,
                                color=para_dict['trait_color'],
                                transparent=True,
                                opacity=para_dict['trait_opacity']))
                self.orbit_data_list_full_bz_tube[-1].stop()
                tmp.stop()
                self.orbit_data_list_full_bz[-1].stop

                orbit_region_id = self.tree.query(orbit_data)[1]
                orbit_in_bz = (orbit_region_id == gamma_region_id)
                orbit_data_in_bz = orbit_data[orbit_in_bz]
                orbit_data = orbit_data_in_bz
                self.orbit_data_list_not_full_bz.append(self.scene.mlab.points3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2],
                        scale_factor=para_dict['trait_width']*4,
                        color=para_dict['trait_color'],
                        transparent=True,
                        opacity=para_dict['trait_opacity']))
                self.orbit_data_list_not_full_bz[-1].stop()
            else:
                print("The orbit of frequency "+str(para_dict['freqs'][row])+" kT not in first bz")
        
        del self.orbit_data_list_pm
        self.orbit_data_list_pm = self.move_pm_orbit_data(bcell, para_dict)
        for row in select_row:
            orbit_data = self.orbit_data_list_pm[row]
            if len(orbit_data) != 0:
                tmp = self.scene.mlab.plot3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2])
                tmp.stop()
                self.orbit_data_list_full_pm_tube.append(self.scene.mlab.points3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2], opacity=0))
                self.orbit_data_list_full_pm_tube[-1].mlab_source.dataset.lines = tmp.mlab_source.dataset.lines
                tmp = tube(self.orbit_data_list_full_pm_tube[-1], tube_radius=para_dict['trait_width'])
                self.orbit_data_list_full_pm.append(surface(tmp,
                                color=para_dict['trait_color'],
                                transparent=True,
                                opacity=para_dict['trait_opacity']))
                self.orbit_data_list_full_pm_tube[-1].stop()
                tmp.stop()
                self.orbit_data_list_full_pm[-1].stop()

                orig_orbit_data = np.dot(orbit_data, np.linalg.inv(bcell))
                more_in_l = np.c_[orig_orbit_data[:, 0] < 0, orig_orbit_data[:, 1] < 0, orig_orbit_data[:, 2] < 0]
                more_in_l = np.sum(more_in_l, axis=1)
                orig_in_pm = orig_orbit_data[more_in_l < 1]
                orig_orbit_data = orig_in_pm
                more_in_l = np.c_[orig_orbit_data[:, 0] > 1, orig_orbit_data[:, 1] > 1, orig_orbit_data[:, 2] > 1]
                more_in_l = np.sum(more_in_l, axis=1)
                orig_in_pm = orig_orbit_data[more_in_l < 1]
                orbit_data = np.dot(orig_in_pm, bcell)

                self.orbit_data_list_not_full_pm.append(self.scene.mlab.points3d(orbit_data[:, 0], orbit_data[:, 1], orbit_data[:, 2],
                            scale_factor=para_dict['trait_width']*4,
                            color=para_dict['trait_color'],
                            transparent=True,
                            opacity=para_dict['trait_opacity']))
                self.orbit_data_list_not_full_pm[-1].stop()
            else:
                print("The orbit of frequency "+str(para_dict['freqs'][row])+" kT not in first bz")
                
    # 将轨道移至初基元胞中
    def move_pm_orbit_data(self, bcell, para_dict):

        new_orbit_data = []
        bz_max_list = []
        orbit_data_list = para_dict['orbit_data']
        for orbit_data in orbit_data_list:
            orbit_data = np.asarray(orbit_data)
            if para_dict['ang_flag']:
                orbit_data = orbit_data * CONVAU2ANG
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
        if para_dict['full_fs_full_trait']:
            bz_max_list = np.asarray(bz_max_list)
            self.pm_size_full = np.max(bz_max_list, axis=0)
        return new_orbit_data          

    # 将轨道移至布里渊区中
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
        orbit_data_list = para_dict['orbit_data']
        for orbit_data in orbit_data_list:
            orbit_data = np.asarray(orbit_data)
            if para_dict['ang_flag']:
                orbit_data = orbit_data * CONVAU2ANG
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
                print(tree.query(np.sum(new_data, axis=0))[1])
        if para_dict['full_fs_full_trait']:
            bz_max_list = np.asarray(bz_max_list)
            self.bz_size_full = np.max(bz_max_list, axis=0)
        return new_orbit_data

    # the layout of the dialog screated
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True  # We need this to resize with the parent widget
                )


class MayaviQWidget(QtGui.QWidget):
    camera = QtCore.pyqtSignal(float, float, float)
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

    def calcFsMayavi(self, bcell, ebands, efermi, para_dict):
        self.visualization.update_calc_fs(bcell, ebands, efermi, para_dict)
    
    def axesMayavi(self, para_dict):
        self.visualization.update_axes(para_dict)

    def magMayavi(self, para_dict):
        self.visualization.update_mag(para_dict)
    
    def traitCalcMayavi(self, bcell, para_dict):
        self.visualization.update_calc_trait(bcell, para_dict)
    
    def traitMayavi(self, para_dict):
        self.visualization.update_trait(para_dict)

    def fsMayavi(self, para_dict):
        self.visualization.update_fs(para_dict)
    
    def seperateCalcFsMayavi(self):
        self.visualization.calc_seperate_fermi_surface()
    
    def seperateCalcTraitMayavi(self, para_dict):
        self.visualization.calc_seperate_orbit(para_dict)
    
    def mousePressEvent(self, event):
        view = self.visualization.update_camera()
        self.camera.emit(view[0], view[1], view[2])
    
    def updateCameraMayavi(self, v1, v2, v3):
        self.visualization.scene.mlab.view(azimuth=v1, elevation=v2, distance=v3, focalpoint=np.array([0.5, 0.5, 0.5]))
    
    def bgMayavi(self, para_dict):
        self.visualization.update_bg(para_dict)
    
    def magSelectMayavi(self, bcell, mag):
        self.visualization.update_select_mag(bcell, mag)
    
    def edgeMayavi(self, para_dict):
        self.visualization.update_edge(para_dict)
    
    def gvCalcMayavi(self, bcell, ebands):
        self.visualization.update_calc_gv(bcell, ebands)
    
    def gvMayavi(self, para_dict):
        self.visualization.update_gv(para_dict)
    
    def synSliceMayavi(self, para_dict, mid_points):
        self.visualization.update_synSlice(para_dict, mid_points)
    
    def magSelectEndMayavi(self):
        self.visualization.update_select_mag_end()


