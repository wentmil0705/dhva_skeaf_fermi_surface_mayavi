__author__ = "Guanzhang Liu"
__maintainer__ = "Guanzhang Liu"
__email__ = "MG21340029@smail.nju.edu.cn"
__date__ = "October 23, 2022"

'''
slice_ds.py为ui_slice.py界面对应功能的设计文件
'''

from ui_py_files.ui_slice import *
from PyQt5 import QtWidgets, QtCore
from share import *
from scipy.interpolate import interpn
import re
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scipy.spatial import Voronoi

class mySlice(QtWidgets.QMainWindow, Ui_slice_window):
    _midSignal = QtCore.pyqtSignal(list)    # 传递截面位置信息

    def __init__(self, parent=None):
        super(mySlice, self).__init__(parent)
        self.setupUi(self)
        self.BXSF_FILE_PATH_list = []
        self.MESH_list = []
        self.EBANDS_list = []
        self.BCELL = None
        self.EFERMI = 0.
        self.BXSF_list = []
        self.ORBIT_DATA = []
        self.ORBIT_PRED = []
        self.bz = np.array([1,1,1])
        self.THETA = 0
        self.PHI = 0
        self.sc = None
        self.region_list = []
        self.min_p = None
        self.max_p = None

        self.slice_actionimport_bxsf_files.triggered.connect(self.showBxsf)
        self.slice_actionimport_results_orbitoutlines_files.triggered.connect(self.showOrbOut)
        self.Theta_line.editingFinished.connect(self.showSc)
        self.Phi_line.editingFinished.connect(self.showSc)
        self.bxsf_checkableComboBox.activated.connect(self.showSlice)
        self.orbits_ComboBox.activated.connect(self.showSlice)
        self.BZ1_line.setText(str(self.bz[0]))
        self.BZ2_line.setText(str(self.bz[1]))
        self.BZ3_line.setText(str(self.bz[2]))
        self.BZ_button.clicked.connect(self.showBZ)
        self.mid_slider.setMinimum(0)
        self.mid_slider.setMaximum(100)
        self.mid_slider.setSingleStep(1)
        self.mid_slider.valueChanged.connect(self.showMidSlice)
        self.Theta_line.setText(str(0.))
        self.Phi_line.setText(str(0.))
        self.showSc()
    
    # 显示滑块所选位置切面
    def showMidSlice(self):        
        try:
            mid_points = self.mid_slider.value() / 100 * (self.max_p - self.min_p) + self.min_p
            self.calcSlice(mid_points)
            if self.actionsync.isChecked():
                self._midSignal.emit(list(mid_points))
        except:
            pass
    
    # 布里渊区的获取
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
    
    # 更新参数函数
    def showBZ(self):
        if self.BCELL is not None:
            if isFloat(self.BZ1_line.text()) & isFloat(self.BZ2_line.text()) & isFloat(self.BZ3_line.text()):
                self.bz = np.array([int(self.BZ1_line.text()), int(self.BZ2_line.text()), int(self.BZ3_line.text())])
            bz_max = np.max(self.bz)
            p, l, f = self.get_brillouin_zone_3d(self.BCELL)
            px, py, pz = np.tensordot(
                self.BCELL,
                np.mgrid[-1:bz_max+1, -1:bz_max+1, -1:bz_max+1],
                axes=[0, 0]
            )
            points = np.c_[px.ravel(), py.ravel(), pz.ravel()]
            self.tree = cKDTree(points)
            gamma_region_id = self.tree.query([0, 0, 0])[1]
            self.region_list = [gamma_region_id]
            bz_list = []
            for i in range(self.bz[0]):
                for j in range(self.bz[1]):
                    for k in range(self.bz[2]):
                        bz_list.append([i, j, k])
            bz_list = bz_list[1:]

            dm = int(np.cbrt(np.max(self.tree.indices) + 1))
            test = np.reshape(np.arange(int(dm**3)), (dm, dm, dm))
            self.all_l = [l]
            for bz_array in bz_list:
                new_l = []
                for part in l:
                    new_part = np.dot(part, np.linalg.inv(self.BCELL)) + bz_array
                    new_part = np.dot(new_part, self.BCELL)
                    new_l.append(new_part)
                self.all_l.append(new_l)
                new_id_loc = np.asarray(np.where(test == gamma_region_id))[:, 0] + bz_array
                self.region_list.append(test[new_id_loc[0]][new_id_loc[1]][new_id_loc[2]])
            
            l_list = []
            for p_l in self.all_l:
                for x in p_l:
                    l_list += list(map(lambda m: list(m), x))

            l_list = np.array(l_list)
            l_list_sc = np.dot(l_list, np.linalg.inv(self.sc))
            min_p = l_list_sc[np.argmin(l_list_sc[:,2])]
            max_p = l_list_sc[np.argmax(l_list_sc[:,2])]
            max_p = np.r_[min_p[0:2], max_p[2]] 
            self.min_p = np.dot(min_p, self.sc)
            self.max_p = np.dot(max_p, self.sc)

    # 显示切面形状  
    def showSlice(self):        
        self.showBZ()
        cur_id = self.orbits_ComboBox.currentIndex()
        if cur_id != -1:
            orbit_data = self.ORBIT_DATA[cur_id]
            mid_points = np.mean(orbit_data, axis=0)
            self.calcSlice(mid_points)
            if self.actionsync.isChecked():
                self._midSignal.emit(list(mid_points))
        else:
            self.showMidSlice()
    
    # 计算切面形状
    def calcSlice(self, mid_points):
        bz_max = np.max(self.bz)
        mid_points_sc = np.dot(mid_points, np.linalg.inv(self.sc*np.sqrt(np.sum(self.BCELL**2, axis=1))))
        x = np.linspace(-1*bz_max, 1*bz_max, 80*bz_max*2)
        y = np.linspace(-1*bz_max, 1*bz_max, 80*bz_max*2)
        X,Y = np.meshgrid(x,y)
        XY = np.concatenate((X[:,:,None], Y[:,:,None]), axis=-1).reshape(-1,2)
        XYZ = np.c_[XY, [mid_points_sc[2]]*len(XY)]
        slice_kgrid = np.dot(XYZ, self.sc*np.sqrt(np.sum(self.BCELL**2, axis=1)))
        find_energy_list = np.dot(slice_kgrid, np.linalg.inv(self.BCELL))
        interpol_list = np.fromiter(map(lambda x: x-(x // 1) if ((x > 1) | (x < 0)) else x, find_energy_list[:, 0]), dtype=float)
        interpol_list = np.c_[interpol_list, np.fromiter(map(lambda x: x-(x // 1) if ((x > 1) | (x < 0)) else x, find_energy_list[:, 1]),dtype=float)]
        interpol_list = np.c_[interpol_list, np.fromiter(map(lambda x: x-(x // 1) if ((x > 1) | (x < 0)) else x, find_energy_list[:, 2]),dtype=float)]
                
        self.slice_widget.initSlice()
        colors = ['brown', 'orangered', 'c', 'purple', 'steelblue']
        for r in self.region_list:
            contour_2d_id = self.tree.query(slice_kgrid)[1]
            contour_2d_bz = np.array(contour_2d_id == r)

            contour_2d = np.dot(slice_kgrid[contour_2d_bz], np.linalg.inv(self.sc * np.sqrt(np.sum(self.BCELL**2, axis=1))))[:,0:2]
            if len(contour_2d) != 0:
                for i, ebands in enumerate(self.EBANDS_list):
                    mesh = self.MESH_list[i]
                    ebands = np.array(ebands).reshape(mesh[0], mesh[1], mesh[2])
                    values = np.pad(ebands, (0,1), mode='wrap')
                    length = np.array(values.shape, dtype=int)
                    new_edge = 1 + 1 / (length-1)
                    points = (np.linspace(0, new_edge[0], length[0]), np.linspace(0, new_edge[1], length[1]), np.linspace(0, new_edge[2], length[2]))
                    res = interpn(points, values, interpol_list, method='linear', fill_value=True, bounds_error=False) 
                    new_res = res[contour_2d_bz]
                    self.slice_widget.plotSlice(contour_2d, new_res, self.EFERMI, colors[i])
        if self.line_checkbox.isChecked():
            try:
                h_slice = np.r_[self.sc[2], -np.dot(mid_points, self.sc[2])]
                for p_l in self.all_l:
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
                        points_3d = np.dot(points, np.linalg.inv(self.sc * np.sqrt(np.sum(self.BCELL**2, axis=1))))    
                        points_2d = points_3d[:,0:2]
                        hull = ConvexHull(points_2d)

                        self.slice_widget.plotEdge(hull, points_2d)
            except:
                pass

    
    # 更新磁场方向
    def showSc(self):
        if isFloat(self.Theta_line.text()) & isFloat(self.Phi_line.text()):
            self.THETA = float(self.Theta_line.text())
            self.PHI = float(self.Phi_line.text())
            self.sc = self.ruc2sc(self.PHI, self.THETA)
            try:
                self.showBZ()
                self.showSlice()
            except:
                pass
    
    # 对输入的bxsf文件进行处理
    def showBxsf(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', DEFAULT_PATH)

        if fname[0] != '':
            bcell, ebands, efermi, mesh = self.interpol(fname[0])
            if self.BCELL is not None:
                if (np.sum(self.BCELL != bcell) != 0) | (np.sum(self.EFERMI != efermi) != 0):
                    QtWidgets.QMessageBox.critical(self, "Error", "Multi-bxsf files are not consistent!")
            else:
                self.BCELL = bcell
                self.EFERMI = efermi
            self.EBANDS_list.append(ebands)
            self.MESH_list.append(mesh)
            if self.BXSF_list == []:
                self.BXSF_list = ['All']
                self.bxsf_checkableComboBox.addItem('All')
            self.BXSF_FILE_PATH_list.append(fname[0])
            self.BXSF_list.append(fname[0].split('/')[-1].replace('.bxsf', ''))
            self.bxsf_checkableComboBox.addItem(self.BXSF_list[-1])
            self.showBZ()

    # 对输入的轨道文件进行处理
    def showOrbOut(self):
        fname = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', DEFAULT_PATH)

        if fname[0] != '':
            orbit_data, orbit_pred_freqs, theta, phi = self.read_au_Ang_file(fname[0])
            if self.THETA is not None:
                if self.THETA != theta:
                    QtWidgets.QMessageBox.critical(self, "Error", "Mag direction are not consistent!")
                else:
                    self.THETA = theta
                    self.PHI = phi
                    self.Theta_line.setText(str(theta))
                    self.Phi_line.setText(str(phi))
                    self.sc = self.ruc2sc(self.PHI, self.THETA)
            orbit_data = self.move_bz_orbit_data(self.bz, self.BCELL, orbit_data)
            self.ORBIT_DATA += orbit_data          
            self.ORBIT_PRED += orbit_pred_freqs
            for i in orbit_pred_freqs:                      
                self.orbits_ComboBox.addItem(i)

    # 插值处理
    def interpol(self, file_path):
        with open(file_path, 'r') as filereader:
            bxsf_orig = filereader.read()
        fmt_1 = re.findall(r'BANDGRID_3D_BANDS[\d\s\.\+\-]+BAND', bxsf_orig)
        fmt_2 = re.findall(r'BEGIN_BANDGRID_3D[\d\.\s\+\-]+BAND', bxsf_orig)
        fmt = fmt_1 if fmt_1 != [] else fmt_2
        fmt = fmt[0].split('\n')
        mesh = np.fromiter(map(lambda x: int(x), fmt[2].split()), dtype=int)
        bcell = fmt[4:7]
        bcell = np.array(list(map(lambda line: list(map(lambda x: float(x), line.split())), bcell)))

        ## 读取所有采样点的能量值
        fermi_ebands3d_uc = re.findall(r'BAND:\s+[\d+\.\s\-\+e]+', bxsf_orig)[0]
        fermi_ebands3d_uc = list(map(lambda x: float(x), fermi_ebands3d_uc.split()[2:]))
        efermi = float(re.findall(r'Fermi Energy:\s+[\d\.\+\-]+', bxsf_orig)[0].split()[-1])

        return bcell, fermi_ebands3d_uc, efermi, mesh
    
    # 磁场方向计算处理
    def ruc2sc(self, PHI, THETA):
        phi = np.radians(PHI)
        theta = np.radians(THETA)
        s = np.sin(phi)
        t = np.cos(phi)
        u = 1 - np.cos(phi)
        v = np.sin(theta)
        w = np.cos(theta)

        sc_to_ruc = np.array([[v*v*u+t, -v*w*u, -w*s],
                            [-v*w*u, w*w*u+t, -v*s],
                            [w*s, v*s, t]
                            ])
        
        kVectors_xyz = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        sc = np.dot(sc_to_ruc, kVectors_xyz)
        return sc
    
    # 轨道文件读取
    def read_au_Ang_file(self, file_path):
        with open(file_path, 'r') as filereader:
            f = filereader.read()

        # 处理没有极值频率的情况（此时文件为空）
        try:                     
            orbit_data = re.findall(r'kx[\sE\-\.\dkyz\+]+', f)
            orbit_data  = list(map(lambda data: list(map(lambda x: list(map(float, x.split())), data.split('\n')[1:-1])), orbit_data))
            # self.ORBIT_DATA += orbit_data
            orbit_pred_freqs = list(map(lambda x: x+'('+self.BXSF_list[-1]+')', re.findall(r'Freq\(kT, average of all copies\)\s+\=\s+(\d+\.\d+)', f)))
            if 'Ang' in file_path:
                orbit_data = list(map(lambda x: np.array(x)*CONVAU2ANG, orbit_data))
            orbit_data = list(map(lambda x: x/(2*np.pi), orbit_data))
            theta = float(re.findall(r'Theta\(deg\)\s+=\s+(\d+\.\d+)', f)[0])
            phi = float(re.findall(r'Phi\(deg\)\s+=\s+(\d+\.\d+)', f)[0])
            return orbit_data, orbit_pred_freqs, theta, phi

        except:
            QtWidgets.QMessageBox.warning(self, "Warning", "au/Ang file has problems, maybe because of empty results.")

    # 将轨道迁移至第一布里渊区
    def move_bz_orbit_data(self, bz, bcell, orbit_data_list):

        new_orbit_data = []
        bz_max = np.max(bz)
        px, py, pz = np.tensordot(
            bcell,
            np.mgrid[-1:bz_max+1, -1:bz_max+1, -1:bz_max+1],
            axes=[0, 0]
        )
        points = np.c_[px.ravel(), py.ravel(), pz.ravel()]
        tree = cKDTree(points)
        for orbit_data in orbit_data_list:
            orig_data = np.dot(orbit_data, np.linalg.inv(bcell))
            orbit_id = tree.query(np.dot(np.mean(orig_data, axis=0), bcell))[1]
            dm = int(np.cbrt(np.max(tree.indices) + 1))
            test = np.reshape(np.arange(int(dm**3)), (dm, dm, dm))
            gamma_region_id = tree.query([0, 0, 0])[1]
            moveArray = np.asarray(np.where(test == gamma_region_id))[:, 0] - np.asarray(np.where(test == orbit_id))[:, 0]
            orig_data = orig_data + moveArray
            new_data = np.dot(orig_data, bcell)
            new_orbit_data.append(new_data)
        return new_orbit_data


    




    

    
