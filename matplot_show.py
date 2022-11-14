import matplotlib
matplotlib.use("Qt5Agg")  # 声明使用pyqt5
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg  # pyqt5的画布
import matplotlib.pyplot as plt
# matplotlib.figure 模块提供了顶层的Artist(图中的所有可见元素都是Artist的子类)，它包含了所有的plot元素
from matplotlib.figure import Figure 
from pyface.qt import QtGui, QtCore
import numpy as np
from scipy.interpolate import splprep, splev
import pandas as pd
import os


class MyMatplotlibFigure(FigureCanvasQTAgg, QtGui.QWidget):
    """
    创建一个画布类，并把画布放到FigureCanvasQTAgg
    """
    def __init__(self, parent=None, dpi=50):
        # 创建一个Figure,该Figure为matplotlib下的Figure，不是matplotlib.pyplot下面的Figure
        self.figs = Figure(dpi=dpi)
        super(MyMatplotlibFigure, self).__init__(self.figs)  # 在父类种激活self.fig，
        # QtGui.QFrame.__init__(self, parent, self.figs)
        self.extremeData = []
        self.orbit = 0
        self.slice = 0
        self.orbit_data = []
        self.row = -1
        self.freq = 0

    def mat_plot_drow(self, t, s):
        """
        用清除画布刷新的方法绘图
        :return:
        """
        self.figs.clf()  # 清理画布，这里是clf()
        self.axes = self.figs.add_subplot(111)  # 清理画布后必须重新添加绘图区
        self.axes.patch.set_facecolor("#01386a")  # 设置ax区域背景颜色
        self.axes.patch.set_alpha(0.5)  # 设置ax区域背景颜色透明度
        self.figs.patch.set_facecolor('#01386a')  # 设置绘图区域颜色
        self.axes.spines['bottom'].set_color('r')  # 设置下边界颜色
        self.axes.spines['top'].set_visible(False)  # 顶边界不可见
        self.axes.spines['right'].set_visible(False)  # 右边界不可见
        # 设置左、下边界在（0，0）处相交
        # self.axes.spines['bottom'].set_position(('data', 0))  # 设置y轴线原点数据为 0
        self.axes.spines['left'].set_position(('data', 0))  # 设置x轴线原点数据为 0
        self.axes.plot(t, s, 'o-r', linewidth=0.5)
        self.figs.canvas.draw()  # 这里注意是画布重绘，self.figs.canvas
        self.figs.canvas.flush_events()  # 画布刷新self.figs.canvas
    
    def initSlice(self):
        self.figs.clf()
        self.axes = self.figs.add_subplot(111)        
        self.axes.set_aspect('equal', adjustable='box')
        self.axes.axis('off')

    def plotSlice(self, contour_2d, new_res, efermi, color):
        self.axes.tricontour(contour_2d[:,0], contour_2d[:,1], new_res, [efermi], colors=color) 
        self.figs.canvas.draw()  # 这里注意是画布重绘，self.figs.canvas
        self.figs.canvas.flush_events()  # 画布刷新self.figs.canvas

    def plotEdge(self, hull, points):
        for count, simplex in enumerate(hull.simplices):
            self.axes.plot(points[simplex, 0], points[simplex, 1], 'k-')  # 绘制边框
            # self.axes.scatter(points[hull.vertices[count], 0], points[hull.vertices[count], 1])
            self.figs.canvas.draw()  # 这里注意是画布重绘，self.figs.canvas
            self.figs.canvas.flush_events()  # 画布刷新self.figs.canvas
        # self.axes.set_xlim(np.min(points[hull.vertices][:,0])-0.2,np.max(points[hull.vertices][:,0])+0.2)
        # self.axes.set_ylim(np.min(points[hull.vertices][:,1])-0.2,np.max(points[hull.vertices][:,1])+0.2)

    def drawOrbit(self, row, orbit_data):
        
        self.row = row
        self.figs.clf()
        self.axes = self.figs.add_subplot(111)

        k_data = np.array(orbit_data[row])
        p1 = k_data[0]
        p2 = k_data[1]
        p3 = k_data[-1]
        a = (p2[1]-p1[1])*(p3[2]-p1[2])-(p3[1]-p1[1])*(p2[2]-p1[2])
        b = (p2[2]-p1[2])*(p3[0]-p1[0])-(p3[2]-p1[2])*(p2[0]-p1[0])
        c = (p2[0]-p1[0])*(p3[1]-p1[1])-(p3[0]-p1[0])*(p2[1]-p1[1])
        n = np.array([a, b, c])
        nz = n / np.sqrt(np.sum(n**2))
        nx = k_data[0] - k_data[-1]
        nx = nx / np.sqrt(np.sum(nx**2))
        ny_x = 1
        ny_z = (nx[0]*nz[1]-nz[0]*nx[1])/(nz[2]*nx[1]-nx[2]*nz[1])
        ny_y = (nz[0]*nx[2]-nx[0]*nz[2])/(nz[2]*nx[1]-nx[2]*nz[1])
        ny = np.array([ny_x, ny_y, ny_z])
        ny = ny / np.sqrt(np.sum(ny**2))
        t_array = np.c_[nx, ny, nz]
        test = np.dot(k_data, t_array)
        # fig = plt.figure(figsize=(2, 2), dpi=300)

        self.axes.set_aspect('equal', adjustable='box')
        # plt.plot(test[:, 0], test[:, 1])

        x = test[:, 0]
        y = test[:, 1]
        tck, u = splprep([x,y], s=0, per=True)
        xi, yi = splev(np.linspace(0, 1, 1000), tck)
        self.axes.fill(xi, yi, facecolor='c')
        # ax.invert_xaxis()
        self.axes.axis('off')
        # plt.savefig('/Users/wentworth/Desktop/test/' + '67_h3_992_smooth.png', dpi=300)
        self.figs.canvas.draw()  # 这里注意是画布重绘，self.figs.canvas
        self.figs.canvas.flush_events()  # 画布刷新self.figs.canvas
        self.orbit_data = np.c_[xi, yi]
    
    def drawExtreme(self, slice_orig, orbit_orig, df, accuracy=1, accuracy2=2):
        if 'orbit_match' in df.columns:
            del df['orbit_match']
        slice = slice_orig
        orbit = orbit_last_match = orbit_orig
        while (slice <= df['slice'].max()-1) & (not np.isnan(orbit_last_match)):
            eval_orbit = pd.Series(list(df.loc[(df['slice'] == slice) & (df['orbit'] == orbit), :].values[0]), list(df.loc[(df['slice'] == slice) & (df['orbit'] == orbit), :].columns))
            orbit_match = []
            for i, line in df[df['slice'] == slice+1].iterrows():
                if (line['avgx'] > (eval_orbit['avgx'] - eval_orbit['stdx'])) & (line['avgx'] < (eval_orbit['avgx'] + accuracy * eval_orbit['stdx'])) \
                        & (line['avgy'] > (eval_orbit['avgy'] - accuracy * eval_orbit['stdy'])) & (line['avgy'] < (eval_orbit['avgy'] + accuracy * eval_orbit['stdy'])) \
                        & (line['maxx'] > (eval_orbit['maxx'] - accuracy2 * eval_orbit['stdx'])) & (line['maxx'] < (eval_orbit['maxx'] + accuracy2 * eval_orbit['stdx'])) \
                        & (line['maxy'] > (eval_orbit['maxy'] - accuracy2 * eval_orbit['stdy'])) & (line['maxy'] < (eval_orbit['maxy'] + accuracy2 * eval_orbit['stdy'])) \
                        & (line['minx'] > (eval_orbit['minx'] - accuracy2 * eval_orbit['stdx'])) & (line['minx'] < (eval_orbit['minx'] + accuracy2 * eval_orbit['stdx'])) \
                        & (line['miny'] > (eval_orbit['miny'] - accuracy2 * eval_orbit['stdx'])) & (line['miny'] < (eval_orbit['miny'] + accuracy2 * eval_orbit['stdx'])):
                        orbit_match.append(line['orbit'])
            if len(orbit_match) == 0:
                orbit_last_match = np.nan
                break
            elif len(orbit_match) == 1:
                df.loc[((df['slice'] == slice) & (df['orbit'] == orbit)), 'orbit_match'] = orbit_match[0]
                orbit_last_match = orbit_match[0]
            else:
                print('activate')
                Bmin = 10000
                orbit_best_match = 0
                for i in orbit_match:
                    B = (line['avgx'] - eval_orbit['avgx']) ** 2 + (line['avgy'] - eval_orbit['avgy']) ** 2 + \
                        (line['maxx'] - eval_orbit['maxx']) ** 2 + (line['maxy'] - eval_orbit['maxy']) ** 2 + \
                        (line['minx'] - eval_orbit['minx']) ** 2 + (line['miny'] - eval_orbit['miny']) ** 2
                    if B < Bmin:
                        Bmin = B
                        orbit_best_match = i
                df.loc[((df['slice'] == slice) & (df['orbit'] == orbit)), 'orbit_match'] = orbit_best_match
                orbit_last_match = orbit_best_match
            slice = slice + 1
            orbit = orbit_last_match

        slice = slice_orig
        orbit = orbit_orig
        while (slice >= 2) & (not np.isnan(orbit)):
            eval_orbit = pd.Series(list(df.loc[(df['slice'] == slice) & (df['orbit'] == orbit), :].values[0]), list(df.loc[(df['slice'] == slice) & (df['orbit'] == orbit), :].columns))
            orbit_match = []
            for i, line in df[df['slice'] == slice-1].iterrows():
                if (line['avgx'] > (eval_orbit['avgx'] - accuracy * eval_orbit['stdx'])) & (line['avgx'] < (eval_orbit['avgx'] + eval_orbit['stdx'])) \
                        & (line['avgy'] > (eval_orbit['avgy'] - eval_orbit['stdy'])) & (line['avgy'] < (eval_orbit['avgy'] + eval_orbit['stdy'])) \
                        & (line['maxx'] > (eval_orbit['maxx'] - 2 * eval_orbit['stdx'])) & (line['maxx'] < (eval_orbit['maxx'] + 2 * eval_orbit['stdx'])) \
                        & (line['maxy'] > (eval_orbit['maxy'] - 2 * eval_orbit['stdy'])) & (line['maxy'] < (eval_orbit['maxy'] + 2 * eval_orbit['stdy'])) \
                        & (line['minx'] > (eval_orbit['minx'] - 2 * eval_orbit['stdx'])) & (line['minx'] < (eval_orbit['minx'] + 2 * eval_orbit['stdx'])) \
                        & (line['miny'] > (eval_orbit['miny'] - 2 * eval_orbit['stdx'])) & (line['miny'] < (eval_orbit['miny'] + 2 * eval_orbit['stdx'])):
                        orbit_match.append(line['orbit'])
            if len(orbit_match) == 0:
                orbit_last_match = np.nan
                break
            elif len(orbit_match) == 1:
                df.loc[((df['slice'] == slice-1) & (df['orbit'] == orbit_match[0])), 'orbit_match'] = orbit
                orbit_last_match = orbit_match[0]
            else:
                print('activate')
                Bmin = 10000
                orbit_best_match = 0
                for i in orbit_match:
                    B = (line['avgx'] - eval_orbit['avgx']) ** 2 + (line['avgy'] - eval_orbit['avgy']) ** 2 + \
                        (line['maxx'] - eval_orbit['maxx']) ** 2 + (line['maxy'] - eval_orbit['maxy']) ** 2 + \
                        (line['minx'] - eval_orbit['minx']) ** 2 + (line['miny'] - eval_orbit['miny']) ** 2
                    if B < Bmin:
                        Bmin = B
                        orbit_best_match = i
                df.loc[((df['slice'] == slice-1) & (df['orbit'] == orbit_best_match)), 'orbit_match'] = orbit
                orbit_last_match = orbit_best_match
            slice = slice - 1
            orbit = orbit_last_match

        slice = slice_orig
        orbit = orbit_orig
        total = pd.DataFrame()
        slice_now = slice_latter = slice
        orbit_now = orbit_latter = orbit
        while slice_now >= 1:
            if len(df[(df['slice'] == slice_now-1) & (df['orbit_match'] == orbit_now)]) == 0:
                break
            slice_now = slice_now - 1
            orbit_now = df[(df['slice'] == slice_now) & (df['orbit_match'] == orbit_now)]['orbit'].values[0]
            total = pd.concat([total, df[(df['slice'] == slice_now) & (df['orbit'] == orbit_now)]])
        while slice_latter <= df['slice'].max():
            total = pd.concat([total, df[(df['slice'] == slice_latter) & (df['orbit'] == orbit_latter)]])
            orbit_latter = df[(df['slice'] == slice_latter) & (df['orbit'] == orbit_latter)]['orbit_match'].values[0]
            slice_latter = slice_latter + 1
            if np.isnan(orbit_latter):
                break
        
        self.figs.clf()
        self.axes = self.figs.add_subplot(111)

        self.axes.plot(np.arange(total['slice'].min(), total['slice'].max()+1), total.sort_values(by=['slice', 'orbit'])['freq'].values)
        self.axes.scatter(slice, df[(df['slice'] == slice) & (df['orbit'] == orbit)]['freq'], c='red')
        self.freq = df[(df['slice'] == slice) & (df['orbit'] == orbit)]['freq'].values[0]

        self.figs.canvas.draw()  # 这里注意是画布重绘，self.figs.canvas
        self.figs.canvas.flush_events()  # 画布刷新self.figs.canvas

        self.extremeData = np.c_[np.arange(total['slice'].min(), total['slice'].max()+1), total.sort_values(by=['slice', 'orbit'])['freq'].values]
        self.slice = slice
        self.orbit = orbit
        return total, slice, orbit
    
    def saveExtremeData(self, file_path):
        if len(self.extremeData) != 0:
            f = '## Extreme value data file generated by Skeaf_demo\n## Slice ' + str(self.slice) + ', orbit '+str(self.orbit)+', freq'+str(self.freq)+'\n## slice,freq\n\n'
            saveData = self.extremeData
            with open(os.path.join(file_path, 'Ext_'+str(self.freq)+'.csv'), 'w') as filewriter:
                filewriter.write(f)
            with open(os.path.join(file_path, 'Ext_'+str(self.freq)+'.csv'), 'a') as filewriter:
                np.savetxt(filewriter, saveData, delimiter=',')
    
    def saveOrbitData(self, file_path, freq_list):
        if self.row != -1:
            f = '## Orbit data file generated by Skeaf_demo\n## Freq '+str(freq_list[self.row])+'\n## x,y\n\n'
            saveData = self.orbit_data
            with open(os.path.join(file_path, 'Orb_'+str(freq_list[self.row])+'.csv'), 'w') as filewriter:
                filewriter.write(f)
            with open(os.path.join(file_path, 'Orb_'+str(freq_list[self.row])+'.csv'), 'a') as filewriter:
                np.savetxt(filewriter, saveData, delimiter=',') 

