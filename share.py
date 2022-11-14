
__author__ = "Guanzhang Liu"
__maintainer__ = "Guanzhang Liu"
__email__ = "MG21340029@smail.nju.edu.cn"
__date__ = "October 23, 2022"

'''share.py为多个类共用调用的一些类与函数'''

DEFAULT_PATH = '/'
CALC_FILE = 'calc'
CONVAU2ANG = 0.529177209
'''
存储值类，解决pyqt中许多函数调用时，只能更改值，而无法记录值的问题
'''
import numpy as np

class saveVal(object):
    
    def __init__(self, v=0.):
        self.val = v
    
    def savVal(self, v):    # 保存值
        self.val = v
    
    def getVal(self):       # 返回值
        return self.val

'''
判断是否为数字函数
x：输入的object
'''
def notFloat(x):
    try:
        float(x)
        return False
    except:
        return True

def isFloat(x):
    try:
        float(x)
        return True
    except:
        return False

def cosangle(x, y):
    Lx=np.sqrt(x.dot(x))
    Ly=np.sqrt(y.dot(y))
    #相当于勾股定理，求得斜线的长度
    cos_angle=x.dot(y)/(Lx*Ly)
    #求得cos_sita的值再反过来计算，绝对长度乘以cos角度为矢量长度，初中知识。。
    angle=np.arccos(cos_angle)
    angle2=angle*360/2/np.pi
    #变为角度
    return angle2
    #x.dot(y) =  y=∑(ai*bi)

## 计算法向量（输入值为平面内两个向量，返回值为平面法向量）
def cal_normal(vec1, vec2):
    vec3_x = vec1[1] * vec2[2] - vec1[2] * vec2[1]
    vec3_y = vec1[2] * vec2[0] - vec1[0] * vec2[2]
    vec3_z = vec1[0] * vec2[1] - vec1[1] * vec2[0]
    return np.array([vec3_x, vec3_y, vec3_z])

## 计算投影向量（输入值为待求向量与平面法向量，返回值为待求向量的在平面的投影向量）
def shadow_vec(h, normal):
    h_shadow = h - (np.dot(h, normal) / np.sum(normal ** 2)) * normal
    return h_shadow


## 非转角测试部分
### 1.计算基矢b1与b2所构成平面的法向量
def calc_angle(bcell, h):
    kVectors_xyz = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    normal = cal_normal(kVectors_xyz[0], kVectors_xyz[1])

    ### 2.从而计算出h在b1与b2构成平面的投影向量h_shadow 
    h_shadow = shadow_vec(h, normal)

    ### 3.h_shadow与b3的角度为theta
    theta = cosangle(h_shadow, kVectors_xyz[0])

    ### 4.h与b3的角度为phi 
    phi = cosangle(h, kVectors_xyz[2])

    return theta, phi

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