# -*- coding: utf-8 -*-
###LCEM:Light condition estimation model###
"""""
Virsion: 4.0
Developer: Wang Bin
Latest modification time: 2019-6-16
"""""

import numpy as np
from laspy.file import File
import pandas as pd
import time

inFile = File("E:\\实验数据\\林下光环境\\light_test\\LCEM\\LCEM V 4\\DATA\\"
              "3rd_正式数据DS0.1_转换为Las.las", mode = "r")

inFile_topo = File("E:\\实验数据\\林下光环境\\light_test\\LCEM\\LCEM V 4\\DATA\\"
                   "Topo_points.las", mode = "r")

inFile_light_feild = np.loadtxt("E:\\实验数据\\林下光环境\\light_test\\LCEM\\LCEM V 4\\DATA\\"
                                "光环境验证数据_正式_海拔.csv", dtype = np.float, delimiter = ',', skiprows= 1)
np.set_printoptions(suppress=True)
#===========================================================#
start_time = time.time()
def scaled_x_dimension(las_file):
    x_dimension = las_file.X
    scale = las_file.header.scale[0]
    offset = las_file.header.offset[0]
    return (x_dimension*scale + offset)
X = scaled_x_dimension(inFile)

def scaled_y_dimension(las_file):
    y_dimension = las_file.Y
    scale = las_file.header.scale[1]
    offset = las_file.header.offset[1]
    return (y_dimension*scale + offset)
Y = scaled_y_dimension(inFile)

def scaled_z_dimension(las_file):
    z_dimension = las_file.Z
    scale = las_file.header.scale[2]
    offset = las_file.header.offset[2]
    return (z_dimension*scale + offset)
Z = scaled_z_dimension(inFile)

def read_classification(las_file):
    class_number = las_file.raw_classification
    return (class_number)
CLASSNUM = read_classification(inFile)

def coords():
    coords_data = np.vstack((X, Y, Z)).transpose()
    return (coords_data)
all_data = coords ()
end_time = time.time()
print("1-读取雷达数据：", end_time - start_time, "S|||点个数：", len(all_data))
#===========================================================#

#===========================================================#
#
def extract_topography_points(las_File):
    J_ETP = CLASSNUM == 2
    kept_points_ETP = las_File[J_ETP]
    return (kept_points_ETP)

end_time = time.time()
print("2-抽取地形点：", end_time - start_time, "S|||点个数：", len(ETP))
#==============================>#==============================>#==============================>
Number_of_point = 0
#==============================>#==============================>#==============================>
def extract_plot_points(las_File, Number_of_point, box_length):
    TCS_X = inFile_light_feild[Number_of_point][0]
    TCS_Y = inFile_light_feild[Number_of_point][1]
    B_X_max = TCS_X + box_length
    B_X_min = TCS_X - box_length
    B_Y_max = TCS_Y + box_length
    B_Y_min = TCS_Y - box_length
    J_X = np.logical_and(las_File[:, 0] <= B_X_max,
                         las_File[:, 0] >= B_X_min)
    J_Y = np.logical_and(las_File[:, 1] <= B_Y_max,
                         las_File[:, 1] >= B_Y_min)
    J_XY = np.logical_and(J_X, J_Y)
    kept_points_EPP = las_File[J_XY]
    return (kept_points_EPP)
EPP = extract_plot_points(all_data, Number_of_point, 50)



end_time = time.time()
print("3-有效范围筛选后：", end_time - start_time, "S|||点个数：", len(EPP))
#===========================================================#
#np.savetxt("E:\\实验数据\\林下光环境\\light_test\\LCEM\\LCEM V 4\\results\\"
#           "EPP_77.csv", EPP, delimiter=',')
#===========================================================#
#
def Merge_data(EPP_in, ETP_in):
    MD = np.vstack((EPP_in, ETP_in))
    return (MD)
Merge_data_fin = Merge_data(EPP, ETP)
end_time = time.time()
print("4-数据合并后：", end_time - start_time, "S|||点个数：", len(Merge_data_fin))
#===========================================================#

#===========================================================#
#transfer coordinate system(TCS)
def TCS(Number_of_point):
    TCS_X = inFile_light_feild[Number_of_point][0]
    TCS_Y = inFile_light_feild[Number_of_point][1]
    TCS_Z = inFile_light_feild[Number_of_point][5]
    #TCS_X = (X.max() + X.min()) / 2
    #TCS_Y = (Y.max() + Y.min()) / 2
    #TCS_Z = (Z.max() + Z.min()) / 2
    oringin_point = [TCS_X, TCS_Y, TCS_Z]
    TCS_result = Merge_data_fin - oringin_point
    return (TCS_result)
TCS_O = TCS(Number_of_point)
end_time = time.time()
print("5-坐标转移后：", end_time - start_time, "S|||点个数：", len(TCS_O))
#===========================================================#

#===========================================================#
#filter points which greater than height of obtain.
# Filter out points below the height of the observation point.(FBO)
def FBO():
    height_OP  = TCS_O[:, 2]
    GT0 = height_OP >= 0
    GT0_kept = TCS_O[GT0]
    return (GT0_kept)
FBO = FBO()
end_time = time.time()
print("6-高度过滤后：", end_time - start_time, "S|||点个数：", len(FBO))
#===========================================================#

#===========================================================#
#
def CCCtoSPC(point_kept):
    X_tra = point_kept[:, 0]
    Y_tra = point_kept[:, 1]
    Z_tra = point_kept[:, 2]
    r = np.sqrt(np.sum((point_kept) ** 2, axis=1))
    theta = np.arccos(Z_tra / (r))
    phi = np.arctan(Y_tra / X_tra)
    # print("r:", r, "theta:", theta, "phi:", phi)
    for index in range(len(point_kept)):
        if point_kept[index][0] >= 0:
            phi[index] = phi[index]
        else:
            phi[index] = phi[index] + np.pi
    for index_exchange in range(len(r)):
        r[index_exchange] = 100  # 投影球面半径无论是多大都对结果不产生影响,只是设定一个值用于创建球面
    # -------------------------------------------------------------
    ass_r = np.vstack((theta, phi, r)).transpose()
    #data_Frame_filter = pd.DataFrame(ass_r, columns=['theta', 'phi', 'r'])
    #filter_duplicate_points = data_Frame_filter.drop_duplicates(keep='first')
    #FDP_array = np.array(filter_duplicate_points)
    #球面过滤的意义不大，并不能过滤很多重复点，主要原因可能是由于带点在空间中并不占据空间所致
    return (ass_r)
CCCtoSPC_FIN = CCCtoSPC(FBO)
end_time = time.time()
print("7-直角坐标系转球极坐标系(已球面过滤)：", end_time - start_time, "S|||点个数：", len(CCCtoSPC_FIN))
#===========================================================#

#===========================================================#
#
def createSphereFishnet(kept_point_on_Sphere):
    theta_Origin = np.pi / 2
    phi_Origin = 0
    r_Origin = 100
    listNumbers1 = []
    listNumbers1_labletheta = []
    list_tra1_CSF= []
    listNumbers2 = []
    listNumbers2_labletheta = []
    list_tra2_CSF = []
    for index_J_theta in range(90):
        for index_J_phi1 in range(180):
            theta_polygon_min = theta_Origin - (index_J_theta * (np.pi / 180))
            phi_polygon_min = phi_Origin + (index_J_phi1 * (np.pi / 180))
            theta_polygon_max = theta_Origin - ((index_J_theta + 1) * (np.pi / 180))
            phi_polygon_max = phi_Origin + ((index_J_phi1 + 1) * (np.pi / 180))
            # ===============>
            theta_J1 = np.logical_and(kept_point_on_Sphere[:, 0] <= theta_polygon_min,
                                     kept_point_on_Sphere[:, 0] >= theta_polygon_max)
            phi_J1 = np.logical_and(kept_point_on_Sphere[:, 1] >= phi_polygon_min,
                                   kept_point_on_Sphere[:, 1] <= phi_polygon_max)
            thetaphi_J1 = np.logical_and(theta_J1, phi_J1)
            keep_points1 = kept_point_on_Sphere[thetaphi_J1]
            listNumbers1 = np.append(listNumbers1, len(keep_points1))
            listNumbers1_labletheta = np.append(listNumbers1_labletheta, index_J_theta)
            list_tra1_CSF = np.vstack((listNumbers1_labletheta, listNumbers1)).transpose()
        for index_J_phi2 in range(180):
            theta_polygon_min = theta_Origin - (index_J_theta * (np.pi / 180))
            phi_polygon_min = phi_Origin - (index_J_phi2 * (np.pi / 180))
            theta_polygon_max = theta_Origin - ((index_J_theta + 1) * (np.pi / 180))
            phi_polygon_max = phi_Origin - ((index_J_phi2 + 1) * (np.pi / 180))
            # ===============>
            theta_J2 = np.logical_and(kept_point_on_Sphere[:, 0] <= theta_polygon_min,
                                     kept_point_on_Sphere[:, 0] >= theta_polygon_max)
            phi_J2 = np.logical_and(kept_point_on_Sphere[:, 1] <= phi_polygon_min,
                                   kept_point_on_Sphere[:, 1] >= phi_polygon_max)
            thetaphi_J2 = np.logical_and(theta_J2, phi_J2)
            keep_points2 = kept_point_on_Sphere[thetaphi_J2]
            listNumbers2 = np.append(listNumbers2, len(keep_points2))
            listNumbers2_labletheta = np.append(listNumbers2_labletheta, index_J_theta)
            list_tra2_CSF = np.vstack((listNumbers2_labletheta, listNumbers2)).transpose()
        end_time = time.time()
        print("8-Sphere_processing: ", (((index_J_theta + 1) / (90)) * 100), '%|||用时：', end_time - start_time)
    list_Numbers_all = np.vstack((list_tra1_CSF, list_tra2_CSF))
    return (list_Numbers_all)
CSF = createSphereFishnet(CCCtoSPC_FIN)
#np.savetxt("E:\\实验数据\\林下光环境\\light_test\\LCEM\\LCEM V 4\\results\\"
           #"CSF.csv", CSF, delimiter=',')
#===========================================================#

#===========================================================#
#
def Binarization(voxel_statistic):
    for index_Binarization in range(len(voxel_statistic)):
        if voxel_statistic[index_Binarization][1] > 1:
            voxel_statistic[index_Binarization][1] = 1
        else:
            voxel_statistic[index_Binarization][1] = 0

    for index_Binarization_degree in range(len(voxel_statistic)):
        if voxel_statistic[index_Binarization_degree][0] < 45:
            voxel_statistic[index_Binarization_degree][1] = 1
        else:
            voxel_statistic[index_Binarization_degree][1] = voxel_statistic[index_Binarization_degree][1]
    return (voxel_statistic)
Binarization_statistic = Binarization(CSF)
#===========================================================#

#===========================================================#
#
def Area_circle(h, R):
    S = 2 * np.pi * R * h
    return (S)

def height(a1, a2, R):
    h1 = np.sin(a1) * R
    h2 = np.sin(a2) * R
    H = h2-h1
    return (H)

def Area_circle_list():
    Area_circle_list = []
    list_Area_crown = []
    labels = []
    for index_a in range(89):
        theta_a1 = 0 + ((index_a) * (np.pi / 180))
        theta_a2 = 0 + ((index_a + 1) * (np.pi / 180))
        H = height(theta_a1, theta_a2, 100)
        S = Area_circle(H, 100)
        Area_circle_list = np.append(Area_circle_list, S)
        labels = np.append(labels, index_a)
    list_circle_area = np.vstack((labels, Area_circle_list)).transpose()
    sum_circle = list_circle_area[:, 1]
    Area_crown = (2 * np.pi * 10000) - (np.sum(sum_circle))
    list_Area_crown = [89, Area_crown]
    list_circle_area_all = np.vstack((list_circle_area, list_Area_crown))
    return (list_circle_area_all)
Area_circle_list_fin = Area_circle_list()
#===========================================================#

#===========================================================#
#
def openness_calc():
    for index_CAP in range(len(Binarization_statistic)):
        I = int(Binarization_statistic[index_CAP][0])
        #print(I)
        J = Binarization_statistic[index_CAP][1]
        #print(J)
        Binarization_statistic[index_CAP][1] = (Area_circle_list_fin[I][1] * J) / 360
    Area_cell_sum = np.sum(Binarization_statistic[:, 1])
    CAP = ((2 * np.pi * 10000) - Area_cell_sum) / (2 * np.pi * 10000)
    return (CAP)
CAP = openness_calc()
end_time = time.time()
print("9-计算完成时间：", end_time - start_time, "S|||Canopy_Openness：", CAP)
#===========================================================#
"""""
def SPCtoCCC(FDP_points):
    theta_S = FDP_points[:, 0]
    phi_S = FDP_points[:, 1]
    r_S = FDP_points[:, 2]
    X_to_CCC = (r_S) * (np.sin(theta_S)) * (np.cos(phi_S))
    Y_to_CCC = (r_S) * (np.sin(theta_S)) * (np.sin(phi_S))
    Z_to_CCC = (r_S) * (np.cos(theta_S))
    coods_toCCC = np.vstack((X_to_CCC, Y_to_CCC, Z_to_CCC)).transpose()
    return (coods_toCCC)
coods_toCCC = SPCtoCCC(CCCtoSPC_FIN)
#print("转换回直角坐标:", len(coods_toCCC))
np.savetxt("E:\\实验数据\\林下光环境\\light_test\\LCEM\\LCEM V 4\\results\\"
           "C2C_77.csv", coods_toCCC, delimiter=',')
"""""
#===========================================================#
