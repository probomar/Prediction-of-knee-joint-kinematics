import numpy as np
import os
import pyvista as pv
import pandas as pd
import math as m
from scipy.spatial.transform import Rotation as R


# create direction in parent folder
def create_direction(path_parent, direction):
    path = os.path.join(path_parent, direction)
    if not os.path.exists(path):
        os.makedirs(path)


# Set folder for results
folder_name1 = 'Results'

# Set model of ligaments
# ligament = 'linear'
ligament = 'non-linear'
folder_name = folder_name1 + '_' + ligament

step = 5  # 'time' step
z_step = 0.02
fiy_step = 0.05


# Optimalization parameters
Fr = 10  # N Residual force
Mr = 5000  # Nmm Residual torque

step_F0 = 0.0005  # Initial force step
step_M0 = 0.00000005  # Initial torque step
alpha = 1.1  # Step magnification factor
beta = 0.5  # Step reduction factor
step_F_min = 1e-04  # Minimal force step
n_max = 50  # Maximum number of iterations

# Mechanical properties
E_cartilage = 10.35  # MPa Young elastic modul of cartilage
ny_cartilage = 0.3  # Poisson's ratio of cartilage
E_meniscus = 10.35  # MPa Young elastic modul of meniscus
ny_meniscus = 0.3  # Poisson's ratio of cartilage
k_bone = 100000  # N/mm^3  Bone stiffness

# Ligaments
epsL = 0.03  # Linear range threshold ligament strain

# Ligaments stiffness (1 - linear, 2 -  non-linear)
k1ACLa = 83.15
k2ACLa = 22.48
k1ACLp = 83.15
k2ACLp = 26.27

k1PCLa = 125.0
k2PCLa = 31.26
k1PCLp = 60.0
k2PCLp = 19.29

k1LCL = 72.22
k2LCL = 10.0

k1MCLa = 91.25
k2MCLa = 10.0
k1MCLo = 27.86
k2MCLo = 5.0
k1MCLd = 21.07
k2MCLd = 5.0

# create folders and file for results
create_direction(folder_name, 'flex')
create_direction(folder_name + '/stress', 'data')
create_direction(folder_name + '/stress', 'figure')
create_direction(folder_name + '/stress_cartilage', 'data')
create_direction(folder_name + '/stress_cartilage', 'figure')
create_direction(folder_name + '/stress_meniscus', 'data')
create_direction(folder_name + '/stress_meniscus', 'figure')
create_direction(folder_name + '/stress_bone', 'data')
create_direction(folder_name + '/stress_bone', 'figure')
create_direction(folder_name + '/cartilage_femoral', 'data')
create_direction(folder_name + '/cartilage_femoral', 'figure')
create_direction(folder_name + '/cartilage_tibial', 'data')
create_direction(folder_name + '/cartilage_tibial', 'figure')
create_direction(folder_name + '/distance_cartilage0', 'data')
create_direction(folder_name + '/distance_cartilage0', 'figure')
create_direction(folder_name + '/distance_meniscus0', 'data')
create_direction(folder_name + '/distance_meniscus0', 'figure')
create_direction(folder_name, 'ligaments')
create_direction(folder_name, 'results')

cor0 = folder_name + '/cor0.csv'
if os.path.exists(cor0):
    os.remove(cor0)

cor = folder_name + '/cor.csv'
if os.path.exists(cor):
    os.remove(cor)

F_M_file = folder_name + '/F_M.csv'
if os.path.exists(F_M_file):
    os.remove(F_M_file)

coor_f = folder_name + '/coor_femur.csv'
if os.path.exists(coor_f):
    os.remove(coor_f)

fileACLa = folder_name + '/ligaments/ACLa.csv'
fileACLp = folder_name + '/ligaments/ACLp.csv'
filePCLa = folder_name + '/ligaments/PCLa.csv'
filePCLp = folder_name + '/ligaments/PCLp.csv'
fileLCL = folder_name + '/ligaments/LCL.csv'
fileMCLa = folder_name + '/ligaments/MCLa.csv'
fileMCLo = folder_name + '/ligaments/MCLo.csv'
fileMCLd = folder_name + '/ligaments/MCLd.csv'
ligaments_files = [fileACLa, fileACLp, filePCLa, filePCLp, fileLCL, fileMCLa, fileMCLo, fileMCLd]

for j in range(len(ligaments_files)):
    if os.path.exists(ligaments_files[j]):
        os.remove(ligaments_files[j])

# set visualisation
pv.global_theme.show_edges = True

# coordinate of femur
coor_femur = pv.Box(bounds=(0, 1, 0, 1, 0, 1))
coor_femur_transform = np.array([[m.cos(m.radians(171.783943631)), m.cos(m.radians(92.495895645)),
                                m.cos(m.radians(97.822782154)), 1.766255221],
                                [m.cos(m.radians(87.588247867)), m.cos(m.radians(177.464638325)),
                                m.cos(m.radians(89.218470512)), 7.144509290],
                                [m.cos(m.radians(97.849429969)), m.cos(m.radians(89.554676924)),
                                m.cos(m.radians(7.862211319)), 0.022166126],
                                [0, 0, 0, 1]])

# import stl models
femur = pv.read('Models-reduce/Segmentation_Model_95_femur.stl')
cartilage = pv.read('Models-reduce/Segmentation_Model_84_lateral_tibial_cartilage.stl') + \
            pv.read('Models-reduce/Segmentation_Model_85_medial_tibial_cartilage.stl')
tibia = pv.read('Models/Segmentation_Model_96_tibia.stl')
femoral_cartilage = pv.read('Models-reduce/Segmentation_Model_83_femoral_cartilage.stl')
tibial_cartilage = pv.read('Models-reduce/Segmentation_Model_84_lateral_tibial_cartilage.stl') + \
                    pv.read('Models-reduce/Segmentation_Model_85_medial_tibial_cartilage.stl')
meniscus = pv.read('Models-reduce/Segmentation_Model_81_lateral_meniscus.stl') + \
            pv.read('Models-reduce/Segmentation_Model_82_medial_meniscus.stl')

flex = femur
flex_cartilage = femoral_cartilage
full_flex = flex + flex_cartilage
full_tibia = tibia + tibial_cartilage + meniscus

#  ligament
#  attachment
ACLaf0 = np.array(femur.points[28])
ACLpf0 = np.array(femur.points[1142])

PCLaf0 = np.array(femur.points[907])
PCLpf0 = np.array(femur.points[169])

LCLf0 = np.array(femur.points[552])

MCLaf0 = np.array(femur.points[507])
MCLof0 = np.array(femur.points[173])
MCLdf0 = np.array(femur.points[1080])

#  input flexion kinematics and external forces
dfmot = pd.read_csv('input_motion.csv', sep=';')
dfext = pd.read_csv('input_external_force.csv', sep=';')
motion = dfmot.to_numpy()
external_forces = dfext.to_numpy()
