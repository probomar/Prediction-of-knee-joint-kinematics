from input import *
import tools as t

# creation coordinate system of femur
coor_femur.transform(coor_femur_transform)

# initialization step (fi = 0)
t.initialization()
# simulation
for i in range(int(motion.shape[0]/step)):
    t.update_scene(i)
