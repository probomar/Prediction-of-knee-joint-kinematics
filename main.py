from input import *
import tools as t
import forces as f

if os.path.exists(folder_name):
    print('Are you sure you want to overwrite saved data? Y/n')
    Yn = input()
else:
    Yn = 'y'

if Yn == 'y':
    # creation coordinate system of femur
    coor_femur.transform(coor_femur_transform)
    _, _, _, _, _, _, _, _, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, _, _, _, _, _, _, _, _ = f.resultant_force(0)
    t.flex_plot2(0, 0, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, axis=None)
    # initialization step (fi = 0)
    t.initialization()
    # simulation
    for i in range(int(motion.shape[0] / step)):
        t.update_scene(i)

else:
    print('Rename "folder_name" in "input.py" please.')
