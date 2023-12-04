import matplotlib.pyplot as plt
import forces as f
from input import *
import time as tim


# first step without rotation (full extension)
def initialization():
    start = tim.time()

    global fii
    fii = motion[0, 1]

    global fi_step
    fi_step = fii

    print('i =', 0)
    print('fi =', fii)

    F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus = \
        rolling(0, fii)
    flex_plot(fii, 0, F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus)

    end = tim.time()
    time = end - start
    print('time =', time)
    x = np.array(position())
    print('position =', x, '\n\n')


# rotate + optimization of position
def update_scene(i):
    start = tim.time()

    global fii
    fii = motion[(i + 1)*step, 1]

    global fi_step
    fi_step = fii - motion[i*step, 1]

    p1_1 = np.array(flex.points[14])
    p2_1 = np.array(flex.points[25])

    print('i =', i+1)
    print('fi =', fii)

    F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus = \
        rolling(i + 1, fii)

    p1_2 = np.array(flex.points[14])
    p2_2 = np.array(flex.points[25])

    axis = actual_axis_of_rotation(p1_1, p2_1, p1_2, p2_2)

    flex_plot(fii, i + 1, F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus, axis)

    end = tim.time()
    time = end - start
    print('time =', time)
    x = np.array(position())
    print('position =', x, '\n', '\n')


# rotation around the folding edge
def rolling(i, fii):
    pm, pl = before_rotation()
    flex.rotate_vector(vector=pm - pl, angle=-fi_step, point=pm, inplace=True)
    flex_cartilage.rotate_vector(vector=pm - pl, angle=-fi_step, point=pm, inplace=True)
    coor_femur.rotate_vector(vector=pm - pl, angle=-fi_step, point=pm, inplace=True)
    F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus = f.force_equilibrium(i, fii)
    return F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus


# find contact point if it doesn't exist use last
def contact_point(last_center=np.array([0, 0, 0])):
    collision, ncol = (tibia + tibial_cartilage).collision((flex + flex_cartilage), generate_scalars=True)
    if ncol != 0:
        contact_volumes = (tibia + tibial_cartilage).boolean_intersection(flex + flex_cartilage)
        if contact_volumes.n_cells == 0:
            # centr = np.array([0, 0, 0])
            centr = last_center
        else:
            threshed = contact_volumes.threshold()
            bodies = threshed.split_bodies()
            volume = np.empty((0, 1), int)
            center = np.empty((0, 3), int)
            for i, body in enumerate(bodies):
                volume = np.append(volume, body.volume)
                center = np.append(center, [np.array(body.center)], axis=0)
            ind = np.argmax(volume)
            centr = center[ind]
    else:
        centr = last_center
    return centr


# return medial and lateral contact point for folding edge
def before_rotation():
    cor = open(cor0, 'a')
    a = 0

    while True:
        z = 0
        fiy = 0
        while True:
            collision, ncol = (tibia + tibial_cartilage).collision(flex + flex_cartilage)

            if ncol == 0:
                print('z')
                transform = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, -z_step], [0, 0, 0, 1]])
                flex.transform(transform)
                flex_cartilage.transform(transform)
                coor_femur.transform(transform)
                z -= z_step
            else:
                break

        bodies, points_medial, points_lateral = contact_volume()

        if (not points_medial.any()) and (a == 1):
            print('medial, z')
            transform = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, -z_step], [0, 0, 0, 1]])
            flex.transform(transform)
            flex_cartilage.transform(transform)
            coor_femur.transform(transform)
            z -= z_step
            a = 0
        elif (not points_medial.any()) and (a != 1):
            print('medial')
            flex.rotate_y(- fiy_step, inplace=True)
            flex_cartilage.rotate_y(- fiy_step, inplace=True)
            coor_femur.rotate_y(- fiy_step, inplace=True)
            fiy -= fiy_step
            a = -1
        elif (not points_lateral.any()) and (a == -1):
            print('lateral, z')
            transform = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, -z_step], [0, 0, 0, 1]])
            flex.transform(transform)
            flex_cartilage.transform(transform)
            coor_femur.transform(transform)
            z -= z_step
            a = 0
        elif (not points_lateral.any()) and (a != -1):
            print('lateral')
            flex.rotate_y(fiy_step, inplace=True)
            flex_cartilage.rotate_y(fiy_step, inplace=True)
            coor_femur.rotate_y(fiy_step, inplace=True)
            fiy += fiy_step
            a = 1
        else:
            break

    print('dz=', z, ', fiy=', fiy)

    cor.write(str('dz = '))
    cor.write(str(z))
    cor.write(str(', dfiy = '))
    cor.write(str(fiy))
    cor.write('\n')
    cor.close()

    pm = np.mean(points_medial, axis=0)
    pl = np.mean(points_lateral, axis=0)

    return pm, pl


# return contact volumes, medial and lateral contact points
def contact_volume():
    contact_volumes = (tibia + tibial_cartilage).boolean_intersection(flex + flex_cartilage)
    threshed = contact_volumes.threshold(0.001, invert=True)
    bodies = threshed.split_bodies()

    points_medial = np.empty((0, 3), int)
    points_lateral = np.empty((0, 3), int)

    for i in bodies.keys():
        point = np.array(bodies[i].center)
        if point[0] < 0:
            points_medial = np.append(points_medial, [np.array(bodies[i].center)], axis=0)
        else:
            points_lateral = np.append(points_lateral, [np.array(bodies[i].center)], axis=0)
    return contact_volumes, points_medial, points_lateral


# find actual axis of rotation
def actual_axis_of_rotation(p1_1, p2_1, p1_2, p2_2):
    p1 = np.mean([p1_1, p1_2], axis=0)
    p2 = np.mean([p2_1, p2_2], axis=0)

    n1 = p1_1 - p1_2
    n2 = p2_1 - p2_2

    d1 = -n1 @ p1
    d2 = -n2 @ p2

    a = (n2[0] * n1[2] - n1[0] * n2[2]) / (n1[1] * n2[2] - n2[1] * n1[2])
    b = (n1[2] * d2 - n2[2] * d1) / (n1[1] * n2[2] - n2[1] * n1[2])

    x1 = 30
    y1 = a * x1 + b
    z1 = ((- n1[0] - a * n1[1]) * x1 - b * n1[1] - d1) / n1[2]

    x2 = -30
    y2 = a * x2 + b
    z2 = ((- n1[0] - a * n1[1]) * x2 - b * n1[1] - d1) / n1[2]

    pA = np.array([x1, y1, z1])
    pB = np.array([x2, y2, z2])
    axis = [pA, pB]

    axis_of_rotation = pd.DataFrame(np.append(pA, pB).reshape(1, 6))
    axis_of_rotation.to_csv(cor, index=False, mode='a', header=False)

    return axis


# return position of femun from frame (x, y, z)
def position():
    flex_position = flex.center
    x0 = flex_position[0]
    y0 = flex_position[1]
    z0 = flex_position[2]
    return x0, y0, z0


# tool for plotting ligaments and its generated force
def plot_lig(p2, lig):
    F_lig = lig[0]
    lig1 = lig[1]
    lig2 = lig[2]
    lig_name = lig[3]
    lig_force_name = lig[4]

    p2.add_mesh(pv.Line(lig1, lig2), color='springgreen', line_width=2, name=lig_name)

    if np.linalg.norm(F_lig) != 0:
        F_lig_direction = F_lig / np.linalg.norm(F_lig)
        line = pv.Line(lig2, lig2 - F_lig)
        tip = pv.Cone(center=lig2 - F_lig_direction * 5, direction=F_lig_direction, height=10, radius=2)
        p2.add_mesh(line + tip, color='seagreen', line_width=5, name=lig_force_name)


# plot bones, cartilages, meniscus, ligaments and forces, use in flex_plot
def plot_flex_view(p2, F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus, axis=None):
    p2.add_mesh(flex, style='wireframe', color='linen')
    p2.add_mesh(tibia, style='wireframe', color='linen')
    p2.add_mesh(femoral_cartilage, style='wireframe', color='gold')
    p2.add_mesh(tibial_cartilage, style='wireframe', color='gold')
    p2.add_mesh(meniscus, style='wireframe', color='orange')
    line = pv.Line(axis[0], axis[1])
    p2.add_mesh(line, color='darkslateblue', line_width=5, name='axis')

    F_direction = -F / np.linalg.norm(F)
    line = pv.Line(soa, soa + F)
    tip = pv.Cone(center=soa - F_direction * 5, direction=F_direction, height=10, radius=2)
    p2.add_mesh(line + tip, color='darkred', line_width=5, name='force')

    if np.linalg.norm(N) != 0:
        N_direction = N / np.linalg.norm(N)
        line = pv.Line(soaN, soaN - N)
        tip = pv.Cone(center=soaN - N_direction * 5, direction=N_direction, height=10, radius=2)
        p2.add_mesh(line + tip, color='goldenrod', line_width=5, name='N')

    F_mus_direction = F_mus / np.linalg.norm(F_mus)
    line = pv.Line(soa_mus, soa_mus - F_mus)
    tip = pv.Cone(center=soa_mus - F_mus_direction * 5, direction=F_mus_direction, height=10, radius=2)
    p2.add_mesh(line + tip, color='coral', line_width=5, name='F_mus')

    plot_lig(p2, ACLa)
    plot_lig(p2, ACLp)
    plot_lig(p2, PCLa)
    plot_lig(p2, PCLp)
    plot_lig(p2, LCL)
    plot_lig(p2, MCLa)
    plot_lig(p2, MCLo)
    plot_lig(p2, MCLd)


def plot_flex_view2(p2, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, axis=None):
    p2.add_mesh(flex, style='wireframe', color='linen')
    p2.add_mesh(tibia, style='wireframe', color='linen')
    p2.add_mesh(femoral_cartilage, style='wireframe', color='gold')
    p2.add_mesh(tibial_cartilage, style='wireframe', color='gold')
    p2.add_mesh(meniscus, style='wireframe', color='orange')
    line = pv.Line(axis[0], axis[1])
    p2.add_mesh(line, color='darkslateblue', line_width=5, name='axis')

    plot_lig(p2, ACLa)
    plot_lig(p2, ACLp)
    plot_lig(p2, PCLa)
    plot_lig(p2, PCLp)
    plot_lig(p2, LCL)
    plot_lig(p2, MCLa)
    plot_lig(p2, MCLo)
    plot_lig(p2, MCLd)


# plot front and side view of knee with its soft issues and forces
def flex_plot(fi_act, time, F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus, axis=None):

    if axis is None:
        axis = [np.array([0, 0, 0]), np.array([0, 0, 0])]

    if time < 10:
        flex_pic = folder_name + '/flex/flex000' + str(time) + '.png'
    elif time < 100:
        flex_pic = folder_name + '/flex/flex00' + str(time) + '.png'
    elif time < 1000:
        flex_pic = folder_name + '/flex/flex0' + str(time) + '.png'
    else:
        flex_pic = folder_name + '/flex/flex' + str(time) + '.png'

    text1 = str(round(motion[(time+1)*step, 2], 1)) + ' % cyklu chůze'
    print(text1)
    text2 = 'fi = ' + str(round(-fi_act, 2)) + '°'

    p2 = pv.Plotter(off_screen=True, shape=(1, 2))
    p2.background_color = 'white'

    p2.subplot(0, 0)
    p2.add_text(text2, font_size=20, color='black')
    p2.add_text(text1, font_size=20, position='lower_left', color='black')
    plot_flex_view(p2, F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus, axis)
    p2.camera.position = (0, -500, 10)
    p2.camera.focal_point = (0, 100, 10)

    p2.subplot(0, 1)
    plot_flex_view(p2, F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus, axis)
    p2.camera.position = (-500, 0, 10)
    p2.camera.focal_point = (100, 0, 10)

    p2.show(screenshot=flex_pic, window_size=(1200, 1200), title=(text1+text2))


def flex_plot2(fi_act, time, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, axis=None):
    if axis is None:
        axis = [np.array([0, 0, 0]), np.array([0, 0, 0])]

    if time < 10:
        flex_pic = folder_name + '/flex/flex000' + str(time) + '.png'
    elif time < 100:
        flex_pic = folder_name + '/flex/flex00' + str(time) + '.png'
    elif time < 1000:
        flex_pic = folder_name + '/flex/flex0' + str(time) + '.png'
    else:
        flex_pic = folder_name + '/flex/flex' + str(time) + '.png'

    text1 = str(round(motion[(time+1)*step, 2], 1)) + ' % cyklu chůze'
    print(text1)
    text2 = 'fi = ' + str(round(-fi_act, 2)) + '°'

    p2 = pv.Plotter(off_screen=True, shape=(1, 2))
    p2.background_color = 'white'

    p2.subplot(0, 0)
    p2.add_text(text2, font_size=20, color='black')
    p2.add_text(text1, font_size=20, position='lower_left', color='black')
    plot_flex_view2(p2, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, axis)
    p2.camera.position = (0, -500, 10)
    p2.camera.focal_point = (0, 100, 10)

    p2.subplot(0, 1)
    plot_flex_view2(p2, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, axis)
    p2.camera.position = (-500, 0, 10)
    p2.camera.focal_point = (100, 0, 10)

    p2.show(screenshot=flex_pic, window_size=(1200, 1200), title=(text1+text2))


# plot femur with color map
def color_plot(fi_act, time, data, name):
    if time < 10:
        file2 = folder_name + '/' + name + '/data/' + name + '000' + str(time) + '.csv'
        picture = folder_name + '/' + name + '/figure/' + name + '000' + str(time) + '.png'
    elif time < 100:
        file2 = folder_name + '/' + name + '/data/' + name + '00' + str(time) + '.csv'
        picture = folder_name + '/' + name + '/figure/' + name + '00' + str(time) + '.png'
    elif time < 1000:
        file2 = folder_name + '/' + name + '/data/' + name + '0' + str(time) + '.csv'
        picture = folder_name + '/' + name + '/figure/' + name + '0' + str(time) + '.png'
    else:
        file2 = folder_name + '/' + name + '/data/' + name + str(time) + '.csv'
        picture = folder_name + '/' + name + '/figure/' + name + str(time) + '.png'

    text1 = str(round(motion[(time + 1) * step, 2], 1)) + ' % cyklu chůze'
    text2 = 'fi = ' + str(round(-fi_act, 2)) + '°'

    if os.path.exists(file2):
        os.remove(file2)

    df_stress = pd.DataFrame(data)
    df_stress.to_csv(file2, index=False, mode='a', header=False)

    mesh1 = flex.extract_cells(range(flex.n_cells))
    mesh1.cell_data['Normal stress [MPa]'] = data

    p1 = pv.Plotter(off_screen=True)
    p1.background_color = 'white'

    p1.add_text(text1, font_size=10, position='lower_left', color='black')
    p1.add_text(text2, font_size=10, color='black')
    sarg = dict(color='black')
    p1.add_mesh(mesh1, scalars='Normal stress [MPa]', clim=[0.0000000001, 10], below_color='grey', above_color='red',
                reset_camera=True, scalar_bar_args=sarg)  # , style='wireframe', color='linen')
    p1.view_yx()
    p1.set_viewup([0, -1, 0])
    p1.show(screenshot=picture, window_size=(960, 720))


def angle(vector_1, vector_2):
    unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
    unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    ang = np.arccos(dot_product) * 180 / m.pi
    if m.isnan(ang):
        ang = 0
    return ang


def plane_distance(xyz0, xyz1, n):
    d = - n[0] * xyz0[0] - n[1] * xyz0[1] - n[2] * xyz0[2]
    e = abs((n[0] * xyz1[0] + n[1] * xyz1[1] + n[2] * xyz1[2] + d))
    f = (m.sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]))
    distance = e / f
    return distance


def vector(coor):
    x = coor.points[1] - coor.points[0]
    y = coor.points[2] - coor.points[0]
    z = coor.points[4] - coor.points[0]
    return x, y, z


def cos(vec1, vec2):
    cosine = m.cos(m.radians(angle(vec1, vec2)))
    return cosine


def transform_matrix(coor1, coor2):
    x1, y1, z1 = vector(coor1)
    x2, y2, z2 = vector(coor2)

    r = coor2.points[0] - coor1.points[0]
    r_length = np.linalg.norm(r)

    x0 = r_length * m.cos(m.radians(angle(x1, r)))
    y0 = r_length * m.cos(m.radians(angle(y1, r)))
    z0 = r_length * m.cos(m.radians(angle(z1, r)))

    xx = cos(x1, x2)
    yy = cos(y1, y2)
    zz = cos(z1, z2)
    xy = cos(x1, y2)
    xz = cos(x1, z2)
    yx = cos(y1, x2)
    yz = cos(y1, z2)
    zx = cos(z1, x2)
    zy = cos(z1, y2)

    transform = np.array([[xx, xy, xz, x0],
                          [yx, yy, yz, y0],
                          [zx, zy, zz, z0],
                          [0, 0, 0, 1]])
    return transform


def motions(coor, transform01):
    R1 = np.empty([0])
    R2 = np.empty([0])
    R3 = np.empty([0])
    T1 = np.empty([0])
    T2 = np.empty([0])
    T3 = np.empty([0])

    for i in range(coor.shape[0]):
        x0 = np.array([1, 0, 0])
        y0 = np.array([0, 1, 0])
        z0 = np.array([0, 0, 1])

        xyz2 = coor[i, 0:3]
        point_x = coor[i, 3:6]
        point_y = coor[i, 6:9]
        point_z = coor[i, 9:12]

        x2 = point_x - xyz2
        y2 = point_y - xyz2
        z2 = point_z - xyz2

        xx = cos(x0, x2)
        yy = cos(y0, y2)
        zz = cos(z0, z2)
        xy = cos(x0, y2)
        xz = cos(x0, z2)
        yx = cos(y0, x2)
        yz = cos(y0, z2)
        zx = cos(z0, x2)
        zy = cos(z0, y2)

        transform02 = np.array([[xx, xy, xz, xyz2[0]],
                                [yx, yy, yz, xyz2[1]],
                                [zx, zy, zz, xyz2[2]],
                                [0, 0, 0, 1]])

        transform = np.dot(np.linalg.inv(transform02), transform01)
        rotation_matrix = transform[0:3, 0:3]
        r = R.from_matrix(rotation_matrix)
        R1i, R2i, R3i = r.as_euler('xyz', degrees=True)
        R1 = np.append(R1, R1i)
        R2 = np.append(R2, R2i)
        R3 = np.append(R3, R3i)

        A = transform[0, 3]
        B = transform[1, 3]
        C = transform[2, 3]

        cR1 = m.cos(m.radians(R1[i]))
        sR1 = m.sin(m.radians(R1[i]))
        cR2 = m.cos(m.radians(R2[i]))
        sR2 = m.sin(m.radians(R2[i]))

        T3 = np.append(T3, (C / sR1 - B / cR1) / cR2 / (cR1 / sR1 + sR1 / cR1))
        T1 = np.append(T1, A - T3[i] * sR2)
        T2 = np.append(T2, B / cR1 + T3[i] * sR1 * cR2 / cR1)

    return R1, R2, R3, T1, T2, T3


def stress_time(name, fold_name):
    df = pd.read_csv(fold_name + '/cor.csv', sep=',', header=None)
    cor_lig_lin = df.to_numpy()
    num = cor_lig_lin.shape[0]

    normal_vector_femur = np.array(flex.cell_normals)
    points_femur = np.array(femur.cell_centers().points)
    area_femur = flex.compute_cell_sizes(length=False, volume=False).cell_data["Area"]
    area_femur = area_femur.reshape([area_femur.shape[0], 1])

    N_lenght_m = np.empty([0])
    N_lenght_l = np.empty([0])
    stress_m = np.empty([0])
    stress_l = np.empty([0])

    for i in range(num+1):
        print(i)
        if i < 10:
            file_name = '000' + str(i) + '.csv'
        elif i < 100:
            file_name = '00' + str(i) + '.csv'
        elif i < 1000:
            file_name = '0' + str(i) + '.csv'
        else:
            file_name = str(i) + '.csv'

        dfstress_meniscus = pd.read_csv(fold_name + name + file_name, header=None)
        stress = dfstress_meniscus.to_numpy()

        N_m = np.array([0, 0, 0])
        N_l = np.array([0, 0, 0])
        soaN_m = np.array([0, 0, 0])
        soaN_l = np.array([0, 0, 0])
        s_m = np.array([0])
        s_l = np.array([0])

        for j in range(stress.shape[0]):
            if points_femur[j, 0] > 0:
                Nj_m = np.array(stress[j] * area_femur[j] * normal_vector_femur[j])
                N_m, _, soaN_m = f.result_of_forces_and_moments(N_m, soaN_m, Nj_m, points_femur[j])
                if s_m < stress[j]:
                    s_m = stress[j]
            else:
                Nj_l = np.array(stress[j] * area_femur[j] * normal_vector_femur[j])
                N_l, _, soaN_l = f.result_of_forces_and_moments(N_l, soaN_l, Nj_l, points_femur[j])
                if s_l < stress[j]:
                    s_l = stress[j]

        N_lenght_m = np.append(N_lenght_m, np.linalg.norm(N_m))
        N_lenght_l = np.append(N_lenght_l, np.linalg.norm(N_l))
        stress_m = np.append(stress_m, np.linalg.norm(s_m))
        stress_l = np.append(stress_l, np.linalg.norm(s_l))

    return N_lenght_m, N_lenght_l, stress_m, stress_l


def residual_time(force):
    force = np.split(force, np.where(np.diff(force[:, 0]))[0] + 1)

    F = np.empty([0])
    M = np.empty([0])

    for i in range(len(force)):
        a = force[i]
        Fi = a[a.shape[0] - 1, 3]
        Mi = a[a.shape[0] - 1, 4]
        F = np.append(F, Fi)
        M = np.append(M, Mi)

    return F, M


def plot_kinematics(mx, my, myerr, xlin, ylin, xnonlin, ynonlin, title, ylabel, name, ymin, ymax):
    plt.plot(mx, my, label='experimentální data')
    plt.fill_between(mx, my - myerr, my + myerr, alpha=0.2)
    plt.plot(xlin, ylin, '^', label='simulace-lin. vazy')
    plt.plot(xnonlin, ynonlin, 'x', label='simulace-nonlin. vazy')
    plt.xlabel('% cyklu chůze')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.ylim(ymin, ymax)
    plt.xlim(10, 50)
    plt.savefig(name, dpi=300)
    plt.show()

