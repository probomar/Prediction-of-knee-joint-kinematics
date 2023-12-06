from input import *
import tools as t


# return normal force from cartilages and meniscus
def normal_force():
    normal_vector_femur = np.array(flex.cell_normals)
    points_femur = np.array(femur.cell_centers().points)
    area_femur = flex.compute_cell_sizes(length=False, volume=False).cell_data["Area"]
    area_femur = area_femur.reshape([area_femur.shape[0], 1])

    N_cartilage = np.empty([0, 3])
    stress_cartilage = np.empty([0, 1])
    distance_cartilage0 = np.empty([0, 1])
    cartilage_femoral = np.empty([0, 1])
    cartilage_tibial = np.empty([0, 1])
    N_bone = np.empty([0, 3])
    stress_bone = np.empty([0, 1])
    N_meniscus = np.empty([0, 3])
    stress_meniscus = np.empty([0, 1])
    distance_meniscus0 = np.empty([0, 1])

    for n in range(points_femur.shape[0]):
        point_femur = points_femur[n]
        pA = point_femur - 10000000 * normal_vector_femur[n]
        pB = point_femur + 10000000 * normal_vector_femur[n]
        points_tibia = tibia.ray_trace(pA, pB)[0]

        # measurement of the thickness of the femoral cartilage
        points_femoral_cartilage = flex_cartilage.ray_trace(pA, pB)[0]
        if np.shape(points_femoral_cartilage)[0] >= 2:
            d = []
            for k in range(np.shape(points_femoral_cartilage)[0]):
                d.append(np.linalg.norm(points_femoral_cartilage[k, :] - point_femur))
            ind = [d.index(x) for x in sorted(d)]
            if d[ind[0]] > 5:
                femoral_cartilage_distance = 0
            else:
                fc1 = points_femoral_cartilage[ind[0], :]
                fc2 = points_femoral_cartilage[ind[1], :]
                femoral_cartilage_distance = np.linalg.norm(fc1 - fc2)
        else:
            femoral_cartilage_distance = 0

        cartilage_femoral = np.append(cartilage_femoral, [[femoral_cartilage_distance]], axis=0)

        # meassurement of meniscus
        points_meniscus = meniscus.ray_trace(pA, pB)[0]
        if np.shape(points_meniscus)[0] >= 2:
            # measurement of the thickness of meniscus
            d = []
            for k in range(np.shape(points_meniscus)[0]):
                d.append(np.linalg.norm(points_meniscus[k, :] - point_femur))
            ind = [d.index(x) for x in sorted(d)]

            m1 = points_meniscus[ind[0], :]
            m2 = points_meniscus[ind[1], :]
            meniscus_distance = np.linalg.norm(m1 - m2)
        else:
            meniscus_distance = 0

        distance_meniscus0 = np.append(distance_meniscus0, [[meniscus_distance]])

        # measurement of the thickness of the tibial cartilage
        points_tibial_cartilage = tibial_cartilage.ray_trace(pA, pB)[0]
        if np.shape(points_tibial_cartilage)[0] >= 2:
            if not points_tibia.any():
                tibial_cartilage_distance = 0
                stress_bone = np.append(stress_bone, [[0]], axis=0)
                N_bone = np.append(N_bone, [[0, 0, 0]], axis=0)


            else:
                d = []
                for k in range(np.shape(points_tibia)[0]):
                    d.append(np.linalg.norm(points_tibia[k, :] - point_femur))
                ind = [d.index(x) for x in sorted(d)]
                point_tibia = points_tibia[ind[0]]

                # bones penetration?
                ft = np.linalg.norm(point_femur - points_tibia[ind[1]])
                tt = np.linalg.norm(point_tibia - points_tibia[ind[1]])

                if ft < tt:
                    print('bone')
                    stress_bone = np.append(stress_bone, [[k_bone * (tt - ft)]], axis=0)
                    N_bone = np.append(N_bone, [k_bone * (tt - ft) * area_femur[n] * normal_vector_femur[n]], axis=0)
                else:
                    stress_bone = np.append(stress_bone, [[0]], axis=0)
                    N_bone = np.append(N_bone, [[0, 0, 0]], axis=0)

                # measurement of the thickness of the tibial cartilage
                d = []
                for k in range(np.shape(points_tibial_cartilage)[0]):
                    d.append(np.linalg.norm(points_tibial_cartilage[k, :] - point_tibia))
                ind = [d.index(x) for x in sorted(d)]
                if d[ind[0]] > 5:
                    tibial_cartilage_distance = 0
                else:
                    tc1 = points_tibial_cartilage[ind[0], :]
                    tc2 = points_tibial_cartilage[ind[1], :]
                    tibial_cartilage_distance = np.linalg.norm(tc1 - tc2)
        else:
            tibial_cartilage_distance = 0
            stress_bone = np.append(stress_bone, [[0]], axis=0)
            N_bone = np.append(N_bone, [[0, 0, 0]], axis=0)

        cartilage_tibial = np.append(cartilage_tibial, [[tibial_cartilage_distance]], axis=0)

        cartilage_distance = femoral_cartilage_distance + tibial_cartilage_distance
        distance_cartilage0 = np.append(distance_cartilage0, [[cartilage_distance]], axis=0)

        # normal force from meniscus
        if not points_meniscus.any():
            stress_meniscus = np.append(stress_meniscus, [[0]], axis=0)
            N_meniscus = np.append(N_meniscus, [[0, 0, 0]], axis=0)
        else:
            min_distance = np.linalg.norm(point_femur - m2)
            if min_distance < meniscus_distance:
                stress_meniscus = np.append(stress_meniscus, [
                    [- (E_meniscus * (1 - ny_meniscus) / (1 + ny_meniscus) / (1 - 2 * ny_meniscus))
                     * (min_distance - meniscus_distance)
                     / meniscus_distance]], axis=0)
                N_meniscus = np.append(N_meniscus, [
                    (E_meniscus * (1 - ny_meniscus) / (1 + ny_meniscus) / (1 - 2 * ny_meniscus))
                    * (min_distance - meniscus_distance) / meniscus_distance * area_femur[n] *
                    normal_vector_femur[n]], axis=0)
            else:
                stress_meniscus = np.append(stress_meniscus, [[0]], axis=0)
                N_meniscus = np.append(N_meniscus, [[0, 0, 0]], axis=0)

        # normal force from cartilages
        distance = np.empty([0])
        if not points_tibia.any():
            stress_cartilage = np.append(stress_cartilage, [[0]], axis=0)
            N_cartilage = np.append(N_cartilage, [[0, 0, 0]], axis=0)
        else:
            for j in range(points_tibia.shape[0]):
                distance = np.append(distance, np.linalg.norm(points_femur[n, :] - points_tibia[j, :]))

            min_distance = np.min(distance)

            if min_distance < cartilage_distance:
                stress_cartilage = np.append(stress_cartilage, [
                    [- (E_cartilage * (1 - ny_cartilage) / (1 + ny_cartilage) / (1 - 2 * ny_cartilage))
                     * (min_distance - cartilage_distance)
                     / cartilage_distance]], axis=0)
                N_cartilage = np.append(N_cartilage, [
                    (E_cartilage * (1 - ny_cartilage) / (1 + ny_cartilage) / (1 - 2 * ny_cartilage))
                    * (min_distance - cartilage_distance)
                    / cartilage_distance * area_femur[n]
                    * normal_vector_femur[n]], axis=0)
            else:
                stress_cartilage = np.append(stress_cartilage, [[0]], axis=0)
                N_cartilage = np.append(N_cartilage, [[0, 0, 0]], axis=0)

    N_forces_cartilage = N_cartilage.reshape([N_cartilage.shape[0], 3])
    N_forces_meniscus = N_meniscus.reshape([N_meniscus.shape[0], 3])
    N_force_bone = N_bone.reshape([N_bone.shape[0], 3])

    N_cartilage = N_forces_cartilage[0]
    N_meniscus = N_forces_meniscus[0]
    N_bone = N_force_bone[0]
    soaN_cartilage = points_femur[0]
    soaN_meniscus = points_femur[0]
    soaN_bone = points_femur[0]

    # resultant forces from cartilage, meniscus, bone
    for j in range(N_forces_cartilage.shape[0] - 1):
        N_cartilage, MN_cartilage, soaN_cartilage = result_of_forces_and_moments(N_cartilage, soaN_cartilage,
                                                                                 N_forces_cartilage[j + 1],
                                                                                 points_femur[j + 1])
        N_meniscus, MN_meniscus, soaN_meniscus = result_of_forces_and_moments(N_meniscus, soaN_meniscus,
                                                                              N_forces_meniscus[j + 1],
                                                                              points_femur[j + 1])

        N_bone, MN_bone, soaN_bone = result_of_forces_and_moments(N_bone, soaN_bone, N_force_bone[j + 1],
                                                                  points_femur[j + 1])
    # resultant force
    N, MN, soaN = result_of_forces_and_moments(N_cartilage, soaN_cartilage, N_meniscus, soaN_meniscus, N_bone,
                                               soaN_bone)

    stress = stress_cartilage + stress_meniscus + stress_bone
    return N, soaN, stress, stress_cartilage, distance_cartilage0, cartilage_femoral, cartilage_tibial, stress_bone, \
        stress_meniscus, distance_meniscus0


def moment_of_force(force, site_of_act, point=np.array([0, 0, 0])):
    if force.all() == 0:
        moment = np.array([0, 0, 0])
    else:
        param = (point - site_of_act @ force) / (force @ force)
        arm = site_of_act + force * param - point
        moment = np.array([- force[1] * arm[2] + force[2] * arm[1],
                           - force[2] * arm[0] + force[0] * arm[2],
                           - force[0] * arm[1] + force[1] * arm[0]])

    return moment


def site_of_action(force, moment):
    F = [[0, force[2], - force[1]], [- force[2], 0, force[0]], [force[1], - force[0], 0]]
    Ff = np.linalg.pinv(F)
    Mm = np.reshape(moment, (3, 1))
    X = np.dot(Ff, Mm)
    soa = np.ravel(X)
    return soa


def result_of_forces_and_moments(force1, site_of_action1,
                                 force2=np.array([0, 0, 0]), site_of_action2=np.array([0, 0, 0]),
                                 force3=np.array([0, 0, 0]), site_of_action3=np.array([0, 0, 0]),
                                 force4=np.array([0, 0, 0]), site_of_action4=np.array([0, 0, 0]),
                                 force5=np.array([0, 0, 0]), site_of_action5=np.array([0, 0, 0]),
                                 force6=np.array([0, 0, 0]), site_of_action6=np.array([0, 0, 0]),
                                 force7=np.array([0, 0, 0]), site_of_action7=np.array([0, 0, 0]),
                                 force8=np.array([0, 0, 0]), site_of_action8=np.array([0, 0, 0]),
                                 force9=np.array([0, 0, 0]), site_of_action9=np.array([0, 0, 0]),
                                 force10=np.array([0, 0, 0]), site_of_action10=np.array([0, 0, 0])):
    force = force1 + force2 + force3 + force4 + force5 + force6 + force7 + force8 + force9 + force10

    moment1 = moment_of_force(force1, site_of_action1)
    moment2 = moment_of_force(force2, site_of_action2)
    moment3 = moment_of_force(force3, site_of_action3)
    moment4 = moment_of_force(force4, site_of_action4)
    moment5 = moment_of_force(force5, site_of_action5)
    moment6 = moment_of_force(force6, site_of_action6)
    moment7 = moment_of_force(force7, site_of_action7)
    moment8 = moment_of_force(force8, site_of_action8)
    moment9 = moment_of_force(force9, site_of_action9)
    moment10 = moment_of_force(force10, site_of_action10)

    moment = moment1 + moment2 + moment3 + moment4 + moment5 + moment6 + moment7 + moment8 + moment9 + moment10

    soa = site_of_action(force, moment)

    return force, moment, soa


def ligament_force(lig0, lig1, lig2, k1lig, k2lig, name1, name2):
    lig = lig2 - lig1
    lig_length = np.linalg.norm(lig)
    lig_length0 = lig0
    dlig_length = lig_length - lig_length0
    nlig = lig / lig_length

    if ligament == 'linear':
        if lig0 < lig_length:
            F_lig = - k1lig * dlig_length * nlig
        else:
            F_lig = np.array([0, 0, 0])

    else:
        eps = dlig_length / lig0
        if (eps >= 0) and (eps <= 2 * epsL):
            F_lig = - k2lig * (dlig_length ** 2) * nlig
        elif eps > 2 * epsL:
            F_lig = - k1lig * (lig_length - (1 + epsL) * lig_length0) * nlig
        else:
            F_lig = np.array([0, 0, 0])
    ligament_value = [F_lig, lig1, lig2, name1, name2, lig0, lig_length]
    return ligament_value


def external_force(i):
    force = - 1 * np.array([-external_forces[i, 3], -external_forces[i, 1], external_forces[i, 2]])
    moment = - 1 * np.array([-external_forces[i, 6], -external_forces[i, 4], external_forces[i, 5]]) * 1000
    soa = site_of_action(force, moment)
    return force, moment, soa


def best_optimization(F_best, M_best, position_best, F, M):
    if F < F_best and M < M_best:
        F_best = F
        M_best = M
        position_best = coor_femur
    return F_best, M_best, position_best


def force_equilibrium(i, fii):
    n = 0
    step_F = step_F0
    step_M = step_M0

    F, M, soa, center, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = resultant_force(i)

    F_length = np.linalg.norm(F)
    F_best = F_length
    F1 = F_length

    M_length = np.linalg.norm(M)
    M_best = M_length
    M1 = M_length
    position_best = coor_femur

    # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
    #       (' ', ' ', step_F, 'F =', F_length, 'N', 'M =', M_length/1000, 'Nm'))

    # F_M = pd.DataFrame(np.append(i, [fii, n, F_length, M_length]).reshape(1, 5))
    # F_M.to_csv(F_M_file, index=False, mode='a', header=False)

    t.print_forces(' ', ' ', step_F, F_length, M_length)
    t.save_forces(i, fii, n, F_length, M_length)

    while True:

        if (F_length < Fr) and (M_length < Mr):
            break

        n += 1

        if M_length > Mr/5:
            M_length_old = M_length

            rotate_angle = M_length * step_M * 180 / m.pi
            flex.rotate_vector(vector=M, angle=rotate_angle, point=center, inplace=True)
            flex_cartilage.rotate_vector(vector=M, angle=rotate_angle, point=center, inplace=True)
            coor_femur.rotate_vector(vector=M, angle=rotate_angle, point=center, inplace=True)

            F, M, _, center, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = resultant_force(i, center)
            F_length = np.linalg.norm(F)
            M_length = np.linalg.norm(M)
            F2 = F_length
            M2 = M_length

            if (F2 > 2*F1) or (M2 > 2*M1):
                print('rotate (F2 > 2*F1) or (M2 > 2*M1)')
                flex.rotate_vector(vector=M, angle=-rotate_angle/2, point=center, inplace=True)
                flex_cartilage.rotate_vector(vector=M, angle=-rotate_angle/2, point=center, inplace=True)
                coor_femur.rotate_vector(vector=M, angle=-rotate_angle/2, point=center, inplace=True)
                F, M, _, center, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = resultant_force(i, center)
                F_length = np.linalg.norm(F)
                M_length = np.linalg.norm(M)
                F1 = F_length
                M1 = M_length

            # F_M = pd.DataFrame(np.append(i, [fii, n, F_length, M_length]).reshape(1, 5))
            # F_M.to_csv(F_M_file, index=False, mode='a', header=False)
            t.save_forces(i, fii, n, F_length, M_length)

            if (F_length < Fr) and (M_length < Mr):
                # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
                #       (n, ' ', ' ', 'F =', F_length, 'N', 'M =', M_length/1000, 'Nm'))
                t.print_forces(n, ' ', ' ', F_length, M_length)
                break
            else:
                F_best, M_best, position_best = best_optimization(F_best, M_best, position_best, F_length, M_length)

            if M_length_old > M_length:
                step_M *= alpha
                step_M = round(step_M, 15)
                # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
                #       ('M', 'alpha', step_M, 'F =', F_length, 'N', 'M =', M_length/1000, 'Nm'))
                t.print_forces('M'+str(n), 'alpha', step_M, F_length, M_length)

            elif M_length_old < M_length:
                step_M *= beta * 1.5
                step_M = round(step_M, 15)
                # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
                #       ('M', 'beta', step_M, 'F =', F_length, 'N', 'M =', M_length/1000, 'Nm'))
                t.print_forces('M'+str(n), 'beta', step_M, F_length, M_length)

            else:
                # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
                #       ('M', ' ', step_M, 'F =', F_length, 'N', 'M =', M_length/1000, 'Nm'))
                t.print_forces('M'+str(n), ' ', step_M, F_length, M_length)
        else:
            print('skip M')


        if F_length > Fr/5:
            F_transform = np.array([[1, 0, 0, F[0] * step_F], [0, 1, 0, F[1] * step_F],
                                [0, 0, 1, F[2] * step_F], [0, 0, 0, 1]])
            flex.transform(F_transform, inplace=True)
            flex_cartilage.transform(F_transform, inplace=True)
            coor_femur.transform(F_transform, inplace=True)

            F_new, M_new, _, center, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = resultant_force(i, center)

            F_length_new = np.linalg.norm(F_new)
            M_length_new = np.linalg.norm(M_new)
            F2 = F_length_new
            M2 = F_length_new

            if (F2 > 2*F1) or (M2 > 2*M1):
                print('transform (F2 > 2*F1) or (M2 > 2*M1)')
                F_transform = np.array([[1, 0, 0, -F[0] * step_F/2], [0, 1, 0, -F[1] * step_F/2],
                                        [0, 0, 1, -[2] * step_F]/2, [0, 0, 0, 1]])
                flex.transform(F_transform, inplace=True)
                flex_cartilage.transform(F_transform, inplace=True)
                coor_femur.transform(F_transform, inplace=True)
                F, M, _, center, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = resultant_force(i, center)
                F_length = np.linalg.norm(F)
                M_length = np.linalg.norm(M)
                F1 = F_length
                M1 = M_length

            # F_M = pd.DataFrame(np.append(i, [fii, n, F_length_new, M_length_new]).reshape(1, 5))
            # F_M.to_csv(F_M_file, index=False, mode='a', header=False)
            t.save_forces(i, fii, n, F_length_new, M_length_new)

            if (F_length_new < Fr) and (M_length_new < Mr):
                # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
                #       (n, ' ', ' ', 'F =', F_length_new, 'N', 'M =', M_length_new/1000, 'Nm'))
                t.print_forces(n, ' ', ' ', F_length_new, M_length_new)
                break
            else:
                F_best, M_best, position_best = best_optimization(F_best, M_best, position_best, F_length, M_length)

            if 1 * F_length > F_length_new:
                step_F *= alpha
                step_F = round(step_F, 10)
                # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
                #       (n, 'alpha', step_F, 'F =', F_length_new, 'N', 'M =', M_length_new/1000, 'Nm'))
                t.print_forces(n, 'alpha', step_F, F_length_new, M_length_new)

            elif 1 * F_length < F_length_new:
                step_F *= beta
                step_F = round(step_F, 10)
                # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
                #       (n, 'beta', step_F, 'F =', F_length_new, 'N', 'M =', M_length_new/1000, 'Nm'))
                t.print_forces(n, 'beta', step_F, F_length_new, M_length_new)

            else:
                # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
                #       (n, ' ', step_F, 'F =', F_length_new, 'N', 'M =', M_length_new/1000, 'Nm'))
                t.print_forces(n, '', step_F, F_length_new, M_length_new)
        else:
            print('skip_F')

        F, F_length, M, M_length = F_new, F_length_new, M_new, M_length_new

        if n == n_max:
            # print("{:<3} {:<7} {:<15} {:<4} {:<20} {:<5} {:<4} {:<20} {:<5}".format
            #       (n, ' ', ' ', 'F =', F_length_new, 'N', 'M =', M_length_new / 1000, 'Nm'))
            t.print_forces(n, ' ', ' ', F_length_new, M_length_new)
            # transform back
            coor = coor_femur
            back_transform = t.transform_matrix(coor, position_best)
            flex.transform(back_transform, inplace=True)
            flex_cartilage.transform(back_transform, inplace=True)
            coor_femur.transform(back_transform, inplace=True)
            t.print_forces(' ', ' ', ' ', F_best, M_best)
            F, M, _, center, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = resultant_force(i, center)
            F_length = np.linalg.norm(F)
            M_length = np.linalg.norm(M)
            t.print_forces(' ', ' ', ' ', F_length, M_length)
            break

    F, M, soa, center, N, soaN, stress, stress_cartilage, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, \
        soa_mus, distance_cartilage0, cartilage_femoral, cartilage_tibial, stress_bone, stress_meniscus, \
        distance_meniscus0 = resultant_force(i, center)

    t.color_plot(fii, i, distance_cartilage0, 'distance_cartilage0')
    t.color_plot(fii, i, cartilage_femoral, 'cartilage_femoral')
    t.color_plot(fii, i, cartilage_tibial, 'cartilage_tibial')
    t.color_plot(fii, i, stress, 'stress')
    t.color_plot(fii, i, stress_cartilage, 'stress_cartilage')
    t.color_plot(fii, i, stress_meniscus, 'stress_meniscus')
    t.color_plot(fii, i, stress_bone, 'stress_bone')

    t.color_plot(fii, i, distance_meniscus0, 'distance_meniscus0')

    ligaments = [ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd]

    for j in range(len(ligaments)):
        F_lig = ligaments[j][0]
        lig0 = ligaments[j][5]
        lig_length = ligaments[j][6]
        lig1 = ligaments[j][1]
        lig2 = ligaments[j][2]
        np_print_lig = np.array([np.linalg.norm(F_lig), F_lig[0], F_lig[1], F_lig[2], lig0, lig_length, lig1[0],
                                 lig1[1], lig1[2], lig2[0], lig2[1], lig2[2]]).reshape(1, 12)
        print_lig = pd.DataFrame(np_print_lig)
        print_lig.to_csv(ligaments_files[j], index=False, mode='a', header=False)

    coor0 = coor_femur.points[0]
    x = coor_femur.points[1]
    y = coor_femur.points[2]
    z = coor_femur.points[4]
    print_coor = pd.DataFrame(np.append(coor0, [x, y, z]).reshape(1, 12))
    print_coor.to_csv(coor_f, index=False, mode='a', header=False)

    return F, soa, N, soaN, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, MCLd, F_mus, soa_mus


def resultant_force(i, center=np.array([0, 0, 0])):
    ACLaf = np.array(flex.points[28])
    ACLat = np.array(tibia.points[538])

    ACLpf = np.array(flex.points[1142])
    ACLpt = np.array(tibia.points[1001])

    PCLaf = np.array(flex.points[907])
    PCLat = np.array(tibia.points[1000])

    PCLpf = np.array(flex.points[169])
    PCLpt = np.array(tibia.points[189])

    LCLf = np.array(flex.points[552])
    LCLt = np.array([55, 20, -40])

    MCLaf = np.array(flex.points[507])
    MCLat = np.array(tibia.points[1272])

    MCLof = np.array(flex.points[173])
    MCLot = np.array(tibia.points[1206])

    MCLdf = np.array(flex.points[1080])
    MCLdt = np.array(tibia.points[718])

    ACLa0 = np.linalg.norm(ACLat - ACLaf0) / 1.00
    ACLp0 = np.linalg.norm(ACLpt - ACLpf0) / 1.051

    PCLa0 = np.linalg.norm(PCLat - PCLaf0) / 1.004
    PCLp0 = np.linalg.norm(PCLpt - PCLpf0) / 1.05

    LCL0 = np.linalg.norm(LCLt - LCLf0) / 1.05

    MCLa0 = np.linalg.norm(MCLat - MCLaf0) / 0.94
    MCLo0 = np.linalg.norm(MCLot - MCLof0) / 1.031
    MCLd0 = np.linalg.norm(MCLdt - MCLdf0) / 1.049

    ACLa = ligament_force(ACLa0, ACLat, ACLaf, k1ACLa, k2ACLa, 'ACLa', 'ACLa_force')
    ACLp = ligament_force(ACLp0, ACLpt, ACLpf, k1ACLp, k2ACLp, 'ACLp', 'ACLp_force')
    PCLa = ligament_force(PCLa0, PCLat, PCLaf, k1PCLa, k2PCLa, 'PCLa', 'PCLa_force')
    PCLp = ligament_force(PCLp0, PCLpt, PCLpf, k1PCLp, k2PCLp, 'PCLp', 'PCLp_force')
    LCL = ligament_force(LCL0, LCLt, LCLf, k1LCL, k2LCL, 'LCL', 'LCL_force')
    MCLa = ligament_force(MCLa0, MCLat, MCLaf, k1MCLa, k2MCLa, 'MCLa', 'MCLa_force')
    MCLo = ligament_force(MCLo0, MCLot, MCLof, k1MCLo, k2MCLo, 'MCLo', 'MCLo_force')
    MCLd = ligament_force(MCLd0, MCLdt, MCLdf, k1MCLd, k2MCLd, 'MCLd', 'MCLd_force')

    N, Nsoa, stress, stress_cartilage, distance_cartilage0, cartilage_femoral, cartilage_tibial, stress_bone, \
        stress_meniscus, distance_meniscus0 = normal_force()
    F_mus, M_mus, soa_mus = external_force(i)

    force, moment, soa = result_of_forces_and_moments(ACLa[0], ACLa[2], ACLp[0], ACLp[2],
                                                      PCLa[0], PCLa[2], PCLp[0], PCLp[2],
                                                      LCL[0], LCL[2],
                                                      MCLa[0], MCLa[2], MCLo[0], MCLo[2], MCLd[0], MCLd[2],
                                                      N, Nsoa, F_mus, soa_mus)

    center = t.contact_point(center)
    moment = moment_of_force(force, soa, center)
    moment[0] = 0
    soa = site_of_action(force, moment)

    return force, moment, soa, center, N, Nsoa, stress, stress_cartilage, ACLa, ACLp, PCLa, PCLp, LCL, MCLa, MCLo, \
        MCLd, F_mus, soa_mus, distance_cartilage0, cartilage_femoral, cartilage_tibial, stress_bone, stress_meniscus, \
        distance_meniscus0
