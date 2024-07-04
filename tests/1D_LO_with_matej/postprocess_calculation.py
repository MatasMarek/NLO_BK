import numpy as np
import os, shutil

def delete_folder(folder):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

def save_grid_csv(name, data, path_to_output):
    with open(path_to_output + "/" + name + '.csv', 'w') as file:
        for line in data:
            file.write(str(line))
            file.write('\n')
    return True

def postprocess_run():
    run_name = 'NLO_1D_data_run'

    path_to_output = 'for_matej/' + run_name
    if not os.path.isdir(path_to_output):
        os.mkdir(path_to_output)

    path_to_calculation = '/Users/mmatasadmin/Library/Mobile Documents/com~apple~CloudDocs/code/NLO_BK/output/'
    calculation = np.load(path_to_calculation + run_name + '/calculation.npy', allow_pickle=True).item()

    delete_folder(path_to_output)
    shutil.copyfile(path_to_calculation + run_name + '/input.py', path_to_output + '/input.py')
    shutil.copyfile(path_to_calculation + run_name + '/const.py', path_to_output + '/const.py')

    if not os.path.isdir(path_to_output + '/grids'):
        os.mkdir(path_to_output + '/grids')
    if not os.path.isdir(path_to_output + '/N'):
        os.mkdir(path_to_output + '/N')

    # round the grid values to 2 decimal points
    calculation['grid']['grid_in_Y'] = np.round(calculation['grid']['grid_in_Y'], 2)


    save_grid_csv('y', calculation['grid']['grid_in_Y'], path_to_output + '/grids')
    save_grid_csv('r', calculation['grid']['grid_in_r'], path_to_output + '/grids')
    if calculation['dimensionality_of_N'] >= 2:
        save_grid_csv('b', calculation['grid']['grid_in_b'], path_to_output + '/grids')
    if calculation['dimensionality_of_N'] >= 3:
        save_grid_csv('t', calculation['grid']['grid_in_theta'], path_to_output + '/grids')
    if calculation['dimensionality_of_N'] == 4:
        save_grid_csv('p', calculation['grid']['grid_in_phi'], path_to_output + '/grids')




    for y_ind in range(len(calculation['grid']['grid_in_Y'])):
        if y_ind > calculation['y_ind']:
            break
        y = round(calculation['grid']['grid_in_Y'][y_ind], 2)
        # reformat y such that it is always two digits before and after the decimal point
        y = str(y).split('.')
        y[0] = y[0].zfill(2)
        y[1] = y[1].ljust(2, '0')
        y = '.'.join(y)

        # replace decimal point with letter 'p' in y
        y = str(y).replace('.', 'p')

        print(y)

        with open(path_to_output + "/N/" + y + '.csv', 'w') as file:
            for r_ind in range(len(calculation['grid']['grid_in_r'])):
                if calculation['dimensionality_of_N'] == 1:
                    file.write(str(calculation['N'][y_ind][r_ind]))
                    file.write('\n')
                elif calculation['dimensionality_of_N'] >= 2:
                    for b_ind in range(len(calculation['grid']['grid_in_b'])):
                        if calculation['dimensionality_of_N'] == 2:
                            file.write(str(calculation['N'][y_ind][r_ind][b_ind]))
                            file.write('\n')
                        elif calculation['dimensionality_of_N'] >= 3:
                            for t_ind in range(len(calculation['grid']['grid_in_theta'])):
                                if calculation['dimensionality_of_N'] == 3:
                                    file.write(str(calculation['N'][y_ind][r_ind][b_ind][t_ind]))
                                    file.write('\n')
                                elif calculation['dimensionality_of_N'] == 4:
                                    for p_ind in range(len(calculation['grid']['grid_in_phi'])):
                                        file.write(str(calculation['N'][y_ind][r_ind][b_ind][t_ind][p_ind]))
                                        file.write('\n')


postprocess_run()



