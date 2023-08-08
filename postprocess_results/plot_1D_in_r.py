import matplotlib.pyplot as plt
import os


def plot_1d_in_r(N_in_r, grid_in_r, y, calculation):
    if not os.path.isdir('../output/' + calculation['run_name'] + '/plots/1D_in_r'):
        os.mkdir('../output/' + calculation['run_name'] + '/plots/1D_in_r')
    plt.plot(grid_in_r, N_in_r, label='Y = ' + str(y))
    plt.legend(loc='lower right', shadow=False, frameon=False)
    # set the x-axis to log scale
    title = str(calculation['dimensionality_of_N']) + 'D BK equation at ' + str(calculation['order_of_BK']) + ' order'
    plt.title(title)
    plt.xscale('log')
    plt.ylabel('N(r, Y)')
    plt.xlabel('r [GeV$^{-1}$]')
    plt.savefig('../output/' + calculation['run_name'] + '/plots/1D_in_r/N_in_r_y_' + str(y) + '.pdf', format='pdf')
    plt.show()
    plt.close()
    return