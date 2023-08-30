import imageio
import numpy as np
from pdf2image import convert_from_path
import os

images = []
path = './gif_workplace/'
# for filename in os.listdir(path):
output = '4D_NLO_proton_profile'
# for i in range(7):
#     filename = 'N_in_r_y_' + str(i) + '.0.pdf'
for y in np.linspace(0., 0.1, 3):

    filename = 'N_y_' + str(round(y, 2)) + '.pdf'
    print(filename)

    images_jpg = convert_from_path(path + filename)
    jpg_filename = filename[:-4] + '.jpg'
    images_jpg[0].save(path + jpg_filename, 'JPEG')

    images.append(imageio.imread(path + jpg_filename))
imageio.mimsave(path + output + '.gif', images, format='GIF', fps=0.6)



