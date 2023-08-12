import imageio
from pdf2image import convert_from_path
import os

images = []
path = '/gif_workplace/'
for filename in os.listdir(path):
    print(filename)

    images = convert_from_path(path + filename)
    jpg_filename = filename[:-4] + '.jpg'
    images[0].save(path + jpg_filename, 'JPEG')

    images.append(imageio.imread(path + jpg_filename))
imageio.mimsave(path + filename[:-4] + '.gif', images, fps=3, )



