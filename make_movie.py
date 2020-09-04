import os
import cv2
import glob

# To convert avi to mp4 use ffmpeg:
# ffmpeg -i DAH_Complete_Sorting_1.avi DAH_Complete_Sorting_1.mp4

# params
fps = 2
foldername = "DAH_Lipid_Bilayer_1_persimg"
output_file_name = "DAH_Lipid_Bilayer_1.avi"

fnames = []
for filepath in sorted(glob.glob(path)):
    fname = filepath.split(os.sep)[1]
    basename = fname.split(".")[0]
    fnames.append(int(basename))
fnames.sort()

video = None
vInit = False

for fn in fnames:

    filepath = foldername + os.sep + repr(fn) + ".png"

    if vInit == False:

        img = cv2.imread(filepath)
        height, width, layers = img.shape

        fourcc = cv2.VideoWriter_fourcc('m', 'p', '4', 'v')
        video = cv2.VideoWriter(output_file_name, fourcc, fps, (width, height), True)

        video.write(img)
        vInit = True

    else:

        video.write(cv2.imread(filepath))

if video is not None:

    cv2.destroyAllWindows()
    video.release()
