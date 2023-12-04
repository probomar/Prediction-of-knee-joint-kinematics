import cv2
# from input import *
import os

def generate_video(path, video_name, img_folder):
    os.chdir(path)

    images = [img for img in os.listdir(img_folder)
              if img.endswith('.jpg') or img.endswith('.jpeg') or img.endswith('.png')]
    frame = cv2.imread(os.path.join(img_folder, images[0]))
    height, width, layers = frame.shape
    video = cv2.VideoWriter(video_name, 0, 3, (width, height))

    for image in images:
        video.write(cv2.imread(os.path.join(img_folder, image)))

    cv2.destroyAllWindows()
    video.release()


# path = folder_name + '/flex'
path = 'C:/Users/pc/PycharmProjects/Prediction-of-knee-joint-kinematics/video'
video_name = '_flex.avi'
# video_name = '_stress.avi'
image_folder = '.'

os.chdir(path)

mean_height = 0
mean_width = 0

num_of_img = len(os.listdir('.'))
print('num_of_img', num_of_img)

generate_video(path, video_name, image_folder)
