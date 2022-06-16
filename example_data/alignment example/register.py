# imageJ assisted manual registering

from __future__ import print_function
import cv2
import numpy as np
import pandas as pd
def alignImages(im_fit, im_ref, p_fit, p_ref):
    if ( len(p_ref) != len(p_fit) ): print("error! The amount of landmark coordinates differ!", len(p_ref), len(p_fit))
    # Find homography
    h, mask = cv2.findHomography(p_fit, p_ref, cv2.RANSAC)
    np.savetxt('h.csv', h)
    # Use homography
    height, width, channels = im_ref.shape
    im1Reg = cv2.warpPerspective(im_fit, h, (width, height))

    return im1Reg, h

if __name__ == '__main__':
    refFilename = "imageV51"
    imFilename = "Klinker LA-mapAB2_005"
    outFilename = "aligned"

    print("Reading reference image : ", refFilename)
    im_ref = cv2.imread(refFilename + ".tif", cv2.IMREAD_COLOR)
    p_ref = pd.read_csv(refFilename + '.csv', delimiter=',', usecols=['X', 'Y']).to_numpy()

    print("Reading image to align : ", imFilename);
    im_fit = cv2.imread(imFilename + ".tif", cv2.IMREAD_COLOR)
    p_fit = pd.read_csv(imFilename + '.csv', delimiter=',', usecols=['X', 'Y']).to_numpy()

    #p_fit = np.genfromtxt('123.csv', delimiter=',')
    #print(p_fit)
    print("Aligning images ...")
    # Registered image will be resotred in imReg.
    # The estimated homography will be stored in h.
    imReg, h = alignImages(im_fit, im_ref, p_fit, p_ref)

    # Write aligned image to disk.
    print("Saving aligned image : ", outFilename);
    cv2.imwrite(outFilename + ".tif", imReg)

    # Print estimated homography
    print("Estimated homography : \n",  h)
