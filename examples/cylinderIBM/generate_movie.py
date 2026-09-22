#!/usr/bin/env python

import cv2
import glob
import re


# Natural sorting helper function
def natural_sort_key(s):
    return [
        int(text) if text.isdigit() else text.lower()
        for text in re.split(r'([0-9]+)', s)
    ]


# Output video and frame rate
movie_fname = 'movie.mp4'
frame_rate = 50

# Find and naturally sort all PNG files
png_files = sorted(glob.glob('movie/*.png'), key=natural_sort_key)

if not png_files:
    raise FileNotFoundError("No PNG files found in movie/.")

# Read the first image to determine the video size
first_img = cv2.imread(png_files[0])
if first_img is None:
    raise RuntimeError(f"Cannot read image: {png_files[0]}")

height, width = first_img.shape[:2]

# Create the video writer
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
my_video = cv2.VideoWriter(
    movie_fname, fourcc, frame_rate, (width, height)
)

if not my_video.isOpened():
    raise RuntimeError(f"Cannot create video: {movie_fname}")

try:
    for fname in png_files:
        this_img = cv2.imread(fname)
        if this_img is None:
            raise RuntimeError(f"Cannot read image: {fname}")

        # Resize only when necessary
        if this_img.shape[:2] != (height, width):
            this_img = cv2.resize(this_img, (width, height))

        my_video.write(this_img)
finally:
    my_video.release()

print(f"Saved {movie_fname} ({len(png_files)} frames, {frame_rate} FPS)")
