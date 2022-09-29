This project was used to learn how to create a video player that reads in a folder of XWD files and displays them in a video format.

1. get_image_from_file.c
    This script was an exploration about how to read XWD images from files.

2. growing_sun_images.c
    This script creates a circle that grows and prints out different versions into XWD files.

3. video_player.c
    This is the original script that creates the video executable files found throughout the repository.


Notes on the ./video executable:

THIS ONLY WORKS ON A LINUX ENVIRONMENT

To use the video script call in the console "./video 800 800 prefix_name starting_frame ending_frame 0 40000"

These parameters are
  ./video        --> The script
  800 800        --> The screen size
  prefix_name    --> The folder and XWD file prefix (ex. "SpaceStation/vids/spaceod"). This prefix should be the filepath and generic file name (not the numbers or extension)
  starting_frame --> The file number to start with (frequently 0)
  ending_frame   --> The file number to end on (frequently 119), in general whichever file is at the end of the desired folder is a good choice
  0              --> whether or not to go back and forth, or just loop (most videos loop seamlessly, but you can change this to 1 to go forwards through the video and then rewind back through).
  40000          --> the time (in ms) between each frame (50000 will seem slower, 30000 will seem faster).

Note: you can call this executable from anywhere as long as your prefix_name is an accurate file path from your current location.