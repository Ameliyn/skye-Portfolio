This folder contains all of my advanced graphics projects. They have semi-sensible names. Each folder should have a README detailing exactly what was accomplished in each project.
Some folders will have subfolders with XWD files. To display those files, you must be on a linux environment and then run ./video. Not all scripts in this folder are well
commented. Some of them were used as proof of concept and then adapted into larger scripts and forgotten about.

This folder also contains a number of files required for compilation of the rest of the projects. To compile a project normally you run the C compiler with the tags -lm -lX11.

THIS FOLDER IS MEANT TO BE VIEWED ON A LINUX MACHINE.
If you do not have a linux machine to view this folder, I advise downloading the files into a free Repl.It repository and then recompiling the scripts as necessary.
This will lose some quality, but the general gist of each project can be attained there. There are aspects that can be viewed on any machine: gifs, jpgs



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