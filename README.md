# cflp-ice
A suite of programs for managing experimental information, image data, internal wave detection and ice edge-detection, via GUIs.
---
# Requirements
 1. Matlab 2022a or later.
 2. The following packages/ add-ons:
	- Curve Fitting Toolbox
	- Image Acquisition Toolbox
	- Image Processing Toolbox
	- Parallel Computing Toolbox
	- Signal Processing Toolbox
	
---
# Primary contents:
There are three general programs:
	1. _**app_experiment**_:
		Responsible for managing, designing, navigating experiment directories.
	2. _**app_meltingDetection**_:
		Responsible for viewing and rendering image-data where wavefields are excited through synthetic-schlieren imagery.
	3. _**app_waveDetection**_:
		Responsible for viewing and processing image-data where the subject is a slab of ice and edge-detection is required.

# Secondary contents:
A suite of mostly self-contained classes are used to operate the general programs. The primary classes, respective to the order of the general programs above are:
	1. _**labExperiment**_
		Skeleton class containing basic laboratory information and records (e.g., densimeter and conductivity-temperature probe outputs).
	2. _**edIce**_
		Superclass for ice edge-detection operations and melting rate calculations.
	3. _**iwdata**_
		Superclass for internal wave detection operations and calculations.
		