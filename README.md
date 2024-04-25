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
	- Skeleton class containing basic laboratory information and records (e.g., densimeter and conductivity-temperature probe outputs).
 2. _**edIce**_
	- Superclass for ice edge-detection operations and melting rate calculations.
 3. _**iwdata**_
	- Superclass for internal wave detection operations and calculations.

# Screenshots
![im2](https://github.com/sweetjim/cflp-ice/assets/60765374/695e8baf-ec8c-4842-a883-36d56120e04c)
![im4a](https://github.com/sweetjim/cflp-ice/assets/60765374/da7d7bdb-64de-4e41-8643-acd85eba88ef)
![im4b](https://github.com/sweetjim/cflp-ice/assets/60765374/66271d55-0fa2-4613-bb07-f8bc73fbca62)
![im5](https://github.com/sweetjim/cflp-ice/assets/60765374/ff33eff3-eb8f-48e8-933c-1c4d263689ed)
![im6](https://github.com/sweetjim/cflp-ice/assets/60765374/405a1f65-80ac-4bee-a13a-695622db6b94)
![im7](https://github.com/sweetjim/cflp-ice/assets/60765374/87c76d78-e3d8-4bf0-adec-6da7e0fc2a25)
![im3](https://github.com/sweetjim/cflp-ice/assets/60765374/3c86cebd-a215-4b91-9dc1-ae2f9081754e)	