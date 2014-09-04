#########################################################################
#		MEDICAL IMAGE PROCESSING TOOLKIT 1.0			#
#########################################################################

Written by:

	Alberto Gomez - alberto.gomez@kcl.ac.uk
	Division of Imaging Sciences and Biomedical Engineering
	King's College London, UK

OpenSource code under the BSD-2 license 2010-2013

----------------------------- DISCLAIMER --------------------------------
This software is distributed with no warranties.
This software was designed for research purposes and not for clinical use.
Please, feel free to distribute, modify and improve this software but please 
keep this disclaimer notice and respect authorship.
--------------------------------------------------------------------------

INTRODUCTION

This is a package of matlab classes and functions intended to facilitate research with medical images and triangulated meshes, as a very rudimentary version of image and mesh classes from ITK for matlab.

SUMMARY

The package is organised as follows:

/class_image/
	contains base classes for storing nD scalar and vector images. It also includes a convenience class for echocardiographic images.

/class_mesh/
	contains base classes for 3D triangulated meshes.
	
/IO/
	contains a number of functions for image and mesh input-output. Some of the code has been borrowed from the  ReadData3D_version1 package by  Dirk-Jan Kroon, which can be found in the matlab exchange.

/sources/
	Contains a number of convenient functions to create common geometrical images or meshes: circles, spheres, prisms, arcs, etc.

/visualization/
	Contains functions to visualize meshes and other geometrical figures: circles, planes, points, etc. There is no viewer for the images; there are plenty in the matlab exchange and images can be written to popular formats and viewed with external software.

/transformations/
	Convenience functions to perform common processing operations: resample, transform, and reslice.

/test/ 
	Contains scripts showing the use of the package. 

DOCUMENTATION

A proper documentation is currently unavailable. However, some files in the *test* folder describe most of the functionalities and typical uses.


