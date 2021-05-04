
![](https://raw.github.com/smutao/PyVibMS/master/for-readme/logo.png)

PyVibMS
=========================
A PyMOL plugin for visualizing vibrations in molecules and solids

Updates
-----
May 3, 2021 - Add support for ORCA 4 and Q-Chem 4/5

May 2, 2021 - Add support for Prof. Grimme's xtb program  



User Guide
------

介绍PyVibMS使用的中文文章 请前往 http://bbs.keinsci.com/thread-22835-1-1.html

Our paper on PyVibMS is now published. 🎉 

PyVibMS: a PyMOL plugin for visualizing vibrations in molecules and solids,
Journal of Molecular Modeling, 2020, 26, 290. https://link.springer.com/article/10.1007/s00894-020-04508-z

This paper serves as a detailed manual of PyVibMS and a few typographical errors need to be corrected:
1) Page 2

<img src="https://raw.github.com/smutao/PyVibMS/master/for-readme/misprint/1.png" width="450">

2) Page 3

<img src="https://raw.github.com/smutao/PyVibMS/master/for-readme/misprint/2.png" width="450">

3) Page 5

<img src="https://raw.github.com/smutao/PyVibMS/master/for-readme/misprint/3.png" width="450">


Installation
=========================
1) Please install PyMOL 2.x and download the zip file of this repository

2) Open PyMOL 2.x 

3) Click "Plugin" -> "plugin manager" -> "Install New Plugin" -> "Choose", then choose the "src/\_\_init\_\_.py" file

4) The "PyVibMS" will be installed to PyMOL and show up in the "Plugin" drop-in menu

PyVibMS runs on various operating systems.

Windows

<img src="https://raw.github.com/smutao/PyVibMS/master/gallery/pyvibms_on_win10.png" width="750">

Linux

<img src="https://raw.github.com/smutao/PyVibMS/master/gallery/pyvibms_on_centos.png" width="750">

Mac OS

<img src="https://raw.github.com/smutao/PyVibMS/master/gallery/pyvibms_on_macos.png" width="750">



Supported quantum chemistry programs 
=========

Natively supported
----

* Gaussian 09/16
* ORCA 4
* xtb 
* Q-Chem
* CRYSTAL17
* VASP 

Generically supported
----

* [Adf](http://www.scm.com/)
* [BDF](http://182.92.69.169:7226)
* [CFour](http://www.cfour.de/)
* [Columbus](http://www.univie.ac.at/columbus/)
* [CP2k](http://www.cp2k.org/)
* [Dalton](http://daltonprogram.org/)
* [deMon2k](http://www.demon-software.com/public_html/)
* [Dmol3](http://accelrys.com/)
* [Firefly](http://classic.chem.msu.su/gran/gamess/)
* [Gamess](http://www.msg.chem.iastate.edu/gamess/)
* [Gamess-UK](http://www.cfs.dl.ac.uk/)
* [Molcas](http://www.molcas.org/) and [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas)
* [Molpro](http://www.molpro.net/)
* [Mopac](http://openmopac.net/)
* [NWChem](http://www.nwchem-sw.org/index.php/Main_Page)
* [Pqs](http://www.pqs-chem.com/)
* [Psi4](http://www.psicode.org/)
* [Turbomole](http://www.cosmologic.de/)
* [Quantum ESPRESSO](https://www.quantum-espresso.org/)
* [DFTB+](https://dftbplus.org/)
* and more


Other analysis tools
------
* [LModeA](https://onlinelibrary.wiley.com/doi/full/10.1002/wcms.1480)
* [Phonopy](https://phonopy.github.io/phonopy/)


Tips
========
* When using the open-source version of PyMOL, one may find it difficult to export the vibration animation directly as a GIF or QuickTime movie due to the missing encoder, a work-around is to export the animation as "PNG images" (by clicking "File"->"Export Movie As"->"PNG Images..."). Then one can use a third-party tool like [ezgif](https://ezgif.com/apng-to-gif) to combine these PNG images as animated GIF image.  
Here is an example  
<img src="https://raw.github.com/smutao/PyVibMS/master/gallery/animations/gif-format.gif" width="350">



