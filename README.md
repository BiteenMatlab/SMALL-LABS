# SMALL-LABS

The SMALL-LABS (Single-Molecule Accurate Localization by Local Background Subtraction) algorithm, accurately locates and measures the intensity of single molecules, regardless of the shape or brightness of the background.

The program can also fit single molecules without doing background subtraction. See the User Guide for more details.

Written by Benjamin P Isaacoff at the University of Michigan.

This algorithm is described in B.P. Isaacoff, Y. Li, S.A. Lee and J.S. Biteen, “SMALL-LABS: An algorithm for measuring single-molecule intensity and position in the presence of obscuring backgrounds,” _Biophysical Journal_, in press (**2019**). https://doi.org/10.1016/j.bpj.2019.02.006

## Installation

Download the entire folder and unzip if you downloaded the .zip folder. Change the working directory in Matlab to this folder and call the functions in the Matlab command window as described in the User Guide.

## Usage

See the Quick Start Guide for a quick introduction to using *SMALLLABS_main*. 

See the User Guide for the details. Briefly, the function *SMALLLABS_main* is a wrapper for the other code to perform all of the steps in the correct order. Simply run *SMALLLABS_main* by specifying the directory containing your movies, specify the three required parameters, and any optional parameters, then run it and click to choose the movies you want to fit. Or run the various programs independently.

## Contributing

Please inform us (isaacoff@umich.edu) before making any changes, then follow the directions below: 
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## Credits

All individual programs should have all their individual attributions in place including authors. 

This code was developed with support from the National Science Foundation (NSF grant CHE-1252322).

The development of this code is greatly indebted to the work of David J Rowland (often referred to as DJR in the code). In addition to containing some functions written by him, I’ve borrowed a lot of code snippets from his programs.

SMALL-LABS uses a number of open-source codes and algorithms:

**_TiffStack_** by DR Muir and BM Kampa 

DR Muir and BM Kampa. 2015. FocusStack and StimServer: A new open source MATLAB toolchain for visual stimulation and analysis of two-photon calcium neuronal imaging data, **Frontiers in Neuroinformatics** 8 *85*. DOI:10.3389/fninf.2014.00085

TIFFStack by Dylan Muir is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License.
Based on a work at http://github.com/DylanMuir/TIFFStack

**_saveastiff_** by YoonOh Tak

Copyright (c) 2012, YoonOh Tak
All rights reserved.
Available at https://www.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack

**_bpass_** by John C. Crocker and David G. Grier  

Copyright (c) 1997, John C. Crocker and David G. Grier
Available at http://www.physics.emory.edu/faculty/weeks//idl/

**_MLEwG_** by KI Mortensen, LS Churchman, JA Spudich, H Flyvbjerg

K. I. Mortensen, L. S. Churchman, J. A. Spudich, and H. Flyvbjerg, Nat. Methods **7**, 377 (2010) doi:10.1038/nmeth.1447 

**_gaussFit_** by David J Rowland

David J.Rowland, Julie S.Biteen. Measuring molecular motions inside single cells with improved analysis of single-particle trajectories. Chemical Physics Letters, **674**, 173-178, 2017. DOI:10.1016/j.cplett.2017.02.052

The code 'gaussFit.m' should be considered 'freeware'- and may be distributed freely in its original form when properly attributed.

**_Track_3D2_** by David J Rowland 

David J.Rowland, Julie S.Biteen. Measuring molecular motions inside single cells with improved analysis of single-particle trajectories. Chemical Physics Letters, **674**, 173-178, 2017. DOI:10.1016/j.cplett.2017.02.052

The code 'Track_3D2.m' should be considered 'freeware'- and may be distributed freely in its original form when properly attributed

**_hungarian_** by Yi Cao 

Available at https://www.mathworks.com/matlabcentral/fileexchange/20652-hungarian-algorithm-for-linear-assignment-problems-v2-3

**_gpufit_** by Adrian Przybylski, Björn Thiel, Jan Keller-Findeisen, Bernd Stock, and Mark Bates

Gpufit: An open-source toolkit for GPU-accelerated curve fitting
Adrian Przybylski, Björn Thiel, Jan Keller-Findeisen, Bernd Stock, and Mark Bates
Scientific Reports, vol. 7, 15722 (2017); doi: https://doi.org/10.1038/s41598-017-15313-9

MIT License.
Copyright (c) 2017 Mark Bates, Adrian Przybylski, Björn Thiel, and Jan Keller-Findeisen


## License

                      GNU GENERAL PUBLIC LICENSE
                       Version 3, 29 June 2007

  See LICENSE.txt
