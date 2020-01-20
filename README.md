# SMF Package
The SMF package is a package includes all code in Matlab for the use and application of the passive stochastic matched filter to the detection of Antarctic blue whale (ABW) calls. 

Author: Léa Bouffaut, Ph.D. [Personal website](https://leabouffaut.home.blog) | [Researchgate](https://www.researchgate.net/profile/Lea_Bouffaut)

This work was conducted during my Ph.D. financed by the french Naval Academy (Institut de Recherche de l'Ecole Navale - Brest, France).

The method is fully described in:

L. Bouffaut, R. Dreo, V. Labat, A. Boudraa and G. Barruol 'Passive stochastic matched filter for antarctic blue whale call detection,' in J. Acoust. Soc. Am, 144(2) (2018).  https://doi.org/10.1121/1.5050520

and any use of this material should refer accordingly. Part of the SMF code was initially written by Gregory Julien, PhD. Julien, G. (2012). 'Filtrage Stochastique et amélioration des performances des systèmes de positionnement d'engins sous-marins en milieu bruyant' (Doctoral dissertation, Toulon).

# This package contains:
1. Main Matlab programs (<i>Offline_savefilterbank</i>, <i>Offline_save_Zcall</i>, <i>Online_application</i>)
1. The same programs but as Matlab notebooks (Matlab Live editor) and their pdf in the <i>Matlab_live</i> folder
1. A function folder (<i>Functions</i>)
1. A folder of saved matrices (<i>Offline_saved</i>)
1. A 24h recording with ABW calls ((<i>RR44_2013_D151.wav</i>).

# Data
A small toy dataset with ABW calls at various SNR is provided (<i>RR44_2013_D151.wav</i>). It consists of a 24h record from the Ocean Bottom Seismomter RR44 deployed during the [RHUM-RUM](http://www.rhum-rum.net/en/) experiment, in the western Indian Ocean. All recordings are hosted and can be downloaded from [RESIF](https://www.resif.fr/) web services.

RHUM-RUM on Git: 
- Global project: https://github.com/rhum-rum
- OBS orientation by John-Robert Scholz: https://gitlab.com/johnrobertscholz/ppol
- OBS for ship noise by Alister Trabattoni: https://github.com/atrabattoni/obsea

RHUM-RUM on researchgate: https://www.researchgate.net/project/RHUM-RUM

# How to run the code

