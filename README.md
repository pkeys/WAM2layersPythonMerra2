# WAM2layersPythonMerra2-beta

RELEASE NOTES
The main developer of the original WAM2layers code is Ruud van der Ent (r.j.vanderent@uu.nl) and the MERRA2, Kenya-specific modifications led by Patrick W. Keys (patrick.keys@colostate.edu). Please let Pat know when you intend to use this Merra2 code for research, I might be able to offer some kind of support if it is in my interest as well, but I cannot guarantee to do any kind of troubleshooting for free as I have to work on my own projects as well.

This code Released under the GNU General Public License. Please cite this Github and 
"Van der Ent, R. J. (2014), A new view on the hydrological cycle over continents, Ph.D. thesis, 96 pp, Delft University of Technology, Delft. 
http://dx.doi.org/10.4233/uuid:0ab824ee-6956-4cc3-b530-3245ab4f32be."

Note:
Please find the original WAM2layers code is found here: https://github.com/ruudvdent/WAM2layersPython
The code here is stable and suitable for use with ERA-Interim data, and includes all necessary download scripts. 
Please see all original Readme and associated documentation on the master branch of this code.

Warning:
This code has been tested for the moisture recycling patterns associated for the country of Kenya ONLY. 
Future work may update the code and test for stability at the Global scale, but this has not been done. 

WAM-2layersMerra2 v1.0-beta | 1-7-2021
- This is a release of the WAM2layersPython code from Ruud van der Ent.
- Major changes have been made from the original code, to the Fluxes_and_States_Masterscript_x, including MERRA2-specific data handling.
- Minor changes were made to the preamble portions of the script in Con_E_Recyc_Masterscript_x, Con_E_Recyc_Output_x, Con_P_Recyc_Masterscript_x, Con_P_Recyc_Output_x, and getconstants. These changes were made to accommodate the different data source (MERRA2), and the focus on Kenya as the region of analysis.


