## nd2toTiff
Converts .nd2 to .tiff

# General
The .nd2 acquired from Nikon Elements can have the following properties: single image, scan (x,y), stack (z),  timelapse (T),  multiple channels (c) or a combination scan + stack + channel + Time (x,y,z,c,T). The nd2toTiff can convert all of these, as long as they are less than ~2GB. This code requires an open source package bfmatlab to run (currently using bfmatlab 6.5.1). 

#Required before running:
Get bfmatlab   
1) Go to:   https://www.openmicroscopy.org/bio-formats/downloads/

2) Download bfmatlab.zip
3) Unzip and move bfmatlab folder to your MATLAB folder (make your life easier and move nd2toTiff in the bfmatlab folder).

4) Add bfmatlab path to Matlab as follows:

a) Go to "File->Set Path" from within MATLAB or type 
                            "pathtool" at the MATLAB prompt.
                            
b) Use the "Add" button to add your desired folder(s) to 
                            the MATLAB path.
                            
c) Click "Save" so that this path is used in future 
                            MATLAB sessions.


# Input/Output

#Input
Required
infile     -> filename or full filepath string of the image stack

Optional
outDir   -> filepath for the Output
nDigits  -> max order of magnitude for naming convention
                               ei. 'dapi001'  nDigits = 3

Output
tif          ->    3D tiff files named "NAME#" were # coresponds to the Z stack and
                  takes into account previous existing files. NAME is taken from the wavelength
                  channel (ie DAPI, CY3,...) and changed so that arjlabimagetools can
                  read it. A folder is made for each file read. 

This is an Example of the naming convention as follows:
                  Our_names     arj_names
                  alexa594  ->  alexa
                  Atto647N  ->  cy
                  cy3            ->  tmr
                  700            ->  nir
Add to this list in the "channelMap" if you have a strange naming convention 

 
#Usage

nd2toTiff('Fish_Scan1.nd2');   % read image stack in working directory and output in working directory

nd2toTiff('/path/to/file/Fish_Scan1.nd2'); % read image stack from specific filepath and output in working directory

nd2toTiff('/path/to/file/*.nd2');  % read all nd2 files from specified filepath and output in working directory

#Usage of 'outDir'

nd2toTiff('/path/to/file/Fish_Scan1.nd2','outDir','/path/to/outputfile/') % read image stack from specific filepath and output file to outDir

nd2toTiff('/path/to/file/*.nd2','outDir','/path/to/outputfile/') % read all nd2 files from specified filepath and output file to outDir

#Usage of ‘nDigit'

nd2toTiff('/path/to/file/Fish_Scan1.nd2','outDir','/path/to/outputfile/',’nDigits’,8) % file NAME# will have 8 digits, ei. 'dapi00000001'
