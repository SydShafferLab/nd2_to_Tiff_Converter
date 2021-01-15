% nd2toTiff
% Reads .nd2 files and produces .tif files
% 
% Input
% *|infile|*  - filename or full filepath string of the image stack
% *|outDir|*  - filepath for the Output
% *|nDigits|* - max order of magnitude for naming convention 
%                                ei. 'dapi001'  nDigits = 3
% 
% Output
% *|tiff|* - 3D tiff files named "NAME#" were # coresponds to the Z stack and
%           takes into acount previous files. NAME is taken from the wavelength
%           channel (ie DAPI, CY3,...) and changed so that arjlabimagetools can
%           read it.
%           This is an Example of the naming convention as follows:
%                   Our_names     arj_names
%                   alexa594  ->  alexa
%                   Atto647N  ->  cy
%                   cy3       ->  tmr
%                   700       ->  nir
%          Add to this list in the "channelMap"
% 
% Required:
%            Get bfmatlab   
%                    1)Go to:   https://downloads.openmicroscopy.org/bio-formats/6.0.1/artifacts/
%                    2)Download bfmatlab.zip
%                    3)Unzip and move bfmatlab folder to your MATLAB folder
%                    4)Add bfmatlab path to Matlab     
% Usage
%  >> T1 = nd2toTiff('Fish_Scan1.nd2');               % read image stack in working directory
% 
%  >> T1 = nd2toTiff('/path/to/file/Fish_Scan1.nd2'); % read image stack from specific filepath
% 
%  >> T1 = nd2toTiff('/path/to/file/*.nd2');          % read all nd2 files from specified filepath
% 
% Usage of 'outDir'
% 
%  >> T1 =
%  nd2toTiff('/path/to/file/Fish_Scan1.nd2','outDir','/path/to/outputfile/') 
% 
%  >> T1 = nd2toTiff('/path/to/file/*.nd2','outDir','/path/to/outputfile/')
% 
% Last Update 11/10/2020 Fixed bug in timelapse reader and added some
% print statements
%
% By Raul A. Reyes Hueros (raanreye)



function [] = nd2toTiff(infile,varargin)

% The MAP
channelMap = containers.Map({'Brightfield', 'DAPI','Dapi' 'YFP', 'GFP', 'CY3', 'Cy3', 'cy3', 'A594', 'CY5', 'A647', '700', 'CY7','NIR','OFF','CY3_1'},...
                            {'trans'      , 'dapi', 'dapi','gfp', 'gfp', 'tmr', 'tmr', 'tmr','alexa', 'cy', 'cy'  , 'nir', 'nir','nir','off','gfp'});


% Input check
p = inputParser;
%p.addRequired('Scan_Width', @ischar)
p.addRequired('inDir', @ischar)
p.addParameter('outDir', '', @ischar);
p.addParameter('nDigits', 3, @(x)validateattributes(x,{'numeric'},{'positive','integer'}));

p.parse(infile, varargin{:});


%  ------------------------------------------------------------------------
%                   What Dir to use????  Not the nicest but it works
%  ------------------------------------------------------------------------

% Where to output files
if isempty(p.Results.outDir)       
    outDir = pwd;
else
    outDir = p.Results.outDir;
end

% Where are the files
infile = p.Results.inDir;

% Get Info from infile
[Dir,file_name,c] = fileparts(infile);

if file_name == "*" % for the *.nd2 case
    c = [];
end

if numel(c) == 0                                       % Folder path given
    Files         = ls(Dir);
    file_name_nd2 = string(split(Files));
    TF            = contains(file_name_nd2,'.nd2');
    file_name_nd2 = file_name_nd2(TF);

    file_name = strings(numel(file_name_nd2));
    out_file_name   = strings(numel(file_name_nd2),1);
    out_file_folder = strings(numel(file_name_nd2),1);
    for j = 1:numel(file_name_nd2)
        [~,hold,~]   = fileparts(file_name_nd2(j));
        file_name(j) = hold;
        
        out_file_name(j)   = strcat(outDir,'/',hold,'.nd2');
        out_file_folder(j) = strcat(outDir,'/',hold);
        
        mkdir(out_file_folder(j))
        
        
        %  ----------------------------------------------------------------
        %                  Are there any previous files  
        %  ----------------------------------------------------------------

        % Assign starting number for output files
        if isempty(dir(fullfile(out_file_folder(j), '*.tif')))
             imageCount(j) = 0;
        else
            currentFiles = dir(fullfile(out_file_folder(j), '*.tif'));
            currentFiles = {currentFiles.name};
            currentFileCount = regexp(currentFiles, '\d+', 'match');
            currentFileCount = vertcat(currentFileCount{:});
            currentFileCount = str2num(cell2mat(currentFileCount));
            imageCount(j) = max(currentFileCount); 
        end
        %  ----------------------------------------------------------------
        
    end
    
    outDir = out_file_folder; % update outdir 
    fprintf('\n');
    fprintf('Reading %d Files:  \n',numel(file_name_nd2));
    
else                                                     % File path given
    file_name     = string(file_name);
    file_name_nd2 = string(strcat(file_name,'.nd2'));
    fprintf('\n');
    fprintf('Reading one File:  \n');
    
    out_file_name   = strcat(outDir,'/',file_name,'.nd2');
    out_file_folder = strcat(outDir,'/',file_name);
    outDir          = out_file_folder; % update outdir 
    
    mkdir(out_file_folder)
    
    %  --------------------------------------------------------------------
    %                  Are there any previous files  
    %  --------------------------------------------------------------------

    % Assign starting number for output files
    if isempty(dir(fullfile(outDir, '*.tif')))
         imageCount = 0;
    else
        currentFiles = dir(fullfile(outDir, '*.tif'));
        currentFiles = {currentFiles.name};
        currentFileCount = regexp(currentFiles, '\d+', 'match');
        currentFileCount = vertcat(currentFileCount{:});
        currentFileCount = str2num(cell2mat(currentFileCount));
        imageCount = max(currentFileCount); 
    end
    %  --------------------------------------------------------------------
    
    
end      


%  ------------------------------------------------------------------------
%                          Look Through nd2  
%  ------------------------------------------------------------------------
nDigits = num2str(p.Results.nDigits);
cnt_mult_files = imageCount;


for f = 1:numel(file_name_nd2)
    fprintf('                %s\n',file_name_nd2(f));
    
    cnt_stacks = 0;
    % Read through 
    full_Data_path = fullfile(Dir,file_name_nd2(f));
    reader = bfGetReader(char(full_Data_path));
    omeMeta = reader.getMetadataStore();
    
    Zstacknumb = omeMeta.getImageCount();                    % number of Z stacks
    for i = 1:Zstacknumb % Running through # of z stacks
        
        cnt_stacks = cnt_stacks + 1;
        
        reader.setSeries(i-1)                                % set ith z stack 
        
        if mod(i,25) == 0
            if i == 25
                % Total stack being read
                fprintf('                 Total Stack = %03d\n',Zstacknumb);
            end
            % stack being read
            fprintf('                       Stack = %03d\n',i);
        end
        
        % Getting info from .nd2 file 
        stackSizeC = omeMeta.getPixelsSizeC(i-1).getValue(); % # of wavelength channels
        stackSizeZ = omeMeta.getPixelsSizeZ(i-1).getValue(); % # of Z slices
        
        stackSizeT = omeMeta.getPixelsSizeT(i-1).getValue(); % # of Time slices
        
        % Arrays
        wave_numb = repmat(1:stackSizeC,1,stackSizeZ);     % repeating vector of channel #
        Z_numb    = repmat(1:stackSizeZ,1,stackSizeC);     % repeating vector of Z #

        for ii = 1:stackSizeC  % (# channels in i)
            cnt = 1; %new .tif
            for iii = 1:stackSizeZ  % (# z slices in stack)
                cnt_time = 0; % new time frame
                for iiii = 1:stackSizeT % (# Time frames)
                    
                    if mod(iiii,10) == 0
                        
                        if iiii == 10
                            % Total frames being read
                            fprintf('            Total Frames = %03d\n',stackSizeT);
                        end
                        
                        % Frame being read
                        fprintf('                       Frame = %03d\n',iiii);
                    end
                    
                    % Read plane from series iSeries at Z, C, T coordinates (iZ, iC, iT)
                    iPlane = reader.getIndex(Z_numb(iii) - 1, wave_numb(ii) - 1, iiii - 1) + 1;



    %________________________%%bfGetPlane fails if iPlane is too large________________
    
                    %stack_fig  = bfGetPlane(reader, 1, 1, 1, 1024, 1024);
                    stack_fig  = bfGetPlane(reader, iPlane);
 %____________________________________________________________________________________
 
 
 
 
 
                    %---------------------------------------------------------
                    %                         Save Tiff 
                    %---------------------------------------------------------
                    % Change name 
                    channelName = omeMeta.getChannelName(i-1, ii-1);
                    channelName = channelName.toCharArray';

                    outBaseName = strcat('%s%0', nDigits, 'd.tif'); % set the 00# in file name 
                    outputFileName = fullfile(outDir(f), sprintf(outBaseName,channelMap(channelName),cnt_mult_files(f) + cnt_stacks + cnt_time));
                    
                    cnt_time = cnt_time + 1;

                    if cnt == 1
                        imwrite(stack_fig, outputFileName)
                        cnt = 1 + cnt;
                    else
                        imwrite(stack_fig, outputFileName, 'writemode', 'append')
                        cnt = 1 + cnt;
                    end
                end
            end
        end            
     end
end

fprintf('\n');
fprintf('   %s\n',"Done");
fprintf('\n');
end


