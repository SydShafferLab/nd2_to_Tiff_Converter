function [] = nd2toTiff(infile,varargin)

% The MAP
channelMap = containers.Map({'Brightfield', 'DAPI', 'YFP', 'GFP', 'CY3', 'Cy3', 'cy3', 'A594', 'CY5', 'A647', '700', 'CY7','NIR'},...
                            {'trans'      , 'dapi', 'gfp', 'gfp', 'tmr', 'tmr', 'tmr','alexa', 'cy', 'cy'  , 'nir', 'nir','nir'});


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
    fprintf('Reading %d Files:  \n',numel(file_name_nd2));
    
else                                                     % File path given
    file_name     = string(file_name);
    file_name_nd2 = string(strcat(file_name,'.nd2'));
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
                
                for iiii = 1:stackSizeT % (# Time frames)
                
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

                    outBaseName = strcat('%s%0', nDigits, 'd.tif');
                    outputFileName = fullfile(outDir(f), sprintf(outBaseName,channelMap(channelName),cnt_mult_files(f) + cnt_stacks));

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
