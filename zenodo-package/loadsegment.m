function [segment] = loadsegment(segnum)
%LOADSEGMENT Summary of this function goes here
%   Detailed explanation goes here

%Constants
sampFreq = 4096;

%Load HDF5 files and segment json files
segmentjson = loadjson('segments.json');
hdf5filejson = loadjson('hdf5files.json');

segparams = segmentjson{segnum};
segstart = segparams.GPS_start;
segend = segparams.GPS_end;
files = segparams.files;

%Get File data and load segment
if length(files) == 1
    
    hdf5file = hdf5filejson{files};
    filename = hdf5file.filename;
    GPSstart = hdf5file.GPS_start;
    GPSend = hdf5file.GPS_end;
%     url = hdf5file.url;
    
    %Get segment indices
    startidx = sampFreq*(segstart- GPSstart);
    endidx = sampFreq*(segend- GPSstart);

    %Load data from file and create segment data vector

    filedata = h5read(filename, '/strain/Strain')';
    segment = filedata(startidx+1:endidx);
%     psdnum = filenums;
else
    segment = [];
%     psdnum = files(1);
    for i = 1:length(files)
        hdf5file = hdf5filejson{files(i)};
        filename = hdf5file.filename;
        GPSstart = hdf5file.GPS_start;
        GPSend = hdf5file.GPS_end;

        if GPSend > segstart && GPSstart < segstart && segend > GPSend
            %Get segment indices
            startidx = sampFreq*(segstart- GPSstart);
            endidx = sampFreq*(GPSend -GPSstart);
            %Load data from file and create segment data vector
            filedata = h5read(filename, '/strain/Strain')';
            seg = filedata(startidx:endidx);
            segment = [segment, seg];
        end
        if GPSstart < segend && GPSend > segend && segstart < GPSstart
            %Get segment indices
            startidx = 1;
            endidx = sampFreq*(segend - GPSstart);
            %Load data from file and create segment data vector
            filedata = h5read(filename, '/strain/Strain')';
            seg = filedata(startidx:endidx-1);
            segment = [segment, seg];
        end
    end
end



end

