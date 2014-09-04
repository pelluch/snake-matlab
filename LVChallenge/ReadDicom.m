% Example to read dicom

filename = 'C:/DocsMaracuya/Datos/LVChallenge/ShortAxis/SC-N-9/IM-1031-0001.dcm' ;
[D, info] = ReadData3D(filename);
imlook3d(D);