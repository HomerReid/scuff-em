real_MFull=hdf5read('Full.hdf5','real_MFull');
imag_MFull=hdf5read('Full.hdf5','imag_MFull');
MFull=real_MFull + sqrt(-1)*imag_MFull;
clear real_MFull; clear imag_MFull;
RHSFull=hdf5read('Full.hdf5','RHSFull');
KNFull=hdf5read('Full.hdf5','KNFull');
