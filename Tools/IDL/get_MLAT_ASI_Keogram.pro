; Extract MLAT corresponding to pixels for the keogram of ASI FYKN and GAJO 

; Loading file containing MLAT, and MLON information for FYKN Camera
filename1 = '/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/Raw/Keograms/thg_asf_fykn_aacgm.tplot'

; Loading file containing MLAT, and MLON information for GAKO Camera
filename2 = '/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/Raw/Keograms/thg_asf_gako_aacgm.tplot'

; Loading the tplot files
tplot_restore, filename=filename1
tplot_restore, filename=filename2

;Getting data to easily write into ascii files
get_data, 'thg_asf_fykn_aacgm', data=dFYKN, dlimits=dl
get_data, 'thg_asf_gako_aacgm', data=dGAKO, dlimits=dl
pixels=[1:256]
;Write data to ascii files
header =['Pixels MLAT']
fykndata  ={pixel:pixels, MLAT:dFYKN.MLAT(128,*)}
gakodata  ={pixel:pixels, MLAT:dGAKO.MLAT(128,*)}

;filename to write the ascii file
filename = '/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/Raw/Keograms/FYKN_MLAT.txt'
write_ascii_cmdline, fykndata, filename, header=header, nrec=nrec
print, 'Filename: ',filename
print, 'The number of data records written from the data structure is: ',nrec

filename = '/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/Raw/Keograms/GAKO_MLAT.txt'
write_ascii_cmdline, gakodata, filename, header=header, nrec=nrec
print, 'Filename: ',filename
print, 'The number of data records written from the data structure is: ',nrec


end