
timespan,'2008-3-26/8',5,/h


thm_load_state, probe = 'd', /get_support_data
thm_load_fgm, lev=1, probe=['d'], /get_support_data,coord='gsm',type='calibrated', suffix='_after_gsm', use_eclipse_corrections=2

;tplot_options, 'title', 'THEMIS FGM Eclipse Corrected Data'
;tplot, ['thd_fgl_after_gsm']

filename = '/home/nithin/Documents/git-repos/energy-height-conversion/Data_Mar_08_Event/thd_fgm_corrected.txt'

;Getting data to easily write into ascii files
get_data, 'thd_fgl_after_gsm', data=dB, dlimits=dl

;Write data to ascii files
header =['        Time              Bx [nT]      By [nT]    Bz [nT]   ']
sdata  ={time:time_string(dB.x), Bx:dB.y(*,0), By:dB.y(*,1), Bz:dB.y(*,2)}
write_ascii_cmdline, sdata, filename, header=header, nrec=nrec

print, 'Filename: ',filename
print, 'The number of data records written from the data structure is: ',nrec

end
