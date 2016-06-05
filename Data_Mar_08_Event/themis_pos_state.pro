
; Routine to calculate MLT, L-shell, R_E, MLAT if THEMIS Spacecraft vs. time 

del_data,'*'

;Set day
timespan, '2008-03-26'

probes=['d','e','c']

for i=0, n_elements(probes)-1 do begin
  probe = probes[i]

  ;Get position in GSM coords
  thm_load_state, probe=probe, datatype='pos' , coord='gsm'

  ;Convert units to RE
  tkm2re, 'th'+probe+'_state_pos'

  ;Get data and put into required form ([time, x, y, z])
  get_data, 'th'+probe+'_state_pos_re', data = d, dlimits=dl
  data = transpose( [[d.x],[d.y]] )

  ;Get l-shell value
  shells = calculate_lshell(data)

  ;Store data
  store_data, 'lshell_value_th'+probe, data = {x:d.x ,y:shells }

endfor

;Calculating the seconds of the year
t0 = cnvtime(2008,03,26,00,00,00)

for i=0, n_elements(probes)-1 do begin
  probe = probes[i]

  ;Transform coordinates to  SM from GSM coords
  cotrans,'th'+probe+'_state_pos','th'+probe+'_state_pos_sm',/gsm2sm
  
  ;Convert to Spherical Coordinates
  xyz_to_polar,'th'+probe+'_state_pos_sm'
  
  
  ;Get data
  ; MLT
  get_data, 'th'+probe+'_state_pos_sm_phi', data = d_phi, dlimits=dl1
  get_data, 'th'+probe+'_state_pos_sm_th', data = d_th, dlimits=dl2
  get_data, 'th'+probe+'_state_pos_sm_mag', data = d_mag, dlimits=dl3
  ;Calculate  MLT
  mlt= ((d_phi.y+180)*(24.0/360))
  mlat= (d_th.y)
  RE= (d_mag.y/6378.1)
 
  ;Store data
  store_data, 'mlt_value_th'+probe, data = {x:d.x ,y:mlt }
  store_data, 'mlat_value_th'+probe, data = {x:d.x ,y:mlat }
  store_data, 'RE_value_th'+probe, data = {x:d.x ,y:RE }

endfor

;Plot
tplot, 'lshell_value_th'+probes
tplot, 'mlt_value_th'+probes
tplot, 'mlat_value_th'+probes
tplot, 'RE_value_th'+probes

; Need to convert to asci files

end
