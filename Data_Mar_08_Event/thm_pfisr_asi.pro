;  Create ASI plots
;  By Nithin  Sivadas, 17th July 2016
;  If disp is not set, then it will display as a plot on the screen
;  if it is set to any value, then it will print as postscript
;  Example command: thm_pfisr_asi,show_time='2008-03-26/12:30:00',disp=1,filename='aurora_12_30_UT.ps'
pro thm_pfisr_asi, show_time=show_time, disp=disp, filename=filename


if keyword_set(disp) then begin 
set_plot, 'ps'
device, /portrait, xsize=10, ysize=10, xoffset=0, yoffset=0, /inches, filename=filename
device, /color, bit=8
device, font_size=13
device, font_index=18
device, /Times
endif


 timespan,'2008-3-26/8',5,/h
 
 ;show_time='2008-03-26/11:36:00'
 
 ; Downloading ASCI high resoution data
 thm_asi_create_mosaic,show_time,show=['KIAN','MCGR','FYKN','GAKO','INUV','WHIT'],central_lon=-150,central_lat=65,scale=2e7,projection='AzimuthalEquidistant', disp=disp
 
 ; Downloading THEMIS coordinates
 thm_load_state, probe =['d','e'], coord='gsm', suffix='_gsm'
 
 ; Calculating the foot prints
 ttrace2iono,'thd_state_pos_gsm', newname='thd_ifoot_geo', external_model='t89', par=2.0D,/km,in_coord='gsm',out_coord='geo'
 ttrace2iono,'thc_state_pos_gsm', newname='thc_ifoot_geo', external_model='t89', par=2.0D,/km,in_coord='gsm',out_coord='geo'
 ;ttrace2iono,'the_state_pos_gsm', newname='the_ifoot_geo', external_model='t89', par=2.0D,/km,in_coord='gsm',out_coord='geo'
 
 ; Converting it to lat and long coordinates
 get_data, 'thd_ifoot_geo', data=d
 lon =!radeg*atan(d.y[*,1],d.y[*,0])
 lat =!radeg*atan(d.y[*,2],sqrt(d.y[*,0]^2+d.y[*,1]^2))
 plots,lon,lat, color=FSC_Color('red'), linestyle=2, thick=4
 
 get_data, 'the_ifoot_geo', data=d1
 lon1 =!radeg*atan(d1.y[*,1],d1.y[*,0])
 lat1 =!radeg*atan(d1.y[*,2],sqrt(d1.y[*,0]^2+d1.y[*,1]^2))
 plots,lon1,lat1, color=FSC_Color('black'), linestyle=3, thick=4
 
 ;get_data, 'thc_ifoot_geo', data=d2
 ;lon2 =!radeg*atan(d2.y[*,1],d2.y[*,0])
 ;lat2 =!radeg*atan(d2.y[*,2],sqrt(d2.y[*,0]^2+d2.y[*,1]^2))
 ;plots,lon2,lat2, color=FSC_Color('red')
 
 ; Identifying THEMIS-D location at the current time
 min_diff=min(abs(d.x-time_double(show_time)),index)
 
 ; THEMIS-E
 min_diff_1=min(abs(d1.x-time_double(show_time)),index1)

 ; THEMIS-C
 ;min_diff_2=min(abs(d2.x-time_double(show_time)),index2)
 
 ;xyouts,lon[index]+0.5,lat[index]+0.1,'THEMIS-D',/data,charsize=2, color=FSC_Color('red'), charthick=4
 xyouts,0.65,0.2,'THEMIS-D',/data,charsize=2, color=FSC_Color('red'), charthick=4, /NORMAL
 plots,0.63,0.21,psym=5, symsize=2, color=FSC_Color('red'),thick=5, /NORMAL
 plots,lon[index],lat[index],psym=5, symsize=2, color=FSC_Color('red'),thick=5 

 xyouts,0.65,0.15,'THEMIS-E',/data,charsize=2, color=FSC_Color('black'), charthick=4, /NORMAL
 plots,0.63,0.16,psym=6, symsize=2, color=FSC_Color('black'),thick=5, /NORMAL
 ;xyouts,lon1[index1]+0.5,lat1[index1]+0.1,'THEMIS-E',/data,charsize=2, color=FSC_Color('black'), charthick=4
 plots,lon1[index1],lat1[index1],psym=6, symsize=2, color=FSC_Color('black'), thick=5
 ;plots,lon2[index2],lat2[index2],psym=2, symsize=2, color=FSC_Color('red')
 
 ;  PFSIR Coordinate
 plat=65.12
 plong=-147.48
 
 xyouts,0.65,0.1,'PFISR',/data,charsize=2, color=FSC_Color('green'), charthick=4, /NORMAL
 plots,0.63,0.11,psym=7, symsize=2, color=FSC_Color('green'),thick=5, /NORMAL
 ;xyouts,plong+0.5,plat+0.1,'PFISR',/data,charsize=2, color=FSC_Color('green'), charthick=4
 plots,plong,plat,psym=7,symsize=2, color=FSC_Color('green'), thick=5

if keyword_set(disp) then begin  
device,/close
set_plot,'x'
endif

end