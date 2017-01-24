
;  Create ASI plots
;  By Nithin  Sivadas, 17th July 2016
;  If disp is not set, then it will display as a plot on the screen
;  if it is set to any value, then it will print as postscript
;  Example command: thm_pfisr_asi,show_time='2008-03-26/12:30:00',disp=1,filename='aurora_12_30_UT.ps'
pro thm_pfisr_asi, show_time=show_time, disp=disp, filename=filename


if keyword_set(disp) then begin 
set_plot,'PS'
device, /portrait, xsize=7, ysize=7, xoffset=0, yoffset=0, /inches, filename=filename
device, /color, bit=8
device, font_size=13
device, font_index=18
device, /Times
device, /encapsulated
endif


 timespan,'2008-3-26/8',5,/h
 
 ;show_time='2008-03-26/11:36:00'
 
 ; Downloading ASCI high resoution data
 thm_asi_create_mosaic,show_time,show=['FYKN','GAKO','WHIT'],central_lon=-150,central_lat=65,scale=2e7,projection='AzimuthalEquidistant', disp=disp
 
 ; Downloading THEMIS coordinates
 thm_load_state, probe =['d','e'], coord='gsm', suffix='_gsm'
 
 ; Calculating the foot prints
 ttrace2iono,'thd_state_pos_gsm', newname='thd_ifoot_geo', external_model='t89', par=2.0D,/km,in_coord='gsm',out_coord='geo'
 ;ttrace2iono,'thc_state_pos_gsm', newname='thc_ifoot_geo', external_model='t89', par=2.0D,/km,in_coord='gsm',out_coord='geo'
 ttrace2iono,'the_state_pos_gsm', newname='the_ifoot_geo', external_model='t89', par=2.0D,/km,in_coord='gsm',out_coord='geo'
 
 ; Converting it to lat and long coordinates
 get_data, 'thd_ifoot_geo', data=d
 lon =!radeg*atan(d.y[*,1],d.y[*,0])
 lat =!radeg*atan(d.y[*,2],sqrt(d.y[*,0]^2+d.y[*,1]^2))
 plots,lon,lat, color=FSC_Color('red'), linestyle=2, thick=4
 
 get_data, 'the_ifoot_geo', data=d1
 lon1 =!radeg*atan(d1.y[*,1],d1.y[*,0])
 lat1 =!radeg*atan(d1.y[*,2],sqrt(d1.y[*,0]^2+d1.y[*,1]^2))
 plots,lon1,lat1, color=FSC_Color('black'), linestyle=3, thick=4
 
 ; Creating a filled circle symbol
 A = FINDGEN(17)*(!PI*2/16.)
 USERSYM, COS(A), SIN(A), thick=2
 
 ; PFISR Beam Coordinates
 latPfisr = [64.9360126649696,  65.1825271637059,  65.2514751021475,  65.3796239705052,  65.5478236611008,  65.3340226640070,  65.7131074989978,  65.4951994661193,  65.9134020565136,  65.6592928489163,  65.8717925105221,  65.1299200000000,  65.7887045529653,  65.6122939769094,  65.4248505612056,  65.2900329725758,  65.7036735838334,  65.5095468945262,  65.6598939437277,  65.3415979860585,  65.4676475238508,  65.2978298743273,  65.1868612866175,  65.1344106429329,  65.0042515198086,  65.0532493517570];
 lonPfisr = [-147.690886717139,  -147.815897933303, -148.259938656474, -148.149227758499, -148.002549361192, -147.726474942448, -147.834144805863, -147.640478130197, -147.676333714786, -147.417046331727, -147.308603947527, -147.471040000000, -146.870063065926, -147.034198069652, -147.205996288591, -147.327940321313, -146.467759685847, -146.656241089619, -146.071588104171, -146.802595485301, -146.221479830276, -146.416899927615, -146.922399015966, -146.545064109959, -146.689441971007, -147.142593447949];
 plots,lonPfisr,latPfisr, psym=8,symsize=0.2, color=FSC_Color('green'), thick=0.5
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
 plots,0.63,0.11,psym=8,symsize=0.2, color=FSC_Color('green'), thick=0.5, /NORMAL
 ;xyouts,plong+0.5,plat+0.1,'PFISR',/data,charsize=2, color=FSC_Color('green'), charthick=4
 ;plots,plong,plat,psym=7,symsize=2, color=FSC_Color('green'), thick=5
 LAYOUT=[2,2,1]
 
if keyword_set(disp) then begin  
device,/close
set_plot,'x'
endif

end