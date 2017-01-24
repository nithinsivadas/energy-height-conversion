timespan,'2008-3-26/8',5,/h
get_timespan,t

X = DOUBLE([t(0):t(1):1.0/16])

Y=4*sin(2*!PI*(X-X(0)))

store_data,'thd_eff_xy',data={x:X,y:Y,v:2}
spd_ui_wavelet,'thd_eff_xy','thd_eff_wlt_xy',[t(0),t(1)],maxpoints=2L^20
get_data, 'thd_eff_xy_wv_pow', data=xwvpow

p2= plot(xwvpow.V, xwvpow.Y(10000,*),$
  LAYOUT=[1,2,2],$
  YTITLE='Power ($nT^2/Hz$)',$
  XTITLE='Frequency ($\Hz$)',$
  TITLE = 'Wavelet Transform using spd_ui_wavelet.pro')
df = [xwvpow.V(0:93)-xwvpow.V(1:94), xwvpow.V(93)-xwvpow.V(94)]
total_power=TOTAL(xwvpow.Y(10000,*)*df)
t2 = TEXT(4, 15,'Net Power ='+[STRING(total_power)]+' $nT^2$',/DATA,FONT_SIZE=14, FONT_NAME='Helvetica')
p=plot(X(0:100)-X(0),Y(0:100),LAYOUT=[1,2,1], $
  YTITLE='Amplitude ($nT$)',$
  XTITLE='Time (s)',$
  TITLE='S(t) = 4*sin(2$\pi$t)', /CURRENT)
end