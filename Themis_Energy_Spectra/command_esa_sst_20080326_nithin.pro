timespan,'2008-3-26/8',5,/h
get_timespan,t

thm_load_esa,probe='d',trange=t,level=2
thm_load_sst2,probe='d',trange=t,level=2
thm_load_esa,probe='e',trange=t,level=2
thm_load_sst2,probe='e',trange=t,level=2

store_data,'thd_combined',data='thd_peef_en_eflux thd_psef_en_eflux'
store_data,'the_combined',data='the_peef_en_eflux the_psef_en_eflux'
zlim,'th?_combined',1e2,1e9,1
ylim,'th?_combined',1e2,1e6,1
tplot,'th?_combined'


